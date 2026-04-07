# coding=utf-8
# pzw
# WES转座子插入突变分析流程
# 可分析Alu LINE-1 SVA HERV
# version 0.4 - 代码优化，添加VCF输出功能，移除VEP注释

"""
流程概述：
1. 提取WES数据中的SoftClip Reads，仅选择包含M和S且长度>=36bp的Reads
2. 统计断点信息，筛选支持reads数>=cutoff的断点
3. 将断点reads与TE参考库比对，获取转座子信息
4. 输出结果（TSV和VCF格式）

注：VEP注释应在下游工具中调用，不在本软件内进行
"""

import os
import sys
import subprocess
import pandas as pd
import pysam
import argparse
import configparser
from datetime import datetime

# 常量定义
MIN_SOFTCLIP_LENGTH = 36
MIN_CIGAR_OPERATIONS = 2
VALID_CHROMOSOMES = (
    [f"chr{i}" for i in range(1, 23)] +
    [str(i) for i in range(1, 23)] +
    ["chrX", "chrY", "X", "Y"]
)

# VCF头部模板
VCF_HEADER_TEMPLATE = '''##fileformat=VCFv4.3
##fileDate={date}
##source=TIEA-WES_v0.4
##reference={reference}
##ALT=<ID=INS:ME:ALU,Description="Alu insertion">
##ALT=<ID=INS:ME:LINE1,Description="LINE-1 insertion">
##ALT=<ID=INS:ME:SVA,Description="SVA insertion">
##ALT=<ID=INS:ME:HERV,Description="HERV insertion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info: NAME,START,END,POLARITY">
##INFO=<ID=MINAME,Number=1,Type=String,Description="Mobile element name">
##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of supporting reads">
##INFO=<ID=MAPPING,Number=1,Type=String,Description="TE mapping details">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
'''


def get_abs_path():
    """获取脚本所在目录的绝对路径"""
    return os.path.dirname(os.path.abspath(__file__))


def process_read(read, min_len=MIN_SOFTCLIP_LENGTH):
    """
    检查read是否符合要求：
    - 不是重复或补充read
    - CIGAR只包含M和S两种操作
    - M和S长度均>=min_len
    """
    if read.is_duplicate or read.is_supplementary:
        return None

    cigar = read.cigartuples
    if cigar is None:
        return None

    if len(cigar) != MIN_CIGAR_OPERATIONS:
        return None

    for op in cigar:
        if op[0] not in [0, 4]:  # 只允许M(0)和S(4)
            return None
        if op[1] < min_len:
            return None

    return read


def filter_softclip_reads(input_bam, output_bam, min_len=MIN_SOFTCLIP_LENGTH):
    """筛选符合要求的SoftClip Reads"""
    with pysam.AlignmentFile(input_bam, "rb") as bamfile:
        with pysam.AlignmentFile(output_bam, "wb", template=bamfile) as outbam:
            for read in bamfile:
                result = process_read(read, min_len)
                if result is not None:
                    outbam.write(result)


def get_breakpoints(read):
    """从CIGAR字符串中提取断点位置"""
    cigar = read.cigartuples
    chrom = read.reference_name
    pos = read.reference_start + 1  # 转为1-based

    breakpoints = []
    current_pos = pos

    for op, length in cigar:
        if op == 0:  # M操作
            current_pos += length
        elif op == 4:  # S操作
            if current_pos > pos:  # 软剪切在匹配之后
                breakpoints.append((chrom, current_pos, read.query_name))

    return breakpoints


def count_breakpoints(bam_file):
    """统计断点信息"""
    breakpoint_counts = {}
    with pysam.AlignmentFile(bam_file, "rb") as bamfile:
        for read in bamfile:
            breakpoints = get_breakpoints(read)
            for chrom, pos, read_name in breakpoints:
                key = (chrom, pos)
                if key not in breakpoint_counts:
                    breakpoint_counts[key] = {"reads": 0, "reads_name": set()}
                breakpoint_counts[key]["reads"] += 1
                breakpoint_counts[key]["reads_name"].add(read_name)

    for key in breakpoint_counts:
        breakpoint_counts[key]["reads_name"] = list(breakpoint_counts[key]["reads_name"])

    return breakpoint_counts


def get_breakpoint_df(input_bam, cutoff=10):
    """整理断点信息为DataFrame并筛选"""
    breakpoint_counts = count_breakpoints(input_bam)

    if not breakpoint_counts:
        return pd.DataFrame(columns=["chrom", "pos", "reads", "reads_name"])

    df = pd.DataFrame(list(breakpoint_counts.values()),
                      index=pd.MultiIndex.from_tuples(breakpoint_counts.keys(),
                                                       names=['chrom', 'pos']))
    df.index = df.index.set_levels(df.index.levels[0].astype(str), level=0)
    df.reset_index(inplace=True)
    df = df[['chrom', 'pos', 'reads', 'reads_name']]

    df = df[df['reads'] >= cutoff]
    df = df[df['chrom'].isin(VALID_CHROMOSOMES)]
    df.reset_index(drop=True, inplace=True)

    return df.copy()


def get_soft_clipped_parts(read):
    """提取软剪切部分的序列和质量值"""
    cigar = read.cigartuples
    sequence = read.query_sequence
    qualities = read.query_qualities
    soft_clipped_parts = []

    start_clip = end_clip = None
    current_pos = 0

    for op, length in cigar:
        if op == 4:  # Soft Clip
            if start_clip is None:
                start_clip = current_pos
            end_clip = current_pos + length
            soft_clipped_parts.append((sequence[start_clip:end_clip],
                                       qualities[start_clip:end_clip]))
            start_clip = end_clip = None
        elif op in [0, 1, 2, 7, 8]:
            current_pos += length

    return soft_clipped_parts


def split_breakpoint_2_fastq(df, input_bam, tmp_dir):
    """将断点reads提取为fastq文件（优化：单次遍历BAM）"""
    read_to_breakpoints = {}
    breakpoint_files = {}

    for idx, row in df.iterrows():
        chrom = row['chrom']
        pos = row['pos']
        reads_name = row['reads_name']

        output_path = os.path.join(tmp_dir, f"{chrom}_{pos}")
        os.makedirs(output_path, exist_ok=True)
        fastq_path = os.path.join(output_path, f"{chrom}_{pos}.fastq")
        breakpoint_files[(chrom, pos)] = open(fastq_path, 'w')

        for read_name in reads_name:
            if read_name not in read_to_breakpoints:
                read_to_breakpoints[read_name] = []
            read_to_breakpoints[read_name].append((chrom, pos))

    with pysam.AlignmentFile(input_bam, "rb") as bamfile:
        for read in bamfile:
            if read.query_name in read_to_breakpoints:
                soft_clipped_parts = get_soft_clipped_parts(read)
                for seq, quals in soft_clipped_parts:
                    for chrom, pos in read_to_breakpoints[read.query_name]:
                        fh = breakpoint_files[(chrom, pos)]
                        fh.write(f"@{read.query_name}\n{seq}\n+\n")
                        fh.write("".join(chr(q + 33) for q in quals) + "\n")

    for fh in breakpoint_files.values():
        fh.close()


def run_bwa_alignment(bwa, te_reference, fastq, sam, threads=2):
    """执行BWA比对"""
    try:
        result = subprocess.run(
            [bwa, "mem", "-t", str(threads), te_reference, fastq],
            capture_output=True, text=True, check=True
        )
        with open(sam, 'w') as f:
            f.write(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"BWA alignment failed: {e.stderr}")
        with open(sam, 'w') as f:
            f.write("@HD\tVN:1.0\n")


def te_mapping(row, te_reference, tmp_dir, bwa):
    """比对TE库并统计结果"""
    row["te_dict"] = {"SVA": 0, "L1": 0, "Alu": 0, "HERV": 0}

    chrom = row['chrom']
    pos = row['pos']
    fastq = os.path.join(tmp_dir, f"{chrom}_{pos}/{chrom}_{pos}.fastq")
    sam = os.path.join(tmp_dir, f"{chrom}_{pos}/{chrom}_{pos}.sam")

    run_bwa_alignment(bwa, te_reference, fastq, sam)

    with pysam.AlignmentFile(sam, "r") as samfile:
        for read in samfile:
            gene = read.reference_name
            if gene is None:
                continue
            row["mapping_dict"].setdefault(gene, 0)
            row["mapping_dict"][gene] += 1

            for te_type in ["SVA", "L1", "Alu", "HERV"]:
                if te_type in gene and gene.startswith(te_type):
                    row["te_dict"][te_type] += 1

    return row


def mapping_filter(input_df):
    """过滤无TE比对的断点"""
    if input_df.empty:
        return pd.DataFrame(columns=["chrom", "pos", "reads", "mapping_dict", "te_dict"])
    output_df = input_df[input_df['mapping_dict'].apply(lambda x: len(x) > 0)]
    output_df = output_df[["chrom", "pos", "reads", "mapping_dict", "te_dict"]]
    return output_df


# ========== VCF输出功能 ==========

def get_primary_te_type(te_dict):
    """从te_dict获取主要TE类型"""
    if not te_dict:
        return "UNKNOWN"
    max_count = max(te_dict.values())
    if max_count == 0:
        return "UNKNOWN"
    for te_type, count in te_dict.items():
        if count == max_count:
            return te_type
    return "UNKNOWN"


def generate_vcf_header(reference):
    """生成VCF头部"""
    return VCF_HEADER_TEMPLATE.format(
        date=datetime.now().strftime("%Y%m%d"),
        reference=reference
    )


def generate_vcf_record(row, sample_id):
    """生成单条VCF记录"""
    chrom = row['chrom']
    pos = row['pos']
    reads = row.get('reads', 0)

    te_type = get_primary_te_type(row.get('te_dict', {}))

    te_to_alt = {
        'Alu': '<INS:ME:ALU>',
        'L1': '<INS:ME:LINE1>',
        'SVA': '<INS:ME:SVA>',
        'HERV': '<INS:ME:HERV>'
    }
    alt = te_to_alt.get(te_type, '<INS:ME:UNKNOWN>')

    info_parts = [
        f"SVTYPE=INS",
        f"SVLEN=.",
        f"MINAME={te_type}",
        f"SUPPORT={reads}"
    ]

    mapping_dict = row.get('mapping_dict', {})
    if mapping_dict:
        mapping_str = ",".join([f"{k}:{v}" for k, v in mapping_dict.items()])
        info_parts.append(f"MAPPING={mapping_str}")

    record_id = f"{sample_id}_MEI_{chrom}_{pos}"
    ref = "N"
    qual = "."
    filt = "PASS" if reads >= 10 else "LowSupport"
    info = ";".join(info_parts)

    return f"{chrom}\t{pos}\t{record_id}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}"


def write_vcf_output(df, output_path, sample_id, reference):
    """写入VCF文件"""
    with open(output_path, 'w', encoding='utf-8') as vcf_file:
        header = generate_vcf_header(reference)
        vcf_file.write(header)

        df_sorted = df.sort_values(['chrom', 'pos'])

        for idx, row in df_sorted.iterrows():
            record = generate_vcf_record(row, sample_id)
            vcf_file.write(record + "\n")


# ========== 主流程 ==========

def te_pipe(prefix, input_bam, output_dir, sup_read, min_len,
            reference, rm_tmp, config_file, output_vcf=False):
    """主分析流程"""

    # 加载配置
    if config_file == "AUTO_DETECT":
        config_file = os.path.join(get_abs_path(), "config.ini")
    config = configparser.ConfigParser()
    config.read(config_file)

    bwa = config["software"]["BWA"]
    samtools = config["software"]["SAMTOOLS"]
    te_genome = config["database"]["TE_GENOME"]

    # 构建输出路径
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    tmp_dir = os.path.join(output_dir, "tmp_" + prefix)
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    output_bam = os.path.join(output_dir, prefix + ".softclip.bam")
    output_tsv = os.path.join(output_dir, prefix + ".te.result.txt")

    # Step 1: 筛选SoftClip Reads
    print("Step 1: Filtering SoftClip reads...")
    filter_softclip_reads(input_bam, output_bam, min_len)
    subprocess.run([samtools, "index", output_bam], check=True, capture_output=True)

    # Step 2: 统计断点
    print("Step 2: Counting breakpoints...")
    df = get_breakpoint_df(output_bam, sup_read)

    if df.empty:
        print("No breakpoints found.")
        empty_df = pd.DataFrame(columns=["chrom", "pos", "reads", "mapping_dict", "te_dict"])
        empty_df.to_csv(output_tsv, header=True, index=False, sep="\t")
        if output_vcf:
            vcf_path = os.path.join(output_dir, prefix + ".te.result.vcf")
            write_vcf_output(empty_df, vcf_path, prefix, reference)
        return

    # Step 3: 生成fastq
    print("Step 3: Generating FASTQ files...")
    split_breakpoint_2_fastq(df, output_bam, tmp_dir)

    # Step 4: TE比对
    print("Step 4: Mapping to TE library...")
    df['mapping_dict'] = [{} for _ in range(len(df))]
    df = df.apply(lambda x: te_mapping(x, te_genome, tmp_dir, bwa), axis=1)

    # Step 5: 过滤结果
    print("Step 5: Filtering results...")
    df = mapping_filter(df)

    if df.empty:
        print("No TE insertions detected.")
        empty_df = pd.DataFrame(columns=["chrom", "pos", "reads", "mapping_dict", "te_dict"])
        empty_df.to_csv(output_tsv, header=True, index=False, sep="\t")
        if output_vcf:
            vcf_path = os.path.join(output_dir, prefix + ".te.result.vcf")
            write_vcf_output(empty_df, vcf_path, prefix, reference)
        return

    # 输出TSV
    df.to_csv(output_tsv, header=True, index=False, sep="\t")
    print(f"TSV output: {output_tsv}")

    # 输出VCF（可选）
    if output_vcf:
        print("Step 6: Writing VCF output...")
        vcf_path = os.path.join(output_dir, prefix + ".te.result.vcf")
        write_vcf_output(df, vcf_path, prefix, reference)
        print(f"VCF output: {vcf_path}")

    # 清理临时文件
    if rm_tmp in ["False", "false"]:
        import shutil
        shutil.rmtree(tmp_dir)


def main():
    parser = argparse.ArgumentParser(
        description="TIEA-WES: Transposon Insertion Event Analyzer for WES data",
        prog="TIEA-WES.py",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version", version="Version 0.4")
    parser.add_argument("-p", "--prefix", type=str, required=True,
                        help="Sample prefix/identifier")
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Input BAM file path")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output directory path")
    parser.add_argument("-f", "--ref", type=str, required=True,
                        help="Reference FASTA file path (for VCF header)")
    parser.add_argument("-s", "--support_reads", type=int, default=10,
                        help="Minimum breakpoint support reads cutoff (default: 10)")
    parser.add_argument("-l", "--min_length", type=int, default=36,
                        help="Minimum softclip length (default: 36)")
    parser.add_argument("-t", "--keep_tmp", type=str, default="False",
                        help="Keep temporary files (default: False)")
    parser.add_argument("-c", "--config", type=str, default="AUTO_DETECT",
                        help="Config file path (default: auto-detect)")
    parser.add_argument("--vcf", action="store_true",
                        help="Output results in VCF format")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    te_pipe(
        prefix=args.prefix,
        input_bam=args.input,
        output_dir=args.output,
        sup_read=args.support_reads,
        min_len=args.min_length,
        reference=args.ref,
        rm_tmp=args.keep_tmp,
        config_file=args.config,
        output_vcf=args.vcf
    )


if __name__ == "__main__":
    main()