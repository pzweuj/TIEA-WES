# coding=utf-8
# pzw
# WES转座子插入突变分析流程
# 可分析Alu LINE-1 SVA HERV
# version 2.0 - 全面重构：合并流程、管道解析、双向softclip、MAPQ过滤、断点聚类

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
from multiprocessing import Pool, cpu_count

try:
    import pyfaidx
    HAS_PYFAIDX = True
except ImportError:
    HAS_PYFAIDX = False

# 常量定义
MIN_SOFTCLIP_LENGTH = 36
MIN_CIGAR_OPERATIONS = 2
MIN_MAPQ = 20
CLUSTER_WINDOW = 10
VALID_CHROMOSOMES = (
    [f"chr{i}" for i in range(1, 23)] +
    [str(i) for i in range(1, 23)] +
    ["chrX", "chrY", "X", "Y", "chrM", "MT"]
)

# VCF头部模板
VCF_HEADER_TEMPLATE = '''##fileformat=VCFv4.3
##fileDate={date}
##source=TIEA-WES_v2.0.0
##reference={reference}
##ALT=<ID=INS:ME:ALU,Description="Alu insertion">
##ALT=<ID=INS:ME:LINE1,Description="LINE-1 insertion">
##ALT=<ID=INS:ME:SVA,Description="SVA insertion">
##ALT=<ID=INS:ME:HERV,Description="HERV insertion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=MINAME,Number=1,Type=String,Description="Mobile element name (TE type)">
##INFO=<ID=MEFAMILY,Number=1,Type=String,Description="Best matching TE family (e.g., AluY, L1HS)">
##INFO=<ID=CONFIDENCE,Number=1,Type=String,Description="Classification confidence: High, Medium, Low">
##INFO=<ID=DIR,Number=1,Type=String,Description="Insertion direction: 5_prime or 3_prime">
##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of supporting reads">
##INFO=<ID=AVGSCLEN,Number=1,Type=Float,Description="Average softclip length">
##INFO=<ID=CLUSTERSIZE,Number=1,Type=Integer,Description="Number of breakpoints merged in cluster">
##INFO=<ID=CLUSTERRANGE,Number=1,Type=String,Description="Position range of clustered breakpoints">
##INFO=<ID=MAPPING,Number=.,Type=String,Description="Detailed TE mapping: family:count pairs">
##INFO=<ID=TECOUNT,Number=.,Type=String,Description="TE type counts: Alu:N,L1:N,SVA:N,HERV:N">
##FILTER=<ID=LowSupport,Description="Less than 10 supporting reads">
##FILTER=<ID=LowConfidence,Description="TE type classification confidence is low">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
'''


def get_abs_path():
    """获取脚本所在目录的绝对路径"""
    return os.path.dirname(os.path.abspath(__file__))


def get_reference_from_fasta(fasta_path):
    """
    从FASTA文件路径提取参考基因组名称

    例如：
    /path/to/hg19.fa -> hg19
    /path/to/GRCh38.primary_assembly.genome.fa -> GRCh38
    """
    if not fasta_path:
        return None

    basename = os.path.basename(fasta_path)
    name = basename.split('.')[0]  # 去除扩展名

    # 常见命名规范化
    name_lower = name.lower()
    if 'hg19' in name_lower or 'grch37' in name_lower:
        if 'hg19' in name_lower:
            return 'hg19'
        return 'GRCh37'
    elif 'hg38' in name_lower or 'grch38' in name_lower:
        if 'hg38' in name_lower:
            return 'hg38'
        return 'GRCh38'

    return name


def get_reference_from_bam(bam_file):
    """
    从BAM header提取参考基因组信息

    优先级：
    1. @SQ 行的 AS 字段（assembly name，如 hg38, GRCh37）
    2. @SQ 行的 UR 字段（URI/path）
    3. 根据特征序列名/长度推断

    返回: reference name string
    """
    with pysam.AlignmentFile(bam_file, "rb") as bamfile:
        header = bamfile.header

        # 方法1: 从 @SQ 的 AS 字段获取（最可靠）
        for sq_line in header.get('SQ', []):
            if 'AS' in sq_line:
                return sq_line['AS']

        # 方法2: 从 @SQ 的 UR 字段获取
        for sq_line in header.get('SQ', []):
            if 'UR' in sq_line:
                ur = sq_line['UR']
                if ur:
                    return ur

        # 方法3: 根据特征序列名推断
        sn_list = [sq.get('SN', '') for sq in header.get('SQ', [])]

        # 检查 hs37d5 (GRCh37+decoy)
        if 'hs37d5' in sn_list:
            return 'GRCh37'

        # 检查 hs38d1 (GRCh38+decoy)
        if 'hs38d1' in sn_list or 'hs38DH' in sn_list:
            return 'GRCh38'

        # 检查染色体命名风格
        has_chr_prefix = any(sn.startswith('chr') for sn in sn_list if sn)
        has_MT = 'MT' in sn_list
        has_chrM = 'chrM' in sn_list

        # 根据第一个常染色体的长度判断
        # GRCh37: chr1 = 249250621
        # GRCh38: chr1 = 248956422
        chr1_len = None
        for sq in header.get('SQ', []):
            sn = sq.get('SN', '')
            if sn in ['1', 'chr1']:
                chr1_len = sq.get('LN', 0)
                break

        if chr1_len:
            # GRCh37 chr1 长度约 249M
            # GRCh38 chr1 长度约 248M
            if chr1_len > 249000000:
                if has_chr_prefix:
                    return 'hg19'
                else:
                    return 'GRCh37'
            elif chr1_len < 249000000:
                if has_chr_prefix:
                    return 'hg38'
                else:
                    return 'GRCh38'

        return 'unknown'


def get_reference_info(bam_file):
    """
    获取参考基因组详细信息用于提示用户

    返回: (reference, has_chr_prefix, suggestion)
    """
    with pysam.AlignmentFile(bam_file, "rb") as bamfile:
        header = bamfile.header

        # 检查是否有 chr 前缀
        sn_list = [sq.get('SN', '') for sq in header.get('SQ', [])]
        has_chr_prefix = any(sn.startswith('chr') for sn in sn_list if sn)

        # 获取参考基因组名称
        reference = get_reference_from_bam(bam_file)

        # 生成建议
        if reference != 'unknown':
            if has_chr_prefix:
                suggestion = f"Detected '{reference}' (with chr prefix)"
            else:
                suggestion = f"Detected '{reference}' (no chr prefix)"
        else:
            if has_chr_prefix:
                suggestion = "Reference has 'chr' prefix. Could be hg19 (GRCh37) or hg38 (GRCh38)."
            else:
                suggestion = "Reference has no 'chr' prefix. Could be GRCh37 or GRCh38."

        return reference, has_chr_prefix, suggestion


def process_read_check(read, min_len=MIN_SOFTCLIP_LENGTH, min_mapq=MIN_MAPQ):
    """
    检查read是否符合要求：
    - MAPQ >= min_mapq
    - 不是重复或补充read
    - CIGAR只包含M和S两种操作
    - M和S长度均>=min_len

    返回: (is_valid, has_5prime_softclip, has_3prime_softclip)
    """
    if read.is_duplicate or read.is_supplementary:
        return False, False, False

    if read.mapping_quality < min_mapq:
        return False, False, False

    cigar = read.cigartuples
    if cigar is None:
        return False, False, False

    if len(cigar) != MIN_CIGAR_OPERATIONS:
        return False, False, False

    # 检查只包含M和S
    for op in cigar:
        if op[0] not in [0, 4]:  # 只允许M(0)和S(4)
            return False, False, False
        if op[1] < min_len:
            return False, False, False

    # 判断softclip方向
    has_5prime = cigar[0][0] == 4  # 开头是S（5' softclip）
    has_3prime = cigar[-1][0] == 4  # 末尾是S（3' softclip）

    return True, has_5prime, has_3prime


def get_breakpoints_birectional(read, has_5prime, has_3prime):
    """
    从CIGAR字符串中提取断点位置（支持双向）

    5' softclip (S+M): 断点在匹配起始位置
    3' softclip (M+S): 断点在匹配结束位置

    返回: [(chrom, pos, direction, read_name)]
    """
    cigar = read.cigartuples
    chrom = read.reference_name
    if chrom is None:
        return []

    ref_start = read.reference_start + 1  # 转为1-based
    breakpoints = []
    current_pos = ref_start

    # 计算5'断点（如果有）
    if has_5prime:
        # S+M模式：断点在M操作开始位置
        softclip_len = cigar[0][1]
        # 5' softclip 表示插入发生在参考基因组该位置之前
        breakpoints.append((chrom, ref_start, "5_prime", read.query_name, softclip_len))

    # 遍历CIGAR计算3'断点
    if has_3prime:
        for op, length in cigar:
            if op == 0:  # M操作
                current_pos += length
        # 3' softclip 表示插入发生在参考基因组该位置之后
        softclip_len = cigar[-1][1]
        breakpoints.append((chrom, current_pos, "3_prime", read.query_name, softclip_len))

    return breakpoints


def get_softclip_sequence(read, has_5prime, has_3prime):
    """
    提取softclip部分的序列和质量值（支持双向）

    返回: [(sequence, qualities, direction)]
    """
    cigar = read.cigartuples
    sequence = read.query_sequence
    qualities = read.query_qualities

    softclip_parts = []
    current_pos = 0

    for op, length in cigar:
        if op == 4:  # Soft Clip
            # 判断方向
            if current_pos == 0 and has_5prime:
                direction = "5_prime"
            else:
                direction = "3_prime"

            softclip_parts.append((
                sequence[current_pos:current_pos + length],
                qualities[current_pos:current_pos + length],
                direction
            ))
            current_pos += length
        elif op in [0, 1, 2, 7, 8]:
            current_pos += length

    return softclip_parts


def process_bam_once(input_bam, sup_read, min_len, min_mapq, tmp_dir):
    """
    一次遍历完成：筛选 + 统计断点 + 写FASTQ

    返回: (breakpoint_df, combined_fastq_path)

    替代原来的 Step 1 + Step 2 + Step 3
    """
    breakpoint_counts = {}  # {(chrom, pos, direction): {"reads": 0, "reads_name": set(), "softclip_len": []}}

    combined_fastq = os.path.join(tmp_dir, "combined_breakpoints.fastq")

    total_reads = 0
    passed_mapq = 0
    passed_cigar = 0
    passed_softclip = 0
    fastq_records = 0

    with pysam.AlignmentFile(input_bam, "rb") as bamfile:
        with open(combined_fastq, 'w') as fq_out:
            for read in bamfile:
                total_reads += 1

                # 检查 MAPQ
                if read.is_duplicate or read.is_supplementary:
                    continue

                if read.mapping_quality < min_mapq:
                    continue
                passed_mapq += 1

                # 检查 CIGAR
                cigar = read.cigartuples
                if cigar is None:
                    continue

                if len(cigar) != MIN_CIGAR_OPERATIONS:
                    continue

                # 检查只包含 M 和 S
                valid_cigar = True
                for op in cigar:
                    if op[0] not in [0, 4]:
                        valid_cigar = False
                        break
                    if op[1] < min_len:
                        valid_cigar = False
                        break

                if not valid_cigar:
                    continue
                passed_cigar += 1

                # 判断 softclip 方向
                has_5prime = cigar[0][0] == 4
                has_3prime = cigar[-1][0] == 4

                if not has_5prime and not has_3prime:
                    continue
                passed_softclip += 1

                # 提取断点
                breakpoints = get_breakpoints_birectional(read, has_5prime, has_3prime)

                if not breakpoints:
                    continue

                # 提取softclip序列
                softclip_parts = get_softclip_sequence(read, has_5prime, has_3prime)

                # 处理每个断点
                for chrom, pos, direction, read_name, softclip_len in breakpoints:
                    key = (chrom, pos, direction)

                    # 统计断点
                    if key not in breakpoint_counts:
                        breakpoint_counts[key] = {
                            "reads": 0,
                            "reads_name": set(),
                            "softclip_lens": []
                        }
                    breakpoint_counts[key]["reads"] += 1
                    breakpoint_counts[key]["reads_name"].add(read_name)
                    breakpoint_counts[key]["softclip_lens"].append(softclip_len)

                    # 写入FASTQ - 找到对应方向的softclip序列
                    for seq, quals, sc_dir in softclip_parts:
                        if sc_dir == direction:
                            modified_name = f"{read_name}@{chrom}_{pos}_{direction}"
                            fq_out.write(f"@{modified_name}\n{seq}\n+\n")
                            fq_out.write("".join(chr(q + 33) for q in quals) + "\n")
                            fastq_records += 1
                            break

    print(f"  Total reads scanned: {total_reads}")
    print(f"  Reads passed MAPQ>={min_mapq}: {passed_mapq}")
    print(f"  Reads passed CIGAR check: {passed_cigar}")
    print(f"  Reads with valid softclip: {passed_softclip}")
    print(f"  FASTQ records written: {fastq_records}")

    # 转换为DataFrame并筛选
    if not breakpoint_counts:
        return pd.DataFrame(columns=["chrom", "pos", "direction", "reads", "reads_name", "avg_softclip_len"]), combined_fastq

    df_data = []
    for key, data in breakpoint_counts.items():
        chrom, pos, direction = key
        avg_len = sum(data["softclip_lens"]) / len(data["softclip_lens"]) if data["softclip_lens"] else 0
        df_data.append({
            "chrom": chrom,
            "pos": pos,
            "direction": direction,
            "reads": data["reads"],
            "reads_name": list(data["reads_name"]),
            "avg_softclip_len": avg_len
        })

    df = pd.DataFrame(df_data)

    # 筛选
    df = df[df['reads'] >= sup_read]
    df = df[df['chrom'].isin(VALID_CHROMOSOMES)]
    df.reset_index(drop=True, inplace=True)

    return df.copy(), combined_fastq


def cluster_breakpoints(df, window=CLUSTER_WINDOW):
    """
    断点聚类合并

    合并同一染色体上相近位置（±window bp）的断点
    保留支持reads最多的作为代表

    返回: 聚合后的DataFrame
    """
    if df.empty:
        return df

    # 按染色体和位置排序
    df_sorted = df.sort_values(['chrom', 'pos', 'direction']).copy()

    clustered = []
    current_cluster = []

    for idx, row in df_sorted.iterrows():
        if not current_cluster:
            current_cluster = [row]
            continue

        # 检查是否在同一聚类
        last_row = current_cluster[-1]
        same_chrom = row['chrom'] == last_row['chrom']
        same_direction = row['direction'] == last_row['direction']
        close_pos = abs(row['pos'] - last_row['pos']) <= window

        if same_chrom and same_direction and close_pos:
            current_cluster.append(row)
        else:
            # 处理当前聚类
            if current_cluster:
                clustered_row = _merge_cluster(current_cluster)
                clustered.append(clustered_row)
            current_cluster = [row]

    # 处理最后一个聚类
    if current_cluster:
        clustered_row = _merge_cluster(current_cluster)
        clustered.append(clustered_row)

    return pd.DataFrame(clustered).reset_index(drop=True)


def _merge_cluster(cluster_rows):
    """
    合并一个聚类中的断点

    选择支持reads最多的作为代表位置
    """
    # 选择reads最多的作为代表
    best_row = max(cluster_rows, key=lambda x: x['reads'])

    # 合并所有reads_name
    all_reads = set()
    total_reads = 0
    all_softclip_lens = []

    for row in cluster_rows:
        all_reads.update(row['reads_name'])
        total_reads += row['reads']
        all_softclip_lens.append(row['avg_softclip_len'])

    # 创建合并后的行
    merged = {
        "chrom": best_row['chrom'],
        "pos": best_row['pos'],  # 使用代表位置
        "direction": best_row['direction'],
        "reads": len(all_reads),  # 使用唯一reads数
        "reads_name": list(all_reads),
        "avg_softclip_len": sum(all_softclip_lens) / len(all_softclip_lens) if all_softclip_lens else 0,
        "cluster_size": len(cluster_rows),  # 记录聚类大小
        "cluster_range": f"{min(r['pos'] for r in cluster_rows)}-{max(r['pos'] for r in cluster_rows)}"
    }

    return merged


def parse_bwa_stdout_piped(stdout_text):
    """
    管道解析BWA输出，不写入中间SAM文件

    返回: breakpoint_stats dict
    """
    breakpoint_stats = {}

    for line in stdout_text.splitlines():
        if line.startswith('@'):
            continue

        fields = line.strip().split('\t')
        if len(fields) < 3:
            continue

        read_name = fields[0]
        flag = int(fields[1])
        ref_name = fields[2]

        # 跳过未比对的reads
        if ref_name == '*' or flag & 4:
            continue

        # 从read name提取断点信息
        if read_name is None or "@" not in read_name:
            continue

        parts = read_name.split("@")
        if len(parts) != 2:
            continue

        breakpoint_key = parts[1]  # chrom_pos_direction

        if breakpoint_key not in breakpoint_stats:
            breakpoint_stats[breakpoint_key] = {
                "mapping_dict": {},
                "te_dict": {"SVA": 0, "L1": 0, "Alu": 0, "HERV": 0}
            }

        # 统计比对结果
        gene = ref_name

        breakpoint_stats[breakpoint_key]["mapping_dict"].setdefault(gene, 0)
        breakpoint_stats[breakpoint_key]["mapping_dict"][gene] += 1

        for te_type in ["SVA", "L1", "Alu", "HERV"]:
            if te_type in gene and gene.startswith(te_type):
                breakpoint_stats[breakpoint_key]["te_dict"][te_type] += 1

    return breakpoint_stats


# ========== 旧函数（保留备用） ==========

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
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True
        )
        with open(sam, 'w') as f:
            f.write(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"BWA alignment failed: {e.stderr}")
        with open(sam, 'w') as f:
            f.write("@HD\tVN:1.0\n")


def split_breakpoint_2_fastq_batch(df, input_bam, tmp_dir):
    """
    将所有断点reads合并为一个fastq文件（批量比对优化）
    read name格式: {original_name}@{chrom}_{pos} 用于回溯断点
    """
    read_to_breakpoints = {}

    for idx, row in df.iterrows():
        chrom = row['chrom']
        pos = row['pos']
        reads_name = row['reads_name']

        for read_name in reads_name:
            if read_name not in read_to_breakpoints:
                read_to_breakpoints[read_name] = []
            read_to_breakpoints[read_name].append((chrom, pos))

    combined_fastq = os.path.join(tmp_dir, "combined_breakpoints.fastq")

    with pysam.AlignmentFile(input_bam, "rb") as bamfile:
        with open(combined_fastq, 'w') as fq_out:
            for read in bamfile:
                if read.query_name in read_to_breakpoints:
                    soft_clipped_parts = get_soft_clipped_parts(read)
                    for seq, quals in soft_clipped_parts:
                        for chrom, pos in read_to_breakpoints[read.query_name]:
                            # 添加断点标记后缀
                            modified_name = f"{read.query_name}@{chrom}_{pos}"
                            fq_out.write(f"@{modified_name}\n{seq}\n+\n")
                            fq_out.write("".join(chr(q + 33) for q in quals) + "\n")

    return combined_fastq


def batch_te_mapping(combined_fastq, te_reference, tmp_dir, bwa, threads=4):
    """
    批量比对所有断点到TE库，管道解析结果（不写中间SAM）

    返回: breakpoint_stats dict
    """
    # 执行批量BWA比对（只启动一次，管道解析）
    try:
        result = subprocess.run(
            [bwa, "mem", "-t", str(threads), te_reference, combined_fastq],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True
        )
        # 直接解析stdout，不写SAM文件
        breakpoint_stats = parse_bwa_stdout_piped(result.stdout)

    except subprocess.CalledProcessError as e:
        print(f"BWA alignment failed: {e.stderr}")
        breakpoint_stats = {}

    return breakpoint_stats


# ========== 保留原单断点比对逻辑（备用） ==========

# 并行处理的全局变量
_parallel_te_reference = None
_parallel_tmp_dir = None
_parallel_bwa = None


def _te_mapping_worker(row_dict):
    """并行处理 worker 函数"""
    global _parallel_te_reference, _parallel_tmp_dir, _parallel_bwa

    result = {
        "chrom": row_dict["chrom"],
        "pos": row_dict["pos"],
        "reads": row_dict["reads"],
        "mapping_dict": {},
        "te_dict": {"SVA": 0, "L1": 0, "Alu": 0, "HERV": 0}
    }

    chrom = row_dict['chrom']
    pos = row_dict['pos']
    fastq = os.path.join(_parallel_tmp_dir, f"{chrom}_{pos}/{chrom}_{pos}.fastq")
    sam = os.path.join(_parallel_tmp_dir, f"{chrom}_{pos}/{chrom}_{pos}.sam")

    try:
        subprocess.run(
            [_parallel_bwa, "mem", "-t", "1", _parallel_te_reference, fastq],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True
        )
        with open(sam, 'w') as f:
            subprocess.run(
                [_parallel_bwa, "mem", "-t", "1", _parallel_te_reference, fastq],
                stdout=f, stderr=subprocess.PIPE, universal_newlines=True, check=True
            )
    except subprocess.CalledProcessError:
        with open(sam, 'w') as f:
            f.write("@HD\tVN:1.0\n")

    with pysam.AlignmentFile(sam, "r") as samfile:
        for read in samfile:
            gene = read.reference_name
            if gene is None:
                continue
            result["mapping_dict"].setdefault(gene, 0)
            result["mapping_dict"][gene] += 1

            for te_type in ["SVA", "L1", "Alu", "HERV"]:
                if te_type in gene and gene.startswith(te_type):
                    result["te_dict"][te_type] += 1

    return result


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
        return pd.DataFrame(columns=["chrom", "pos", "direction", "reads", "mapping_dict", "te_dict"])
    output_df = input_df[input_df['mapping_dict'].apply(lambda x: len(x) > 0)]
    # 保留所有可用列
    cols_to_keep = ["chrom", "pos", "direction", "reads", "mapping_dict", "te_dict"]
    if "avg_softclip_len" in output_df.columns:
        cols_to_keep.append("avg_softclip_len")
    if "cluster_size" in output_df.columns:
        cols_to_keep.append("cluster_size")
    if "cluster_range" in output_df.columns:
        cols_to_keep.append("cluster_range")
    output_df = output_df[cols_to_keep]
    return output_df


# ========== VCF输出功能 ==========

def get_primary_te_type_with_confidence(mapping_dict, te_dict):
    """
    改进的TE类型决策：

    算法：
    1. 计算每个大类内的家族一致性（主家族占比）
    2. 综合得分 = 计数 × 一致性
    3. 比较得分，确定置信度

    返回：(te_type, confidence, best_family)

    示例：
    mapping_dict = {"AluY": 8, "AluSx": 2, "L1HS": 10}
    te_dict = {"Alu": 10, "L1": 10, "SVA": 0, "HERV": 0}

    Alu: 计数10, 一致性=8/10=0.8, 得分=8
    L1: 计数10, 一致性=10/10=1.0, 得分=10

    结果: L1 (得分更高)
    """
    if not te_dict or not mapping_dict:
        return "UNKNOWN", "Low", None

    # 统计每个大类的详细信息
    te_stats = {}

    for te_type in ["Alu", "L1", "SVA", "HERV"]:
        total = te_dict.get(te_type, 0)
        if total == 0:
            te_stats[te_type] = {
                "total": 0,
                "consistency": 0,
                "score": 0,
                "best_family": None,
                "best_family_count": 0
            }
            continue

        # 从 mapping_dict 提取该大类的家族分布
        family_counts = {}
        for gene, count in mapping_dict.items():
            # 匹配家族名称（如 AluY, L1HS, SVA_A）
            if gene.startswith(te_type):
                family_counts[gene] = count

        if not family_counts:
            te_stats[te_type] = {
                "total": total,
                "consistency": 0,
                "score": 0,
                "best_family": None,
                "best_family_count": 0
            }
            continue

        # 找出主家族
        best_family = max(family_counts.keys(), key=lambda k: family_counts[k])
        best_family_count = family_counts[best_family]

        # 一致性 = 主家族计数 / 该大类总计数
        consistency = best_family_count / total

        # 综合得分 = 计数 × 一致性
        score = total * consistency

        te_stats[te_type] = {
            "total": total,
            "consistency": consistency,
            "score": score,
            "best_family": best_family,
            "best_family_count": best_family_count
        }

    # 找出得分最高的类型
    max_score = max(s["score"] for s in te_stats.values())
    if max_score == 0:
        return "UNKNOWN", "Low", None

    # 获取所有候选类型（得分 > 0）
    candidates = [(te_type, stats) for te_type, stats in te_stats.items() if stats["score"] > 0]
    candidates.sort(key=lambda x: x[1]["score"], reverse=True)

    best_type = candidates[0][0]
    best_stats = candidates[0][1]

    # 计算置信度
    if len(candidates) == 1:
        confidence = "High"
    else:
        # 比较第一名和第二名的得分差距
        second_score = candidates[1][1]["score"]

        # 得分比例
        score_ratio = best_stats["score"] / (best_stats["score"] + second_score + 0.001)

        # 计数比例（用于辅助判断）
        count_ratio = best_stats["total"] / (best_stats["total"] + candidates[1][1]["total"] + 0.001)

        # 置信度判断
        if score_ratio >= 0.7:
            confidence = "High"
        elif score_ratio >= 0.55:
            confidence = "Medium"
        else:
            confidence = "Low"

    # 如果一致性过低，降低置信度
    if best_stats["consistency"] < 0.5:
        confidence = "Low"
    elif best_stats["consistency"] < 0.7 and confidence == "High":
        confidence = "Medium"

    return best_type, confidence, best_stats["best_family"]


def get_primary_te_type(te_dict):
    """从te_dict获取主要TE类型（兼容旧接口）"""
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


def get_ref_base(fasta_handler, chrom, pos):
    """
    从参考基因组获取指定位置的碱基

    Args:
        fasta_handler: pyfaidx.Fasta 对象或 None
        chrom: 染色体名
        pos: 位置 (1-based)

    Returns:
        碱基字符，如果无法获取则返回 'N'
    """
    if fasta_handler is None:
        return 'N'

    try:
        # pyfaidx 使用 1-based 坐标
        # 处理染色体命名差异
        ref_chrom = chrom

        # 尝试原始名称
        if ref_chrom in fasta_handler:
            seq = fasta_handler[ref_chrom][pos-1:pos]
            return str(seq).upper()
        else:
            return 'N'
    except Exception:
        return 'N'


def generate_vcf_record(row, sample_id, fasta_handler=None):
    """生成单条VCF记录"""
    chrom = row['chrom']
    pos = row['pos']
    reads = row.get('reads', 0)
    direction = row.get('direction', 'unknown')

    mapping_dict = row.get('mapping_dict', {})
    te_dict = row.get('te_dict', {})

    # 使用改进的决策算法
    te_type, confidence, best_family = get_primary_te_type_with_confidence(mapping_dict, te_dict)

    te_to_alt = {
        'Alu': '<INS:ME:ALU>',
        'L1': '<INS:ME:LINE1>',
        'SVA': '<INS:ME:SVA>',
        'HERV': '<INS:ME:HERV>'
    }
    alt = te_to_alt.get(te_type, '<INS:ME:UNKNOWN>')

    # 获取参考碱基
    ref = get_ref_base(fasta_handler, chrom, pos)

    info_parts = [
        f"SVTYPE=INS",
        f"SVLEN=.",
        f"MINAME={te_type}",
        f"SUPPORT={reads}",
        f"CONFIDENCE={confidence}",
        f"DIR={direction}"
    ]

    # 添加最佳家族信息
    if best_family:
        info_parts.append(f"MEFAMILY={best_family}")

    # 添加聚类信息（如果有）
    if 'cluster_size' in row and row['cluster_size'] > 1:
        info_parts.append(f"CLUSTERSIZE={row['cluster_size']}")
        info_parts.append(f"CLUSTERRANGE={row['cluster_range']}")

    # 添加平均softclip长度
    if 'avg_softclip_len' in row:
        info_parts.append(f"AVGSCLEN={row['avg_softclip_len']:.1f}")

    # 添加详细映射信息
    if mapping_dict:
        mapping_str = ",".join([f"{k}:{v}" for k, v in sorted(mapping_dict.items(), key=lambda x: -x[1])])
        info_parts.append(f"MAPPING={mapping_str}")

    # 添加各类型计数
    if te_dict:
        te_str = ",".join([f"{k}:{v}" for k, v in te_dict.items()])
        info_parts.append(f"TECOUNT={te_str}")

    record_id = f"{sample_id}_MEI_{chrom}_{pos}_{direction}"
    qual = "."
    filt = "PASS" if reads >= 10 else "LowSupport"

    # 低置信度时添加 FILTER 标记
    if confidence == "Low":
        if filt == "PASS":
            filt = "LowConfidence"
        else:
            filt = filt + ";LowConfidence"

    info = ";".join(info_parts)

    return f"{chrom}\t{pos}\t{record_id}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}"


def natural_sort_key(chrom):
    """
    自然排序key函数，让染色体按 1,2,3...10,11...X,Y 排序
    而不是字符串排序 1,10,11,2,3...
    """
    # 去除chr前缀
    chrom_clean = chrom.replace('chr', '')

    # 定义排序顺序
    order = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}

    if chrom_clean in order:
        return order[chrom_clean]
    else:
        try:
            return int(chrom_clean)
        except ValueError:
            return 999  # 未知染色体放最后


def write_vcf_output(df, output_path, sample_id, reference, fasta_handler=None):
    """写入VCF文件"""
    with open(output_path, 'w', encoding='utf-8') as vcf_file:
        header = generate_vcf_header(reference)
        vcf_file.write(header)

        # 自然排序
        df_sorted = df.copy()
        df_sorted['_sort_key'] = df_sorted['chrom'].apply(natural_sort_key)
        df_sorted = df_sorted.sort_values(['_sort_key', 'pos'])
        df_sorted = df_sorted.drop(columns=['_sort_key'])

        for idx, row in df_sorted.iterrows():
            record = generate_vcf_record(row, sample_id, fasta_handler)
            vcf_file.write(record + "\n")


# ========== 主流程 ==========

def te_pipe(prefix, input_bam, output_dir, sup_read, min_len,
            rm_tmp, config_file, n_threads, min_mapq, cluster_window, ref_fasta=None):
    """主分析流程（v2.0重构版本）"""

    # 确定参考基因组名称和加载FASTA
    fasta_handler = None
    reference = 'unknown'

    if ref_fasta and os.path.exists(ref_fasta):
        # 从FASTA文件名提取参考名称
        reference = get_reference_from_fasta(ref_fasta)
        print(f"Reference: {reference} (from {os.path.basename(ref_fasta)})")

        # 加载FASTA获取REF碱基
        if HAS_PYFAIDX:
            try:
                import pyfaidx
                fasta_handler = pyfaidx.Fasta(ref_fasta)
                print(f"Loaded reference FASTA for REF bases")
            except Exception as e:
                print(f"Warning: Could not load reference FASTA: {e}")
                print("         REF bases will be 'N'")
        else:
            print("Warning: pyfaidx not installed. REF bases will be 'N'.")
            print("         Install with: pip install pyfaidx")
    else:
        # 从BAM header获取参考基因组信息
        reference, has_chr, suggestion = get_reference_info(input_bam)
        print(f"Reference: {suggestion}")
        if reference == 'unknown':
            print("         VCF header will use 'unknown'. Use -r to specify reference FASTA.")
        if ref_fasta:
            print(f"Warning: Reference FASTA not found: {ref_fasta}")

    # 加载配置
    if config_file == "AUTO_DETECT":
        config_file = os.path.join(get_abs_path(), "config.ini")
    config = configparser.ConfigParser()
    config.read(config_file)

    bwa = config["software"]["BWA"]
    te_genome = config["database"]["TE_GENOME"]

    # 构建输出路径
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    tmp_dir = os.path.join(output_dir, "tmp_" + prefix)
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    output_vcf = os.path.join(output_dir, prefix + ".te.result.vcf")

    # Step 1-3: 合并处理（筛选 + 统计断点 + 写FASTQ）
    print("Step 1: Processing BAM (filter + count breakpoints + write FASTQ)...")
    df, combined_fastq = process_bam_once(input_bam, sup_read, min_len, min_mapq, tmp_dir)

    if df.empty:
        print("No breakpoints found after filtering.")
        empty_df = pd.DataFrame(columns=["chrom", "pos", "direction", "reads", "mapping_dict", "te_dict"])
        write_vcf_output(empty_df, output_vcf, prefix, reference, fasta_handler)
        print(f"VCF output: {output_vcf}")
        return

    print(f"  Found {len(df)} breakpoints before clustering")

    # Step 2: 断点聚类
    print(f"Step 2: Clustering breakpoints (window={cluster_window}bp)...")
    df = cluster_breakpoints(df, cluster_window)
    print(f"  {len(df)} breakpoints after clustering")

    if df.empty:
        print("No breakpoints after clustering.")
        empty_df = pd.DataFrame(columns=["chrom", "pos", "direction", "reads", "mapping_dict", "te_dict"])
        write_vcf_output(empty_df, output_vcf, prefix, reference, fasta_handler)
        print(f"VCF output: {output_vcf}")
        return

    # Step 3: 批量TE比对（管道解析）
    print(f"Step 3: Batch mapping to TE library ({n_threads} threads)...")
    breakpoint_stats = batch_te_mapping(combined_fastq, te_genome, tmp_dir, bwa, n_threads)

    # 将结果合并到DataFrame（使用字典直接赋值，避免apply）
    mapping_dicts = []
    te_dicts = []

    for idx, row in df.iterrows():
        key = f"{row['chrom']}_{row['pos']}_{row['direction']}"
        stats = breakpoint_stats.get(key, {})
        mapping_dicts.append(stats.get("mapping_dict", {}))
        te_dicts.append(stats.get("te_dict", {"SVA": 0, "L1": 0, "Alu": 0, "HERV": 0}))

    df['mapping_dict'] = mapping_dicts
    df['te_dict'] = te_dicts

    # Step 4: 过滤结果
    print("Step 4: Filtering results...")
    df = mapping_filter(df)

    if df.empty:
        print("No TE insertions detected.")
        empty_df = pd.DataFrame(columns=["chrom", "pos", "direction", "reads", "mapping_dict", "te_dict"])
        write_vcf_output(empty_df, output_vcf, prefix, reference, fasta_handler)
        print(f"VCF output: {output_vcf}")
        return

    # Step 5: 输出VCF
    print("Step 5: Writing VCF output...")
    write_vcf_output(df, output_vcf, prefix, reference, fasta_handler)
    print(f"VCF output: {output_vcf}")
    print(f"Total TE insertions detected: {len(df)}")

    # 清理临时文件
    if rm_tmp not in ["False", "false", "FALSE"]:
        import shutil
        shutil.rmtree(tmp_dir)


def main():
    parser = argparse.ArgumentParser(
        description="TIEA-WES: Transposon Insertion Event Analyzer for WES data",
        prog="TIEA-WES.py",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version", version="Version 2.0.0")
    parser.add_argument("-p", "--prefix", type=str, required=True,
                        help="Sample prefix/identifier")
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Input BAM file path")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output directory path")
    parser.add_argument("-r", "--ref_fasta", type=str, default=None,
                        help="Reference genome FASTA file (optional)\n"
                             "  - Provides real REF bases in VCF\n"
                             "  - Sets reference name in VCF header\n"
                             "  - Auto-detected from BAM header if not provided")
    parser.add_argument("-s", "--support_reads", type=int, default=10,
                        help="Minimum breakpoint support reads cutoff (default: 10)")
    parser.add_argument("-l", "--min_length", type=int, default=36,
                        help="Minimum softclip length (default: 36)")
    parser.add_argument("-q", "--min_mapq", type=int, default=20,
                        help="Minimum mapping quality (default: 20)")
    parser.add_argument("-w", "--cluster_window", type=int, default=10,
                        help="Breakpoint clustering window in bp (default: 10)")
    parser.add_argument("-t", "--keep_tmp", type=str, default="False",
                        help="Keep temporary files (default: False)")
    parser.add_argument("-c", "--config", type=str, default="AUTO_DETECT",
                        help="Config file path (default: auto-detect)")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads for BWA alignment (default: 4)")

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
        rm_tmp=args.keep_tmp,
        config_file=args.config,
        n_threads=args.threads,
        min_mapq=args.min_mapq,
        cluster_window=args.cluster_window,
        ref_fasta=args.ref_fasta
    )


if __name__ == "__main__":
    main()