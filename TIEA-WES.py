# coding=utf-8
# pzw
# 20241015
# WES转座子插入突变分析流程
# 可分析Alu LINE-1 SVA HERV
# version 20241017 配置文件版
# version 20241016 初始版本

"""
1，提取出WES数据中的SoftClip Reads，仅选择包含且只包含M和S，并且M和S均在36bp以上的Reads；
2，整理这些Reads的断点信息，只去分析有5条以上Reads支持的断点（待定）；
3，将支持的断点的Reads提取出来，与TE参考库进行比对，获得转座子的信息；
4，获得注释结果是：
    断点 ins 转座子  支持Reads数
"""

import os
import sys
import shutil
import pandas as pd
import pysam
import argparse
import configparser

# 获得脚本路径
def get_abs_path():
    return os.path.dirname(os.path.abspath(__file__))

# 获取符合要求的clip read
def process_read(read, min_len=36):
    # 检查是否为重复reads或补充reads
    if read.is_duplicate or read.is_supplementary:
        return
    
    cigar = read.cigartuples
    # 检查是否有cigar
    if cigar is None:
        return

    # 只允许同时包含2个Cigar值的Reads
    if len(cigar) != 2:
        return

    # 检查Tag组成和长度
    more_than_S_M = False
    length_fail = False
    for op in cigar:
        if op[0] not in [0, 4]:
            more_than_S_M = True
        if op[1] < min_len:
            length_fail = True

    # 只允许同时包含S和M的Reads
    if more_than_S_M:
        return
        
    # 检查长度
    if length_fail:
        return
    
    return read

# 单线程
def single_threads_get_reads(input_bam, output_bam, length=36):
    with pysam.AlignmentFile(input_bam, "rb") as bamfile:
        with pysam.AlignmentFile(output_bam, "wb", template=bamfile) as outbam:
            for read in bamfile:
                result = process_read(read, length)
                if result is not None:
                    outbam.write(result)

# 从CIGAR字符串中提取断点位置
def get_breakpoints(read):
    cigar = read.cigartuples
    chrom = read.reference_name
    # pos = read.reference_start
    # 在0位置后面插入，解析为1位置突变
    pos = read.reference_start + 1
    
    breakpoints = []
    
    current_pos = pos
    for op, length in cigar:
        if op == 0:  # M操作
            current_pos += length
        elif op == 4:  # S操作
            if current_pos > pos:  # 软剪切发生在匹配之后
                breakpoints.append((chrom, current_pos, read.query_name))
            else:  # 软剪切发生在匹配之前
                current_pos += length
    return breakpoints

# 统计断点信息
def count_breakpoints(bam_file):
    breakpoint_counts = {}
    with pysam.AlignmentFile(bam_file, "rb") as bamfile:
        for read in bamfile:
            breakpoints = get_breakpoints(read)
            if breakpoints is not None:
                for chrom, pos, read_name in breakpoints:
                    key = (chrom, pos)
                    if key not in breakpoint_counts:
                        breakpoint_counts[key] = {"reads": 0, "reads_name": set()}
                    breakpoint_counts[key]["reads"] += 1
                    breakpoint_counts[key]["reads_name"].add(read_name)
    
    for key in breakpoint_counts:
        breakpoint_counts[key]["reads_name"] = list(breakpoint_counts[key]["reads_name"])
    return breakpoint_counts

# 整理所有的断点信息形成一个矩阵
def get_breakpoint_df(input_bam, cutoff=10):
    # 统计断点信息
    breakpoint_counts = count_breakpoints(input_bam)
    
    # 将统计结果转换为DataFrame
    df = pd.DataFrame(list(breakpoint_counts.values()), index=pd.MultiIndex.from_tuples(breakpoint_counts.keys(), names=['chrom', 'pos']))
    df.index = df.index.set_levels(df.index.levels[0].astype(str), level=0)  # 确保染色体名称为字符串类型
    df.reset_index(inplace=True)  # 将索引转换为列

    # 重排列顺序
    df = df[['chrom', 'pos', 'reads', 'reads_name']]

    # 筛选出reads数大于cutoff的断点
    df = df[df['reads'] >= cutoff]

    # 只需要特定的染色体
    chrom_list = ["X", "Y", "chrX", "chrY"]
    for i in range(1, 23):
        chrom_list.append("chr" + str(i))
        chrom_list.append(str(i))
    df = df[df['chrom'].isin(chrom_list)]
    df.reset_index(inplace=True)

    return df.copy()

# 获得软剪切部分
def get_soft_clipped_parts(read):
    cigar = read.cigartuples
    sequence = read.query_sequence
    qualities = read.query_qualities
    soft_clipped_parts = []

    start_clip = end_clip = None
    current_pos = 0
    for op, length in cigar:
        if op == 4:  # Soft Clip (S)
            if start_clip is None:
                start_clip = current_pos
            end_clip = current_pos + length
            soft_clipped_parts.append((sequence[start_clip:end_clip], qualities[start_clip:end_clip]))
            start_clip = end_clip = None
        else:
            current_pos += length if op in [0, 1, 2, 7, 8] else 0  # Match/Mismatch (M), Insertion (I), Skipped region from the reference (N), Sequence match (7), Padding (8)
    return soft_clipped_parts

# 将结果拆分为多个子目录，并形成fastq
def split_breakpoint_2_fastq(df, input_bam, tmp_dir):
    for idx, row in df.iterrows():
        chrom = row['chrom']
        pos = row['pos']
        reads_name = row['reads_name']

        # 创建子目录
        output_path = os.path.join(tmp_dir, str(chrom) + "_" + str(pos))
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        # 提取reads
        with pysam.AlignmentFile(input_bam, "rb") as bamfile:
            # 打开FASTQ文件用于写入
            with open(output_path + "/" + str(chrom) + "_" + str(pos) + ".fastq", 'w') as fastqfile:
                # 读取并筛选符合条件的reads
                for read in bamfile:
                    if read.query_name in reads_name:                        
                        soft_clipped_parts = get_soft_clipped_parts(read)
                        for seq, quals in soft_clipped_parts:
                            # 将read转换为FASTQ格式并写入文件
                            fastqfile.write("@{}\n".format(read.query_name))
                            fastqfile.write("{}\n".format(seq))
                            fastqfile.write("+\n")
                            quality_scores = "".join([chr(score + 33) for score in quals])
                            fastqfile.write("{}\n".format(quality_scores))

# 比对TE库
# https://hgdownload.soe.ucsc.edu/hubs/RepeatBrowser2020/hg38reps/hg38reps.fa
def TE_mapping(row, te_reference, tmp_dir, bwa):
    row["mapping_dict"] = {}
    row["te_dict"] = {"SVA": 0, "L1": 0, "Alu": 0, "HERV": 0}
    chrom = row['chrom']
    pos = row['pos']
    fastq = os.path.join(tmp_dir, str(chrom) + "_" + str(pos) + "/" + str(chrom) + "_" + str(pos) + ".fastq")
    sam = os.path.join(tmp_dir, str(chrom) + "_" + str(pos) + "/" + str(chrom) + "_" + str(pos) + ".sam")
    cmd = bwa + " mem -t 2 " + te_reference + " " + fastq + " > " + sam
    os.system(cmd)

    # 整理比对结果
    with pysam.AlignmentFile(sam, "r") as samfile:
        for read in samfile:
            gene = read.reference_name
            if gene is None:
                continue
            row["mapping_dict"].setdefault(gene, 0)
            row["mapping_dict"][gene] += 1

            for j in ["SVA", "L1", "Alu", "HERV"]:
                if j in gene and gene.startswith(j):
                    row["te_dict"][j] += 1

    return row

# 进一步过滤提升效率，仅保留比对到转座子基因的行
def mapping_filter(input_df):
    output_df = input_df[input_df['mapping_dict'].apply(lambda x: len(x) > 0)]
    output_df = output_df[["chrom", "pos", "reads", "mapping_dict", "te_dict"]]
    return output_df

# 注释 使用VEP
def annotation_multi(df, reference, tmp_ensembl_input, output, vep, vep_setting):
    # VEP配置
    singularity = vep_setting["SINGULARITY"]
    image = vep_setting["IMAGE"]
    vep_db = vep_setting["VEP_DB"]
    vep_option = vep_setting["VEP_OPTION"]
    ref_fasta = pysam.FastaFile(reference)
    singularity = True if singularity in ["True", "true", "TRUE"] else False

    # 获得输入输出的路径
    tmp_ensembl_input = os.path.abspath(tmp_ensembl_input)
    input_path = os.path.dirname(tmp_ensembl_input)
    output = os.path.abspath(output)
    output_path = os.path.dirname(output)
    reference = os.path.abspath(reference)
    reference_dir = os.path.dirname(reference)

    # 需要将df的chrom和pos形成txt格式
    ensembl_input = open(tmp_ensembl_input, "w", encoding="utf-8")
    check_list = ["A", "T", "C", "G"]
    for index, row in df.iterrows():
        chrom = row['chrom']
        pos = row['pos']
        ref = ref_fasta.fetch(chrom, pos - 1, pos)
        ref_next = ref_fasta.fetch(chrom, pos, pos + 1)
        base_start = "A"
        for cl in check_list:
            if cl != ref and cl != ref_next:
                base_start = cl
                break
        
        # 替换标签
        base_string = base_start + "ATCGATCGATCG"

        # 其他列可以根据需要添加
        ensembl_input.write(f"{chrom}\t{pos}\t{pos}\t{ref}/{base_string}\t.\n")
    ensembl_input.close()

    # 注释
    cmd = ""
    if singularity:
        cmd += f"singularity exec -B {vep_db}:{vep_db},{input_path}:{input_path},{output_path}:{output_path},{reference_dir}:{reference_dir} {image} "
    cmd += f"{vep} --offline --cache --dir_cache {vep_db} -fa {reference} "
    cmd += f"{vep_option} "
    cmd += f"-i {tmp_ensembl_input} -o {output}"
    print(cmd)
    os.system(cmd)

# 从注释结果中决出转录本
def anno_result_to_filter(input_file, transcript_db):
    df = pd.read_excel(transcript_db, sheet_name=0)
    transcript_dict = {}
    for idx, row in df.iterrows():
        transcript_dict[row["#Gene"]] = row["Transcript"]
    
    # 注释结果
    data = {
        'Location': [],
        'REF_ALLELE': [],
        'Allele': [],
        'Consequence': [],
        'SYMBOL': [],
        'Feature': [],
        'EXON': [],
        'INTRON': [],
        'HGVSc': [],
        'HGVSg': []
    }

    # 读入输入文件
    with open(input_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            Location = line[1]
            REF_ALLELE = line[2]
            Allele = line[3]
            Consequence = line[4]
            Feature = line[7]
            EXON = line[8]
            INTRON = line[9]
            SYMBOL = line[10]
            HGVSc = line[12]
            HGVSg = line[14]

            # 判断基因是否在字典中，不在就不要了
            if SYMBOL not in transcript_dict:
                continue

            # 判断是否是默认转录本，不是就不要了
            if Feature.split(".")[0] != transcript_dict[SYMBOL].split(".")[0]:
                continue

            data["Location"].append(Location)
            data["REF_ALLELE"].append(REF_ALLELE)
            data["Allele"].append(Allele)
            data["Consequence"].append(Consequence)
            data["SYMBOL"].append(SYMBOL)
            data["Feature"].append(Feature)
            data["EXON"].append(EXON)
            data["INTRON"].append(INTRON)
            data["HGVSc"].append(HGVSc)
            data["HGVSg"].append(HGVSg)
    df = pd.DataFrame(data)
    df[["chrom", "pos"]] = df["Location"].str.split(":", expand=True)
    df["MatchKey"] = df["chrom"] + "-" + df["pos"]
    df.drop(["Location", "chrom", "pos"], axis=1, inplace=True)
    return df

# 处理结果df
def handle_result(row, cut_off):
    te_dict = row["te_dict"]
    highest_key_list = []
    highest_reads = 0
    for key in te_dict:
        if te_dict[key] > highest_reads:
            highest_reads = te_dict[key]
    for key in te_dict:
        if highest_reads >= cut_off and highest_reads > 0:
            if te_dict[key] == highest_reads:
                highest_key_list.append(key)
    te_mark = str(row["Allele"])
    hgvsc_result = []
    hgvsg_result = []
    hgvsc = str(row["HGVSc"])
    hgvsg = str(row["HGVSg"])
    for k in highest_key_list:
        if k == "L1":
            k = "LINE-1"
        hgvsc = hgvsc_result.append(hgvsc.replace(te_mark, k) + " Reads:" + str(highest_reads))
        hgvsg = hgvsg_result.append(hgvsg.replace(te_mark, k) + " Reads:" + str(highest_reads))
    row["HGVSc"] = "-" if hgvsc is None else hgvsc.replace(te_mark, "TEs")
    row["HGVSg"] = "-" if hgvsg is None else hgvsg.replace(te_mark, "TEs")
    row["TE_Result_c"] = "-"
    row["TE_Result_g"] = "-"
    if len(highest_key_list) != 0:
        row["TE_Result_c"] = "|".join(hgvsc_result)
        row["TE_Result_g"] = "|".join(hgvsg_result)
    return row
        
# 主流程
def te_pipe(prefix, input_bam, output_dir, sup_read, min_len, result_cutoff, reference, annotate, rm_tmp, config_file):
    ############# 配置 ##############
    if config_file == "AUTO_DETECT":
        config_file = os.path.join(get_abs_path(), "config.ini")
    config = configparser.ConfigParser()
    config.read(config_file)

    software = config["software"]
    database = config["database"]
    vep_setting = config["vep_setting"]

    bwa = software["BWA"]
    samtools = software["SAMTOOLS"]
    vep = software["VEP"]

    te_genome = database["TE_GENOME"]
    transcript = database["TRANSCRIPT"]
    ################################

    ############## 构建输出路径 #################
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    tmp_dir = os.path.join(output_dir, "tmp_" + prefix)
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    output_bam = os.path.join(output_dir, prefix + ".softclip.bam")
    output_xls = os.path.join(output_dir, prefix + ".te.result.txt")
    ensembl_input = os.path.join(tmp_dir, prefix + ".te.ensembl_input.txt")
    ensembl_output = os.path.join(tmp_dir, prefix + ".te.ensembl_anno.txt")
    #############################################

    ################## 获得SoftClip #####################
    single_threads_get_reads(input_bam, output_bam, min_len)
    cmd = samtools + " index " + output_bam
    os.system(cmd)
    df = get_breakpoint_df(output_bam, sup_read)
    split_breakpoint_2_fastq(df, output_bam, tmp_dir)
    df = df.apply(lambda x: TE_mapping(x, te_genome, tmp_dir, bwa), axis=1)
    ####################################################
    
    ################### 结果过滤 ########################
    df = mapping_filter(df)
    ####################################################
    
    ################### 结果注释 ########################
    anno_df = None
    if annotate in ["True", "true"]:
        annotation_multi(df, reference, ensembl_input, ensembl_output, vep, vep_setting)
        anno_df = anno_result_to_filter(ensembl_output, transcript)
        df["MatchKey"] = df["chrom"].astype(str) + "-" + df["pos"].astype(str)
        df_merge = df.merge(anno_df, on="MatchKey", how="left")
        df_merge["EXON"] = df_merge["EXON"].astype(str)
        df_merge["INTRON"] = df_merge["INTRON"].astype(str)
        df_merge = df_merge.rename(columns={"reads": "softclip_reads"})
        df_merge = df_merge.apply(lambda x: handle_result(x, result_cutoff), axis=1)
        df_merge.drop(["MatchKey"], axis=1, inplace=True)
        df_merge.to_csv(output_xls, header=True, index=False, sep="\t")
    else:
        df.to_csv(output_xls, header=True, index=False, sep="\t")
    ####################################################

    # 清理
    if rm_tmp in ["False", "false"]:
        shutil.rmtree(tmp_dir)

def main():
    parser = argparse.ArgumentParser(
        description="使用示例：python3 TE_analysis.py -p sample_id -i input_bam -o output_dir",
        prog="TE_analysis.py",
        usage="python3 TE_analysis.py [-h] -p sample_id -i input_bam -o output_dir",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version", version="Version 0.2 20241017")
    parser.add_argument("-p", "--prefix", type=str, help="sample prefix")
    parser.add_argument("-i", "--input", type=str, help="input bam file")
    parser.add_argument("-o", "--output", type=str, help="output path")
    parser.add_argument("-f", "--ref", type=str, help="reference fasta")
    parser.add_argument("-a", "--anno", type=str, help="annotation, default=True", default="True")
    parser.add_argument("-s", "--sr", type=int, help="breakpoint support reads cutoff, default=10", default=10)
    parser.add_argument("-l", "--length", type=int, help="min softclip length, default=36", default=36)
    parser.add_argument("-r", "--cutoff", type=int, help="insertion support reads cutoff, default=5", default=5)
    parser.add_argument("-t", "--rm_tmp", type=str, help="remove tmp file, default=False", default="False")
    parser.add_argument("-c", "--config", type=str, help="config file, default=base_dir/config.ini", default="AUTO_DETECT")
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    te_pipe(prefix=args.prefix, input_bam=args.input, output_dir=args.output, sup_read=args.sr, min_len=args.length, result_cutoff=args.cutoff, reference=args.ref,  annotate=args.anno, rm_tmp=args.rm_tmp, config_file=args.config)

if __name__ == "__main__":
    main()

# end
