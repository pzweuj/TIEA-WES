# TIEA-WES

## Overview

Transposon Insertion Event Analyzer for Whole Exome Sequencing data (TIEA-WES) is a tool designed to detect transposon insertion events from next-generation sequencing data generated through probe-captured Whole Exome Sequencing (WES). The software is optimized to identify insertions of Alu, LINE-1 (L1), SVA, and HERV elements with high efficiency.

## Requirements

- Python 3.x
- pandas
- pysam
- BWA
- Samtools

## Installation

```bash
git clone https://github.com/pzweuj/TIEA-WES.git
```

## Databases

Prepare the transposon genome reference:

```bash
wget https://hgdownload.soe.ucsc.edu/hubs/RepeatBrowser2020/hg38reps/hg38reps.fa
bwa index hg38reps.fa
```

## Configuration

Edit `config.ini`:

```ini
[software]
BWA = bwa
SAMTOOLS = samtools

[database]
TE_GENOME = hg38reps.fa
```

## Usage

```bash
python TIEA-WES.py -h
```

### Basic usage

```bash
python TIEA-WES.py -p <sample_id> -i <sample.bam> -o <output_dir> -f <ref.fa>
```

### Output VCF format

```bash
python TIEA-WES.py -p <sample_id> -i <sample.bam> -o <output_dir> -f <ref.fa> --vcf
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-p` | Sample prefix/identifier | Required |
| `-i` | Input BAM file | Required |
| `-o` | Output directory | Required |
| `-f` | Reference FASTA (for VCF header) | Required |
| `-s` | Minimum breakpoint support reads | 10 |
| `-l` | Minimum softclip length | 36 |
| `-t` | Keep temporary files | False |
| `-c` | Config file path | Auto-detect |
| `--vcf` | Output VCF format | False |

## Output

### TSV format (default)

Output file: `<prefix>.te.result.txt`

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| pos | Breakpoint position |
| reads | Supporting softclip reads count |
| mapping_dict | TE mapping details (e.g., `{AluY:10, L1HS:3}`) |
| te_dict | TE type summary (e.g., `{Alu:10, L1:3, SVA:0, HERV:0}`) |

### VCF format (--vcf)

Output file: `<prefix>.te.result.vcf`

VCF follows standard format with symbolic alleles:
- `<INS:ME:ALU>` - Alu insertion
- `<INS:ME:LINE1>` - LINE-1 insertion
- `<INS:ME:SVA>` - SVA insertion
- `<INS:ME:HERV>` - HERV insertion

INFO fields:
- `SVTYPE=INS` - Structural variant type
- `SVLEN=.` - Length (unknown)
- `MINAME` - Mobile element type
- `SUPPORT` - Supporting reads count
- `MAPPING` - Detailed mapping info

## Downstream Annotation

VEP annotation can be performed on the VCF output:

```bash
vep -i sample.te.result.vcf -o sample.annotated.vcf --cache --offline
```

## License

GNU General Public License v3.0