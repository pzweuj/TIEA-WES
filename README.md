# TIEA-WES

> ⚠️ **IMPORTANT: v2.0 is a complete rewrite and is NOT compatible with previous versions.**
>
> If you need the legacy version (v0.0.1), download it here:
> https://github.com/pzweuj/TIEA-WES/archive/refs/tags/v0.0.1.zip
>
> **Major changes in v2.0:**
> - Single-pass BAM processing (3x faster)
> - Bidirectional softclip detection (5' and 3')
> - MAPQ filtering and breakpoint clustering
> - Pipeline BWA output parsing (no intermediate files)
> - Simplified parameters (reference auto-detected from BAM)
> - VCF-only output format

## Overview

Transposon Insertion Event Analyzer for Whole Exome Sequencing data (TIEA-WES) is a tool designed to detect transposon insertion events from next-generation sequencing data generated through probe-captured Whole Exome Sequencing (WES). The software is optimized to identify insertions of Alu, LINE-1 (L1), SVA, and HERV elements with high efficiency.

## Requirements

- Python 3.6+
- pandas
- pysam
- BWA
- pyfaidx (optional, for real REF bases in VCF)

## Installation

### Option 1: Git Clone

```bash
git clone https://github.com/pzweuj/TIEA-WES.git
```

### Option 2: Docker (Recommended)

```bash
# Build image
docker build -t tiea-wes:2.0 .

# Or pull from registry (if available)
# docker pull pzweuj/tiea-wes:2.0
```

## Databases

### For Native Installation

Prepare the transposon genome reference:

```bash
wget https://hgdownload.soe.ucsc.edu/hubs/RepeatBrowser2020/hg38reps/hg38reps.fa
bwa index hg38reps.fa
```

### For Docker Installation

TE reference library (hg38reps.fa) is automatically downloaded and indexed during image build. No manual setup required.

## Configuration

Edit `config.ini`:

```ini
[software]
BWA = bwa

[database]
TE_GENOME = hg38reps.fa
```

Note: SAMTOOLS is no longer required in v2.0.

## Usage

### Native Installation

```bash
python TIEA-WES.py -h
```

### Docker Usage

```bash
# Build image (includes downloading and indexing TE reference library)
docker build -t tiea-wes:2.0 .

# Basic run - only need to mount data directory
docker run --rm -v /path/to/data:/data tiea-wes:2.0 \
    python /app/TIEA-WES.py -p sample -i /data/sample.bam -o /data/output

# With reference genome FASTA (for real REF bases and correct header)
docker run --rm -v /path/to/data:/data -v /path/to/ref:/ref tiea-wes:2.0 \
    python /app/TIEA-WES.py -p sample -i /data/sample.bam -o /data/output -r /ref/hg19.fa

# With custom parameters
docker run --rm -v /path/to/data:/data tiea-wes:2.0 \
    python /app/TIEA-WES.py \
        -p sample \
        -i /data/sample.bam \
        -o /data/output \
        -s 10 \
        -l 36 \
        -q 20 \
        -w 10 \
        --threads 4

# Or use the convenience script
./docker-run.sh -d /path/to/data -- \
    -p sample -i /data/sample.bam -o /data/output

# Interactive mode
docker run -it --rm -v /path/to/data:/data tiea-wes:2.0 bash
```

### Basic usage

```bash
python TIEA-WES.py -p <sample_id> -i <sample.bam> -o <output_dir>
```

### With explicit reference

```bash
python TIEA-WES.py -p <sample_id> -i <sample.bam> -o <output_dir> -f <ref.fa>
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-p` | Sample prefix/identifier | Required |
| `-i` | Input BAM file | Required |
| `-o` | Output directory | Required |
| `-r` | Reference genome FASTA file | None |
| `-s` | Minimum breakpoint support reads | 10 |
| `-l` | Minimum softclip length | 36 |
| `-q` | Minimum mapping quality (MAPQ) | 20 |
| `-w` | Breakpoint clustering window (bp) | 10 |
| `-t` | Keep temporary files | False |
| `-c` | Config file path | Auto-detect |
| `--threads` | Number of threads for BWA alignment | 4 |

**The `-r` parameter:**
- Provides real REF bases in VCF output
- Sets reference name in VCF header (extracted from filename)
- If not provided, REF will be 'N' and reference name is auto-detected from BAM header

## Algorithm Improvements (v2.0)

### 1. Single-Pass BAM Processing
- Merged Step 1-3 into one pass through the BAM file
- Filters, counts breakpoints, and writes FASTQ simultaneously
- **3x speedup** on BAM processing

### 2. Pipeline BWA Parsing
- Directly parses BWA stdout without writing intermediate SAM file
- **Reduced I/O overhead**

### 3. Bidirectional Softclip Detection
- Detects both 5' (S+M) and 3' (M+S) softclip patterns
- `5_prime`: insertion occurs before reference position
- `3_prime`: insertion occurs after reference position
- **Increased sensitivity**

### 4. MAPQ Filtering
- Filters reads with MAPQ < threshold (default 20)
- **Reduced false positives from low-quality alignments**

### 5. Breakpoint Clustering
- Merges nearby breakpoints within ±window bp
- Records cluster size and range in VCF
- **Cleaner output, reduced redundant calls**

## Output

Output file: `<prefix>.te.result.vcf`

VCF follows standard format with symbolic alleles:
- `<INS:ME:ALU>` - Alu insertion
- `<INS:ME:LINE1>` - LINE-1 insertion
- `<INS:ME:SVA>` - SVA insertion
- `<INS:ME:HERV>` - HERV insertion

INFO fields:
- `SVTYPE=INS` - Structural variant type
- `SVLEN=.` - Length (unknown)
- `MINAME` - Mobile element type (Alu/L1/SVA/HERV)
- `MEFAMILY` - Best matching TE family (e.g., AluY, L1HS)
- `CONFIDENCE` - Classification confidence (High/Medium/Low)
- `DIR` - Insertion direction (5_prime/3_prime)
- `SUPPORT` - Supporting reads count
- `AVGSCLEN` - Average softclip length
- `CLUSTERSIZE` - Number of breakpoints merged in cluster
- `CLUSTERRANGE` - Position range of clustered breakpoints
- `MAPPING` - Detailed family mapping (e.g., `AluY:8,AluSx:2`)
- `TECOUNT` - TE type counts (e.g., `Alu:10,L1:10,SVA:0,HERV:0`)

FILTER fields:
- `PASS` - Passed all filters
- `LowSupport` - Less than 10 supporting reads
- `LowConfidence` - TE type classification confidence is low

### TE Type Decision Algorithm

When multiple TE types have similar counts, the algorithm uses **family consistency scoring**:

1. Calculate consistency for each TE type: `consistency = best_family_count / total_count`
2. Compute score: `score = total_count × consistency`
3. Select the TE type with highest score
4. Assign confidence level based on score ratio between top candidates

Example:
```
mapping_dict = {AluY:8, AluSx:2, L1HS:10}
te_dict = {Alu:10, L1:10}

Alu: score = 10 × (8/10) = 8.0
L1:  score = 10 × (10/10) = 10.0

Result: L1 (High confidence)
```

## License

GNU General Public License v3.0