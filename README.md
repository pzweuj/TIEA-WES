# TIEA-WES

## Overview

Transposon Insertion Event Analyzer for Whole Exome Sequencing data (TIEA-WES) is a  tool designed to detect transposon insertion events from next-generation sequencing (NGS) data generated through probe-captured Whole Exome Sequencing (WES). The software is optimized to identify insertions of Alu, LINE-1 (L1), SVA, and HERV elements with high efficiency. 

## Requirements

To use TIEA-WES, you need to have the following software installed:

- Python 3.x
- pandas
- pysam
- BWA
- Samtools
- VEP (Optional)

## Installation

Clone the repository:

```bash
git clone https://github.com/pzweuj/TIEA-WES.git
```

## Databases

Before running TIEA-WES, you need to prepare the  transposon genome. Follow these steps:

```bash
wget https://hgdownload.soe.ucsc.edu/hubs/RepeatBrowser2020/hg38reps/hg38reps.fa
# index the fasta
bwa index hg38reps.fa
```

To ensure consistent naming across different datasets, prepare a mapping table that correlates different transcript names or identifiers used in various sources. 
Here is an example of how to create a transcript name mapping table.

```bash
wget https://github.com/pzweuj/TIEA-WES/raw/refs/heads/main/data/Transcript_20240226.xlsx
```

## Usage

### Modifying the Configuration File

[Example](https://github.com/pzweuj/TIEA-WES/blob/main/config.ini)

```config
[software]
BWA = 
SAMTOOLS = 
VEP = 

[database]
TE_GENOME = 
TRANSCRIPT = 

[vep_setting]
SINGULARITY = 
IMAGE = 
VEP_DB = 
VEP_OPTION = 
```

### Run the python script

```bash
python TIEA-WES.py -h
```

To run TIEA-WES on your dataset, execute the following command:

```bash
python TIEA-WES.py -p <sample id> -i <sample.bam> -o <output dir> -f <ref.fa> 
```

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](https://github.com/pzweuj/TIEA-WES/blob/main/LICENSE) file for details.

