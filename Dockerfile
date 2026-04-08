# TIEA-WES Docker Image
# Transposon Insertion Event Analyzer for WES data
# Version: 2.0.0

FROM python:3.10-slim-bullseye

LABEL maintainer="pzweuj"
LABEL version="2.0.0"
LABEL description="TIEA-WES: Transposon Insertion Event Analyzer for WES data"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Shanghai

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    bwa \
    wget \
    curl \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install --no-cache-dir \
    pandas>=1.3.0 \
    pysam>=0.19.0 \
    pyfaidx>=0.7.0

# Create working directories
WORKDIR /app

# Clone TIEA-WES from GitHub
RUN git clone https://github.com/pzweuj/TIEA-WES.git . && \
    rm -rf .git

# Create reference directory and download TE reference library
RUN mkdir -p /app/reference && \
    cd /app/reference && \
    echo "Downloading TE reference library (hg38reps.fa)..." && \
    wget -q https://hgdownload.soe.ucsc.edu/hubs/RepeatBrowser2020/hg38reps/hg38reps.fa && \
    echo "Building BWA index..." && \
    bwa index hg38reps.fa && \
    echo "TE reference ready."

# Update config.ini to use built-in reference
RUN sed -i 's|TE_GENOME = hg38reps.fa|TE_GENOME = /app/reference/hg38reps.fa|' config.ini

# Create data directory for input/output
RUN mkdir -p /data

WORKDIR /data

# Default command
CMD ["python", "/app/TIEA-WES.py", "--help"]

# Usage:
# Build: docker build -t tiea-wes:2.0 .
#
# Run:
#   docker run --rm -v /path/to/data:/data tiea-wes:2.0 \
#       python /app/TIEA-WES.py -p sample -i /data/sample.bam -o /data/output
#
# Or with custom parameters:
#   docker run --rm -v /path/to/data:/data tiea-wes:2.0 \
#       python /app/TIEA-WES.py -p sample -i /data/sample.bam -o /data/output -s 10 --threads 4