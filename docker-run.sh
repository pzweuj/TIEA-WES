#!/bin/bash
# TIEA-WES Docker Runner
# Simplified script to run TIEA-WES in Docker container
# Note: TE reference library is built into the image

set -e

# Default values
IMAGE_NAME="tiea-wes:2.0"
DATA_DIR="${TIEA_DATA:-$(pwd)/data}"

# Help message
usage() {
    echo "Usage: $0 [options] -- <TIEA-WES arguments>"
    echo ""
    echo "Options:"
    echo "  -d, --data DIR        Data directory (default: ./data)"
    echo "  -i, --image NAME      Docker image name (default: tiea-wes:2.0)"
    echo "  -h, --help            Show this help message"
    echo ""
    echo "Note: TE reference library (hg38reps.fa) is built into the image."
    echo ""
    echo "Example:"
    echo "  $0 -d /path/to/data -- -p sample -i /data/sample.bam -o /data/output"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--data)
            DATA_DIR="$2"
            shift 2
            ;;
        -i|--image)
            IMAGE_NAME="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check if there are remaining arguments for TIEA-WES
if [ $# -eq 0 ]; then
    usage
fi

# Run Docker (reference is built-in, only mount data)
docker run --rm \
    -v "${DATA_DIR}:/data" \
    "${IMAGE_NAME}" \
    python /app/TIEA-WES.py "$@"