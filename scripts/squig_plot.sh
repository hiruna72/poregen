#!/bin/bash

# steps
# squig plot

RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

# set -x

OUTPUT_DIR="squig_plots"
[ "${OUTPUT_DIR}" ] || die "OUTPUT_DIR is empty"

REF="short_names_transcripts.fa" #reference
# MAP_SAM="mapq60_mapped.bam"
MAP_BAM="subset_mapq60.bam"
FASTQ="20k.fastq"
SIGNAL="20k.blow5"

CALCULATE_OFFSETS_TOOL="squigualiser calculate_offsets"
REFORM_TOOL="squigualiser reform"

REFORMAT_PAF="20k.reform.paf"
RE_REFORM_K=2
RE_REFORM_M=1
RE_REFORMAT_PAF="${OUTPUT_DIR}/re_reform.paf"
MOVES_BAM="pass_merge.bam"


READ_ID="--read_id 00002bcc-f711-46f5-b192-0a40af8d571d"

create_output_dir() {
	test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
	mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"
}

calculate_offsets_read_id() {
	info "calculating offsets for ${READ_ID}..."
	${CALCULATE_OFFSETS_TOOL} ${READ_ID} -o ${OUTPUT_DIR}/"one_read.pdf" -p ${REFORMAT_PAF} -f ${FASTQ} -s ${SIGNAL} -k 9 --rna || die "calculate_offsets failed"
}

re_reform() {
	info "running re-reform..."
	${REFORM_TOOL} -k ${RE_REFORM_K} -m ${RE_REFORM_M} --bam ${MOVES_BAM} -c -o ${RE_REFORMAT_PAF} || die "re reform failed"
}

set -x
#create_output_dir
calculate_offsets_read_id
# re_reform

info "success"