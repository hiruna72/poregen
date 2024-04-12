#!/bin/bash

# steps
# step1
# step2

set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

RUN_NO="p1"
[ "${RUN_NO}" ] || die "RUN_NO is empty"

#...directories files tools arguments commands clean
OUTPUT_DIR=$RUN_NO
[ "${OUTPUT_DIR}" ] || die "OUTPUT_DIR is empty"
test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"

REF="Homo_sapiens.GRCh38.cdna.all.fa"
F5C="/data/hiruna/f5c-v1.4/f5c_x86_64_linux"
BAM="filtered.bam"
SIGNAL="reads.blow5"
FASTQ="reads.fastq"

# $F5C index reads.fastq --slow5 ${SIGNAL} || die "f5c index failed"
#$F5C eventalign --rna --pore rna004 --secondary=no -g ${REF} -b ${BAM} --slow5 ${SIGNAL} -r ${FASTQ} -o ${OUTPUT_DIR}/a.tsv 2>>${OUTPUT_DIR}/stderr || die "evenalign failed"
#$F5C eventalign --rna --pore r9 --secondary=no -g ${REF} -b ${BAM} --slow5 ${SIGNAL} -r ${FASTQ} -o ${OUTPUT_DIR}/a.tsv 2>>${OUTPUT_DIR}/stderr || die "evenalign failed"

#KMER_MODEL="kmer_models/rna04_ont_9mer_levels_v1_transformed_3.0.txt"
#$F5C eventalign --rna --kmer-model $KMER_MODEL --secondary=no -g ${REF} -b ${BAM} --slow5 ${SIGNAL} -r ${FASTQ} -o ${OUTPUT_DIR}/a.tsv 2>>${OUTPUT_DIR}/stderr || die "evenalign failed"
#KMER_MODEL="kmer_models/pore_model_v0_3.0.txt"
#$F5C eventalign --rna --kmer-model $KMER_MODEL --secondary=no -g ${REF} -b ${BAM} --slow5 ${SIGNAL} -r ${FASTQ} -o ${OUTPUT_DIR}/a.tsv 2>>${OUTPUT_DIR}/stderr || die "evenalign failed"

KMER_MODEL="uhr_prom_rna004_2/transformed_model"
$F5C eventalign --rna --kmer-model $KMER_MODEL --secondary=no -g ${REF} -b ${BAM} --slow5 ${SIGNAL} -r ${FASTQ} -o ${OUTPUT_DIR}/a.tsv 2>>${OUTPUT_DIR}/stderr || die "evenalign failed"


#for file in kmer_models/trained_5mer_models/*.model; do
#	$F5C eventalign --rna --kmer-model $file --secondary=no -g ${REF} -b ${BAM} --slow5 ${SIGNAL} -r ${FASTQ} -o ${OUTPUT_DIR}/a.tsv 2>>${OUTPUT_DIR}/stderr || die "evenalign failed"
#done
#
#KMER_MODEL="kmer_models/r9_DNA_with_RNA_header.model"
#$F5C eventalign --rna --kmer-model $KMER_MODEL --secondary=no -g ${REF} -b ${BAM} --slow5 ${SIGNAL} -r ${FASTQ} -o ${OUTPUT_DIR}/a.tsv 2>>${OUTPUT_DIR}/stderr || die "evenalign failed"

info "$(date)"
info "done"
exit 0