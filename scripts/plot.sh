#!/bin/bash

# steps
# step1
# step2

set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

RUN_NO="plot0"
[ "${RUN_NO}" ] || die "RUN_NO is empty"

#...directories files tools arguments commands clean
OUTPUT_DIR=$RUN_NO
[ "${OUTPUT_DIR}" ] || die "OUTPUT_DIR is empty"
# test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
# mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"

REF="Homo_sapiens.GRCh38.cdna.all.fa"
F5C="/data/hiruna/f5c-v1.4/f5c_x86_64_linux"
BAM="filtered.bam"
SIGNAL="reads.blow5"
FASTQ="reads.fastq"
SQUIG_PLOT_DIR="${OUTPUT_DIR}/squig_plot_reads"
REF_REGION="ENST00000361739.1:380-450"

REFORM_TOOL="squigualiser reform"
REALIGN_TOOL="squigualiser realign"
PLOT_TOOL="squigualiser plot"
PLOT_TRACK_TOOL="squigualiser plot_tracks"
CALCULATE_OFFSETS_TOOL="squigualiser calculate_offsets"

f5c_eventalign() {

#	 $F5C index reads.fastq --slow5 ${SIGNAL} || die "f5c index failed"
	$F5C eventalign --rna --pore rna004 --secondary=no -g ${REF} -b ${BAM} --slow5 ${SIGNAL} -r ${FASTQ} --sam | samtools sort -o ${OUTPUT_DIR}/f5c.bam || die "evenalign failed"
	samtools index ${OUTPUT_DIR}/f5c.bam

	KMER_MODEL="kmer_models/rna04_ont_9mer_levels_v1_transformed_3.0.txt"
	$F5C eventalign --rna --kmer-model $KMER_MODEL --secondary=no -g ${REF} -b ${BAM} --slow5 ${SIGNAL} -r ${FASTQ} --sam | samtools sort -o ${OUTPUT_DIR}/ONT.bam || die "evenalign failed"
	samtools index ${OUTPUT_DIR}/ONT.bam

	KMER_MODEL="uhr_prom_rna004_2/transformed_model"
	$F5C eventalign --rna --kmer-model $KMER_MODEL --secondary=no -g ${REF} -b ${BAM} --slow5 ${SIGNAL} -r ${FASTQ} --sam | samtools sort -o ${OUTPUT_DIR}/poregen.bam || die "evenalign failed"
	samtools index ${OUTPUT_DIR}/poregen.bam

}
plot_eventaligns() {
	mkdir -p "${SQUIG_PLOT_DIR}" || die "Failed creating ${SQUIG_PLOT_DIR}"

	TRACK_COMMAND_FILE=${SQUIG_PLOT_DIR}/track_commands_${FUNCNAME[0]}.txt
	rm -f ${TRACK_COMMAND_FILE}
	echo "num_commands=3" > ${TRACK_COMMAND_FILE}
	echo "plot_heights=*" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup --rna -f ${REF} -s ${SIGNAL} -a ${OUTPUT_DIR}/f5c.bam --region ${REF_REGION} --tag_name F5c_9mer_RNA004 --no_overlap --plot_limit 6 --sig_scale znorm --base_shift 0" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup --rna -f ${REF} -s ${SIGNAL} -a ${OUTPUT_DIR}/ONT.bam --region ${REF_REGION} --tag_name ONT_9mer_RNA004 --no_overlap --plot_limit 6 --sig_scale znorm --base_shift 0" >> ${TRACK_COMMAND_FILE}
	echo "squigualiser plot_pileup --rna -f ${REF} -s ${SIGNAL} -a ${OUTPUT_DIR}/poregen.bam --region ${REF_REGION} --tag_name Poregen_5mer_RNA004 --no_overlap --plot_limit 6 --sig_scale znorm --base_shift -2" >> ${TRACK_COMMAND_FILE}

	cat ${TRACK_COMMAND_FILE}

	TESTCASE="rna004_eventaligns"
	OUTPUT="${OUTPUT_DIR}/testcase_${TESTCASE}"
	${PLOT_TRACK_TOOL} --shared_x -f ${TRACK_COMMAND_FILE} -o ${SQUIG_PLOT_DIR}/${TESTCASE} --tag_name ${TESTCASE} || die "testcase:$TESTCASE failed"

}
#f5c_eventalign
plot_eventaligns


info "$(date)"
info "done"
exit 0