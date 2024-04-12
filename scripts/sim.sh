#!/bin/bash

# steps
# step1
# step2

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

RUN_NO="sim0"
[ "${RUN_NO}" ] || die "RUN_NO is empty"


# REF="Homo_sapiens.GRCh38.cdna.all.fa"
# REF=/genome/gencode.v40.transcripts.fa
REF="ENST00000472952.3.fa"
SIGNAL="reads.blow5"
FASTQ="reads.fastq"
SQUIG_PLOT_DIR="${OUTPUT_DIR}/squig_plot_reads"
REF_REGION="ENST00000361739.1:380-450"
SQUIGULATOR="/data/hiruna/poregen_paper/squigulator-v0.3.0/squigulator"

REMOVE_TMP(){
    rm -f new.blow5 new.fastq a.acc a.log
}

set -x
echo "RNA004 RNA inbuilt rna004-prom default"
${SQUIGULATOR} --seed 100 -x rna004-prom $REF -o new.blow5  2> a.log || die "squigulator failed"
/install/buttery-eel-0.4.2+dorado7.2.13/scripts/eel  -i new.blow5 --config rna_rp4_130bps_sup.cfg --device cuda:all -o new.fastq &>> a.log|| die "eel failed"
identityrna.sh $REF new.fastq > a.acc 2>> a.log || die "identitydna failed"
cat a.acc

echo "RNA004 RNA inbuilt rna-r9-prom"
${SQUIGULATOR} --seed 100 -x rna-r9-prom $REF -o new.blow5  2> a.log || die "squigulator failed"
/install/buttery-eel-0.4.2+dorado7.2.13/scripts/eel  -i new.blow5 --config rna_rp4_130bps_sup.cfg --device cuda:all -o new.fastq &>> a.log|| die "eel failed"
identityrna.sh $REF new.fastq > a.acc 2>> a.log || die "identitydna failed"
cat a.acc

echo "RNA004 RNA ONT 9-mer"
KMER_MODEL="kmer_models/rna04_ont_9mer_levels_v1_transformed_3.0.txt"
${SQUIGULATOR} --seed 100 -x rna004-prom --kmer-model ${KMER_MODEL}  $REF -o new.blow5  2> a.log || die "squigulator failed"
/install/buttery-eel-0.4.2+dorado7.2.13/scripts/eel  -i new.blow5 --config rna_rp4_130bps_sup.cfg --device cuda:all -o new.fastq &>> a.log|| die "eel failed"
identityrna.sh $REF new.fastq > a.acc 2>> a.log || die "identitydna failed"
cat a.acc

#echo "RNA004 r9 DNA 5-mer"
#KMER_MODEL="kmer_models/r9_DNA_with_RNA_header.model"
#${SQUIGULATOR} --seed 100 -x rna004-prom --kmer-model ${KMER_MODEL}  $REF -o new.blow5  2> a.log || die "squigulator failed"
#/install/buttery-eel-0.4.2+dorado7.2.13/scripts/eel  -i new.blow5 --config rna_rp4_130bps_sup.cfg --device cuda:all -o new.fastq &>> a.log|| die "eel failed"
#identityrna.sh $REF new.fastq > a.acc 2>> a.log || die "identitydna failed"
#cat a.acc
#
#echo "RNA004 RNA poregen 5-mer"
#KMER_MODEL="localdataset_repeat/transformed_model"
#${SQUIGULATOR} --seed 100 -x rna004-prom --kmer-model ${KMER_MODEL}  $REF -o new.blow5  2> a.log || die "squigulator failed"
#/install/buttery-eel-0.4.2+dorado7.2.13/scripts/eel  -i new.blow5 --config rna_rp4_130bps_sup.cfg --device cuda:all -o new.fastq &>> a.log|| die "eel failed"
#identityrna.sh $REF new.fastq > a.acc 2>> a.log || die "identityrna failed"
#cat a.acc
#
echo "RNA004 RNA ONT 9-mer stddev adjusted"
KMER_MODEL="kmer_models/rna004.nucleotide.9mer.model.updated"
${SQUIGULATOR} --seed 100 -x rna004-prom --kmer-model ${KMER_MODEL} $REF -o new.blow5  2> a.log || die "squigulator failed"
/install/buttery-eel-0.4.2+dorado7.2.13/scripts/eel  -i new.blow5 --config rna_rp4_130bps_sup.cfg --device cuda:all -o new.fastq &>> a.log|| die "eel failed"
identityrna.sh $REF new.fastq > a.acc 2>> a.log || die "identitydna failed"
cat a.acc

#echo "RNA004 RNA poregen 5-mer stddev adjusted"
#KMER_MODEL="localdataset_repeat/transformed_model_adjusted_stddev"
#${SQUIGULATOR} --seed 100 -x rna004-prom --kmer-model ${KMER_MODEL}  $REF -o new.blow5  2> a.log || die "squigulator failed"
#/install/buttery-eel-0.4.2+dorado7.2.13/scripts/eel  -i new.blow5 --config rna_rp4_130bps_sup.cfg --device cuda:all -o new.fastq &>> a.log|| die "eel failed"
#identityrna.sh $REF new.fastq > a.acc 2>> a.log || die "identityrna failed"
#cat a.acc

echo "RNA004 RNA poregen 5-mer stddev adjusted"
KMER_MODEL="uhr_prom_rna004_2/transformed_model"
${SQUIGULATOR} --seed 100 -x rna004-prom --kmer-model ${KMER_MODEL}  $REF -o new.blow5  2> a.log || die "squigulator failed"
/install/buttery-eel-0.4.2+dorado7.2.13/scripts/eel  -i new.blow5 --config rna_rp4_130bps_sup.cfg --device cuda:all -o new.fastq &>> a.log|| die "eel failed"
identityrna.sh $REF new.fastq > a.acc 2>> a.log || die "identityrna failed"
cat a.acc

echo "RNA004 RNA poregen 5-mer stddev adjusted"
KMER_MODEL="uhr_prom_rna004_2/transformed_model_adjusted_stddev"
${SQUIGULATOR} --seed 100 -x rna004-prom --kmer-model ${KMER_MODEL}  $REF -o new.blow5  2> a.log || die "squigulator failed"
/install/buttery-eel-0.4.2+dorado7.2.13/scripts/eel  -i new.blow5 --config rna_rp4_130bps_sup.cfg --device cuda:all -o new.fastq &>> a.log|| die "eel failed"
identityrna.sh $REF new.fastq > a.acc 2>> a.log || die "identityrna failed"
cat a.acc

REMOVE_TMP

info "$(date)"
info "done"
exit 0