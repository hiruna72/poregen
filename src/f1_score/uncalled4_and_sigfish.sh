#!/bin/bash

REF="Homo_sapiens.GRCh38.cdna.all.fa"
SIGNAL="reads.blow5"
SIGNAL_FAST5="fast5"
FASTQ="reads.fastq"
BAM="sorted_dorado.bam"
POREGEN_MODEL="poregen_RNA004_5mer.model"
SIGFISH=sigfish/sigfish
REGION_FILE="region.file"

./dorado-0.8.2-linux-x64/bin/dorado basecaller -x cuda:all rna004_130bps_fast@v5.1.0 ${SIGNAL_FAST5} --emit-moves --emit-sam --reference ${REF}  --mm2-opts "-x splice -k 14 --secondary=no" | grep -v pi: | samtools sort -o ${BAM}

uncalled4 align -p20 --kit SQK-RNA004 --flowcell FLO-PRO004RA --rna --ref ${REF} --reads ${SIGNAL} --bam-in ${BAM} -o uncalled.bam --bam-f5c
uncalled4 align -p20 --kit SQK-RNA004 --flowcell FLO-PRO004RA -m ${POREGEN_MODEL} --rna --ref ${REF} --reads ${SIGNAL} --bam-in ${BAM} -o uncalled_poregen.bam --bam-f5c

f5c-v1.4/f5c_x86_64_linux eventalign --rna --pore rna004 --secondary=no -g ${REF} -b ${BAM} --slow5 ${SIGNAL} -r ${FASTQ} --sam | samtools sort -o f5c.bam && samtools index f5c.bam

echo "ENST00000343262.8:1-985" > ${REGION_FILE}
cat ${REGION_FILE}
samtools faidx ${REF} --region-file ${REGION_FILE} -o ref.fasta || die "samtools faidx failed"
sed -i "1s/.*/>sim_ref/" ref.fasta

${SIGFISH} dtw -t 8 ref.fasta ${SIGNAL} --sam --rna -q 700 | samtools sort -o sorted_sigfish.bam && samtools index sorted_sigfish.bam
${SIGFISH} dtw -t 8 ref.fasta ${SIGNAL} --sam --rna -q 700 --kmer-model ${POREGEN_MODEL} | samtools sort -o sorted_sigfish_poregen.bam && samtools index sorted_sigfish_poregen.bam

python f1score.py f5c.bam sorted_uncalled.bam --read_limit 1000 --base_shift 0  --threshold 1  --rna
python f1score.py f5c.bam sorted_uncalled_poregen.bam --read_limit 1000 --base_shift -2  --threshold 1  --rna
python f1score.py sorted_sigfish.bam sorted_sigfish_poregen.bam --read_limit 1000 --base_shift -2  --threshold 1  --rna
