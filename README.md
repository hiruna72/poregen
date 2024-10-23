# Poregen

This is a tookit to process ONT raw signal to base alignment.
The following tools are available.

1. `gmove` - collect raw signal event samples for kmers from an event alignment file. The alignment can be signal-to-read or signal-to-reference.

## Help
````

 ./poregen-v0.0.1-linux-x86-64 gmove
Usage: poregen gmove reads.blow5 event_alignment_file output_dir

basic options:
   -k INT                     kmer_size [9]
   -m INT                     move start offset [0]
   -s INT                     kmer start offset [0]
   --scaling INT              scaling [1] (0-no scaling, 1-medmad)
   --margin INT               signal print margin on both sides of the sub signal[0] 
   --sample_limit INT         maximum number of instances to output for a kmer [100] 
   --file_limit INT           maximum number of kmer files to output [500] 
   --kmer_file FILE           kmer file (optional) 
   --index_start INT          1-based closed interval index of start kmer [500] 
   --index_end INT            1-based closed interval index of end kmer [500] 
   --fastq FILE               fastq file (optional - should be provided with .paf) 
   -d                         delimit output files per read
   --max_dur                  maximum move duration allowed for samples [70]
   --min_dur                  maximum move duration allowed for samples [5]
   --pa_min                   minimum pA level a sampling signal should have [40.000]
   --pa_max                  maximum pA level a sampling signal should have [180.000]
   --kmer_pick_margin         distance in bases from an indel when picking a kmer as sample [2]
   --rna                      dataset is rna
   -h                         help
   --verbose INT              verbosity level [3]
   --version                  print version

````
## Quick start

For Linux: download the compiled binaries from the [latest release](https://github.com/hiruna72/poregen/releases).
```sh
wget "https://github.com/hiruna72/poregen/releases/download/v0.0.1/poregen-v0.0.1-linux-x86-64"
./poregen-v0.0.1-linux-x86-64
```

## Compilation and running

```
git clone --recursive https://github.com/hiruna72/poregen.git
cd poregen
mkdir build
make
./poregen --help
```

## Example 
The example command below collects raw signal events for all the 5-mers with a starting offset of 4 (the first 4 events will be skipped). It takes three main arguments: the sequence file (reads.fastq), the signal file (reads.blow5) and the event alignment file (event-alignment.paf). Additionally, each k-mer will have maximum 5000 event samples collected. The raw signals will be converted to pA levels and normalized using Median-Median-Absolute-Difference scaling (med-MAD). The collected event length will be betweeen 19 and 51 signal points.

````
poregen gmove --fastq reads.fastq --k-mer_size 5 -m 4 --sample_limit 5000
        --rna --scaling med-mad --min_dur 19 --max_dur 51 reads.blow5 event-alignment.paf output_dir

````

## Example workflow
### Generating a 5-mer for RNA004 using dorado basecaller's move table
The dataset (the raw signals and the bash scripts) for this workflow is available at [10.5281/zenodo.10966311](https://doi.org/10.5281/zenodo.10966311)

````
# STEP 1 - basecall the dataset to obtain the move table (install buttery-eel-0.4.2 and download ont-dorado-server-7.2.13 and edit the script.sh)
#buttery-eel-0.4.2+ont-dorado-server-7.2.13.script.sh  --moves_out -i ${SIGNAL} --config rna_rp4_130bps_sup.cfg --device cuda:all -o ${SAM} &>> a.log || die "eel failed"

# STEP 2 - the basecalling output sam file has no header; hence add a fake header and convert to bam format. convert bam to fastq and fasta formats.
#echo -e fake_reference'\t'0 > fake_reference.fa.fai
#samtools view ${SAM} -h -t fake_reference.fa.fai -o ${BAM}
#samtools fastq ${BAM} > ${FASTQ} && sed -i '2~4s/N/T/g' ${FASTQ}
#awk 'NR%4==1 {print ">"substr($0, 2)} NR%4==2' ${FASTQ} > ${FASTA}

# STEP 3 - check if the fasta file has enough k-mer (5-mer) coverage
#python count_kmer_freq.py 5 ${FASTA} > ${OUTPUT_DIR}/read_kmer_freq.txt

# STEP 4 - convert the move table to ss format using squigualiser reform
#source squigualiser_venv/bin/activate && squigualiser reform -b ${BAM} --rna -c -o ${REFORM} -k 1 && deactivate

# STEP 5 - collect raw signal event samples for each kmer using poregen gmove program
#${POREGEN} gmove --fastq ${FASTQ} -k 5 --sig_move_offset ${SIG_MOVE_OFFSET} --sample_limit ${SAMPLE_LIMIT} --file_limit 5000 --rna --scaling 1 ${SIGNAL} ${REFORM} ${OUTPUT_DIR}/poregen_output --max_dur ${MAX_DUR} --min_dur ${MIN_DUR} || die "gmove failed"

# STEP 6 - use the collected smaples to calculate median and stddev for each kmer to get a raw k-mer model
#calculate_mean_stddev_all "${OUTPUT_DIR}/poregen_output/dump" "${OUTPUT_DIR}/raw_model" 3.1

# STEP 7 - a) transform the median value to (median*STDDEV)+MEAN where MEAN and STDDEV is the mean and stddev of the all the pA converted signal dataset as a whole.
# b) project stddev values to a custom range [2.5 4]
#apply_transformation "${OUTPUT_DIR}/raw_model" "${OUTPUT_DIR}/transformed_model" 17.569300789355 84.112089074928 2.5 4

# (optional) - set the stddev values to the central 3-mer (XYYYZ) as observed in ONT r9 RNA 5-mer model's stddev values
#set_stddev "${OUTPUT_DIR}/transformed_model" "r9.4_70bps.u_to_t_rna.5mer.template.model" "${OUTPUT_DIR}/transformed_model_adjusted_stddev"

# (optional) - calculate the mean and stddev of the all the pA converted signal dataset as a whole.
#sigtk pa -n ${SIGNAL} | cut -f3 | sed 's/,/\n/g' | datamash mean 1 sstdev 1
# (optional) - process the gmove output to calculate the median dwell time of each k-mer
#calculate_dwell_times_medians "${OUTPUT_DIR}/poregen_output/dump" "${OUTPUT_DIR}/dwell_times"


````

## Acknowledgement
Code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2), [F5c](https://github.com/hasindu2008/f5c), [Nanopolish](https://github.com/jts/nanopolish), and [Sigtk](https://github.com/hasindu2008/sigtk).




