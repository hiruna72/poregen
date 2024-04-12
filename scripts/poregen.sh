#!/bin/bash

# set -x
RED='\033[0;31m' ; GREEN='\033[0;32m' ; NC='\033[0m' # No Color
die() { echo -e "${RED}$1${NC}" >&2 ; echo ; exit 1 ; } # terminate script
info() {  echo ; echo -e "${GREEN}$1${NC}" >&2 ; }
info "$(date)"

RUN_NO="uhr_prom_rna004_2"
[ "${RUN_NO}" ] || die "RUN_NO is empty"
ulimit -n 5000
OUTPUT_DIR=$RUN_NO

SIGNAL="PNXRXX240010_reads_20k.blow5"
SAM="${OUTPUT_DIR}/reads.sam"
BAM="${OUTPUT_DIR}/reads.bam"
REFORM="${OUTPUT_DIR}/reform.paf"
FASTQ="${OUTPUT_DIR}/reads.fastq"
FASTA="${OUTPUT_DIR}/reads.fasta"

POREGEN="poregen-v0.0.1-linux-x86-64"
SIG_MOVE_OFFSET=4
SAMPLE_LIMIT=5000
MIN_DUR=19
MAX_DUR=51

create_output_dir() {
	test -d "${OUTPUT_DIR}" && ask "${OUTPUT_DIR}"
	test -d "$OUTPUT_DIR" && rm -r "$OUTPUT_DIR"
	mkdir "$OUTPUT_DIR" || die "Failed creating $OUTPUT_DIR"
}

calculate_dwell_times_medians() {
    # Check if the correct number of arguments is provided
    if [ "$#" -ne 2 ]; then
        echo "Usage: calculate_dwell_times_medians <directory_path> <output_file>"
        return 1
    fi

    local directory="$1"
    local output_file="$2"

    # Loop through each file in the directory
    for file in "$directory"/*; do
        # Get the file name without the path
        filename=$(basename "$file")

        # Calculate the median of comma counts for the file and append to the output file
        median=$(awk -F';' '{for(i=1;i<=NF;i++) print gsub(",", "", $i)}' "$file" | datamash median 1)
        echo -e "$filename\t$median" >> "$output_file"
    done
}

calculate_mean_stddev_all() {
    # Check if the correct number of arguments is provided
    if [ "$#" -ne 3 ]; then
        echo "Usage: calculate_mean_stddev_all <input_directory> <output_file> <limit>"
        return 1
    fi

    local input_directory="$1"
    local output_file="$2"
    local limit="$3"

    # Clear the output file if it exists
    > "$output_file"

    # Loop through each file in the directory
    for file in "$input_directory"/*; do
        # Check if the file is a regular file
        if [ -f "$file" ]; then
            # Calculate the mean and standard deviation of all values in the file
            mean=$(tr ';,' '\n' < "$file" | tail -n +2 | datamash median 1)
            stddev=$(tr ';,' '\n' < "$file" | tail -n +2 | datamash sstdev 1)

            # Limit the standard deviation if it exceeds the specified limit
            if (( $(echo "$stddev > $limit" | bc -l) )); then
                stddev=$limit
            fi

            # Append the results to the output file
            echo -e "$(basename "$file")\t$mean\t$stddev" >> "$output_file"
        fi
    done
}

apply_transformation() {
    # Check if the correct number of arguments is provided
    if [ "$#" -ne 6 ]; then
        echo "Usage: apply_transformation <input_file> <output_file> <A> <B> <C> <D>"
        return 1
    fi

    local input_file="$1"
    local output_file="$2"
    local A="$3"
    local B="$4"
    local C="$5"
    local D="$6"

    # Clear the output file if it exists
    > "$output_file"

    # Write the header to the output file
    echo -e "#ont_model_name\tnone" >> "$output_file"
    echo -e "#kit\tnone" >> "$output_file"
    echo -e "#strand\ttemplate" >> "$output_file"
    echo -e "#k\t5" >> "$output_file"
    echo -e "#alphabet\tnucleotide" >> "$output_file"
    echo -e "#original_file\tnone" >> "$output_file"
#    echo -e "kmer\tlevel_mean\tlevel_stdv" >> "$output_file"
    echo -e "kmer\tlevel_mean\tlevel_stdv\tsd_mean\tsd_stdv\tweight" >> "$output_file"

    # Calculate the minimum and maximum values of the standard deviation column
    read min_stddev max_stddev <<< $(cut -f3 "$input_file" | datamash min 1 max 1)

    # Read each line from the input file and apply the transformation
    while IFS=$'\t' read -r kmer level_mean level_stdv sd_mean sd_stdv weight; do
        # Calculate the transformed mean using the formula (mean * A) + B
        transformed_level_mean=$(echo "($level_mean * $A) + $B" | bc -l)
#        echo $transformed_level_mean

        # Normalize the standard deviation to the range [C, D]
        normalized_level_stdv=$(echo "($level_stdv - $min_stddev) * ($D - $C) / ($max_stddev - $min_stddev) + $C" | bc -l)

        # Append the transformed results to the output file
        echo -e "$kmer\t$transformed_level_mean\t$normalized_level_stdv" >> "$output_file"
    done < "$input_file"
}

set_stddev() {
    # Check if the correct number of arguments is provided
    if [ "$#" -ne 3 ]; then
        echo "Usage: set_stddev <input_file> <r9_rna_model> <output_file>"
        return 1
    fi

    local input_file="$1"
    local r9_model="$2"
    local output_file="$3"

    head -n 7 $input_file > header
    tail -n +8 $input_file | cut -f 1,2 > temp
    tail -n +8 $r9_model | cut -f 3 | paste temp - > temp1
    cat header temp1 > $output_file
    rm temp
    rm temp1
}

#create_output_dir

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

info "$(date)"
info "done"
exit 0
