# author: Hiruna Samarakoon hiruna@unsw.edu.au
# just the basic implementation. no optimisations!
import pysam
import argparse


def parse_ss_string(ss_string, start_ref_idx, end_ref_idx, start_signal_idx, args):
    """
    Parses an 'ss' encoded string and creates a signal-to-reference mapping.

    Args:
    ss_string (str): The ss encoded string.
    start_ref_idx (int): The starting index of the reference base.
    start_signal_idx (int): The starting index of the signal.

    Returns:
    dict: A mapping of signal positions to reference positions.
    """
    mapping = {}
    signal_pos = start_signal_idx
    ref_pos = end_ref_idx if args.rna else start_ref_idx

    # Ensure the ss string ends with a non-numeric character
    if ss_string[-1].isdigit():
        raise ValueError("Invalid ss string. It should end with a non-numeric character.")
    
    # Parsing loop
    current_number = ""
    for char in ss_string:
        if char.isdigit():
            current_number += char  # Build the number
        else:
            if current_number:
                # Number of elements (could be signal points or bases, depending on context)
                num_elements = int(current_number)
                if char == ',':  # Step to the next reference base
                    for _ in range(num_elements):
                        mapping[signal_pos] = ref_pos
                        signal_pos += 1
                    if args.rna:
                        ref_pos -= 1
                    else:
                        ref_pos += 1
                elif char == 'D':  # Deletion, no signal points to this reference base
                    if args.rna:
                        ref_pos -= num_elements
                    else:
                        ref_pos += num_elements
                elif char == 'I':  # Insertion, signal points are not mapped to any reference
                    for _ in range(num_elements):
                        mapping[signal_pos] = -1  # Insertion, no reference position
                        signal_pos += 1

            # Reset current number after processing
            current_number = ""

    return mapping

def parse_region(region):
    """Parse a region string in the format 'chr:start-end'."""
    try:
        chrom, positions = region.replace(",", "").split(":")
        start, end = map(int, positions.split("-"))
        return chrom, start, end
    except ValueError:
        raise ValueError("Region must be in the format 'chr:start-end'.")

def compare_mappings(signal_pos_1, ref_pos_1, signal_pos_2, ref_pos_2, args):
    # Initialize TP, FP, FN
    TP = FP = TN = FN = 0

    if args.region:
        _, args_ref_start, args_ref_end = parse_region(args.region)

    # Initialize pointers for both mappings
    idx1, idx2 = 0, 0
    len_1, len_2 = len(signal_pos_1), len(signal_pos_2)

    if signal_pos_1[idx1] < signal_pos_2[idx2]:
        # Advance mapping_1 until the first signal position matches mapping_2
        while idx1 < len_1 and signal_pos_1[idx1] < signal_pos_2[0]:
            idx1 += 1
    elif signal_pos_1[idx1] > signal_pos_2[idx2]:
        # Advance mapping_2 until the first signal position matches mapping_1
        while idx2 < len_2 and signal_pos_2[idx2] < signal_pos_1[0]:
            idx2 += 1

    # Now iterate through both mappings
    while idx1 < len_1 and idx2 < len_2:
        sig1, ref1 = signal_pos_1[idx1], ref_pos_1[idx1]
        sig2, ref2 = signal_pos_2[idx2], ref_pos_2[idx2]
        # print(sig1,sig2,ref1,ref2)
        # print(ref1,args_ref_start,args_ref_end)

        if args.region:
            if args_ref_start > ref1 + 1 or args_ref_end < ref1 + 1:
                idx1 += 1
                idx2 += 1
                continue

        if sig1 == sig2:
            if ref1 == -1 and ref2 == -1:
                TN += 1 # Truen Negative: both did not map to any reference position
            elif ref1 == -1 and ref2 != -1:
                FP += 1  # False Positive: signal point should not match to any reference position
            elif ref1 != -1 and ref2 == -1:
                FN += 1  # False Negative: signal point should match to a reference position
            if abs(ref1 - ref2) <= args.threshold:
                TP += 1  # True Positive: both map to the same referene position
            else:
                FP += 1 # both map to different reference positions
        else:
            print("error: sig1 != sig2. check logic", sig1, sig2)
            exit()
        idx1 += 1
        idx2 += 1

    # Return the results
    return TP, FP, TN, FN

def calculate_metrics(TP, FP, TN, FN):
    # Calculate precision
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    # Calculate recall (sensitivity or TPR)
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    # Calculate F1 score
    F1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    # Calculate specificity (TNR)
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0.0
    # Calculate accuracy
    accuracy = (TP + TN) / (TP + FP + TN + FN) if (TP + FP + TN + FN) > 0 else 0.0

    return precision, recall, F1_score, specificity, accuracy

def evaluate_alignments(alignment_1, alignment_2, args):
    ss1, start_ref_idx1, end_ref_idx1, start_signal_idx1 = alignment_1[0], alignment_1[1], alignment_1[2], alignment_1[3]
    ss2, start_ref_idx2, end_ref_idx2, start_signal_idx2 = alignment_2[0], alignment_2[1], alignment_2[2], alignment_2[3]

    mapping_1 = parse_ss_string(ss1, start_ref_idx1, end_ref_idx1, start_signal_idx1, args)
    mapping_2 = parse_ss_string(ss2, start_ref_idx2, end_ref_idx2, start_signal_idx2, args)

    # Sort the mappings
    sorted_mapping_1 = dict(sorted(mapping_1.items()))
    sorted_mapping_2 = dict(sorted(mapping_2.items()))

    # Create four arrays
    signal_pos_1 = list(sorted_mapping_1.keys())  # Keys of mapping_1
    ref_pos_1 = list(sorted_mapping_1.values())  # Values of mapping_1
    signal_pos_2 = list(sorted_mapping_2.keys())  # Keys of mapping_2
    ref_pos_2 = list(sorted_mapping_2.values())  # Values of mapping_2

    TP, FP, TN, FN = compare_mappings(signal_pos_1, ref_pos_1, signal_pos_2, ref_pos_2, args)

    return TP, FP, TN, FN


def load_file_to_dict(file_path, region=None):
    read_dict = {}
    
    # Open the BAM or SAM file
    with pysam.AlignmentFile(file_path, "r") as file:

        # If a region is provided, parse it
        if region:
            chrom, start, end = parse_region(region)
            records = file.fetch(chrom, start, end)
        else:
            records = file

        for record in records:
            # Check if the record is a primary alignment and a forward mapped alignment
            if not record.is_secondary and not record.is_reverse:
                # Extract read ID
                read_id = record.query_name
                
                # Check if both 'ss' and 'si' tags exist, otherwise raise an error
                if not record.has_tag('ss'):
                    raise ValueError(f"'ss' tag not found in record with read ID {read_id}")
                if not record.has_tag('si'):
                    raise ValueError(f"'si' tag not found in record with read ID {read_id}")
                
                # Extract mapping position and auxiliary tags
                mapping_position = (record.reference_start, record.reference_end)
                ss_tag = record.get_tag('ss')
                si_tag = record.get_tag('si')
                
                # Store in the dictionary
                read_dict[read_id] = (mapping_position, ss_tag, si_tag)
    
    return read_dict

def compare_files(args):
    # Load the first file into a dictionary
    dict1 = load_file_to_dict(args.bam1, args.region)
    # Load the second file into another dictionary
    dict2 = load_file_to_dict(args.bam2, args.region)

    TOT_TP = TOT_FP = TOT_TN = TOT_FN = 0
    read_count = 0

    # Compare the two dictionaries
    for read_id in dict1:
        if args.read_id and read_id != args.read_id:
            continue
        if read_id in dict2:
            pos1, ss1, si1_string = dict1[read_id]
            pos2, ss2, si2_string = dict2[read_id]
            # print(f"{read_id} {si1_string} | {si2_string}")
            si1 = tuple(map(int, si1_string.split(',')))
            si2 = tuple(map(int, si2_string.split(',')))
            # print(si1,si2)
            start_ref_idx1 = si1[3] if args.rna else si1[2]
            start_ref_idx2 = si2[3] if args.rna else si2[2]

            end_ref_idx1 = si1[2] if args.rna else si1[3]
            end_ref_idx2 = si2[2] if args.rna else si2[3]

            start_signal_idx1 = si1[0]
            start_signal_idx2 = si2[0]

            if args.rna:
                end_ref_idx2 += args.base_shift
            else:
                start_ref_idx2 += args.base_shift


            alignment_1 = (ss1, start_ref_idx1, end_ref_idx1, start_signal_idx1)
            alignment_2 = (ss2, start_ref_idx2, end_ref_idx2, start_signal_idx2)
            
            TP, FP, TN, FN = evaluate_alignments(alignment_1,alignment_2,args)
            TOT_TP += TP
            TOT_FP += FP
            TOT_TN += TN
            TOT_FN += FN
        # Apply the read limit
        read_count += 1
        if args.read_limit is not None and read_count == args.read_limit:
            break
    precision, recall, F1_score, specificity, accuracy = calculate_metrics(TOT_TP, TOT_FP, TOT_TN, TOT_FN)
    print("TP\tFP\tTN\tFN", TOT_TP, TOT_FP, TOT_TN, TOT_FN, sep="\t")
    print("precision\trecall\tF1_score\tspecificity\taccuracy", f"{precision:.3f}", f"{recall:.3f}", f"{F1_score:.3f}", f"{specificity:.3f}", f"{accuracy:.3f}", sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare BAM/SAM files.")
    parser.add_argument("bam1", help="Path to the first BAM/SAM file.")
    parser.add_argument("bam2", help="Path to the second BAM/SAM file.")
    parser.add_argument("--read_limit", type=int, default=100, help="Limit the number of records to process.")
    parser.add_argument("--base_shift", type=int, default=0, help="Base shift to apply to the second alignment file.")
    parser.add_argument("--read_id", type=str, default=None, help="Specific read ID to compare.")
    parser.add_argument("--region", type=str, default=None, help="Genomic region to filter reads, format: 'chr:start-end'.")
    parser.add_argument('--rna', required=False, action='store_true', help="specify for RNA reads")
    parser.add_argument('--threshold', type=int, default=0, help="margin of error between reference positions allowed")

    args = parser.parse_args()

    compare_files(args)

