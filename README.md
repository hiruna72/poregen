# Poregen

This is a tookit to process ONT raw signal to base alignment.
The following tools are available.

1. `gmove` - collect raw signal event samples for kmers from an event alignment file. The alignment can be signal-to-read or signal-to-reference.

## Help

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

## Compilation and running

```
sudo apt-get install zlib1g-dev   #install zlib development libraries
git clone --recursive https://github.com/hiruna72/poregen.git
cd poregen
make
./poregen --help
```

## Acknowledgement
Code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2), [F5c](https://github.com/hasindu2008/f5c), [Nanopolish](https://github.com/jts/nanopolish), and [Sigtk](https://github.com/hasindu2008/sigtk).




