# Poregen

This is a package of tools to process ONT raw signal to base alignment.
At the moment following tools are available.
1. sigb_formater - reformat move table generated by Guppy/[Dorado](https://github.com/hiruna72/slow5-dorado)/[Buttery-eel](https://github.com/Psy-Fer/buttery-eel) according to signal-base alignment [specification](https://hasindu2008.github.io/f5c/docs/output)
2. kmer_freq - calculate kmer frequency in a fastq file
3. gmove - produce a pore model using move table (in development)

## Compilation and running

```
sudo apt-get install zlib1g-dev   #install zlib development libraries
git clone --recursive https://github.com/hiruna72/poregen.git
cd poregen
make
./poregen --help
```
The commands to install zlib development libraries on some popular distributions:

```
On Debian/Ubuntu : sudo apt-get install zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install zlib-devel
On OS X : brew install zlib
```

## Sigb_formater example
1. Run basecaller
```
guppy_basecaller -c [DNA model] -i [INPUT] --moves_out --bam_out --save_path [OUTPUT]
```
2. Merge passed BAM files to create a single BAM file
```
samtools merge pass/*.bam -o pass_bam.bam
```
3. Reformat move table 
```
sigb_formater -m 1 -k 1 -o output.paf pass_bam.bam
```

## Move table explanation (unconfirmed)
Nanopore basecallers output move arrays in SAM/BAM format. The important fields are listed below.
1. read_id
2. basecalled fastq sequence length
3. basecalled fastq sequence
4. stride used in the neural network (down sampling factor)
5. raw signal length
6. raw signal trim offset
7. move table

An example move array looks like the following,
```
110100010101000101011010101111…
```
The number of ones (1) in the move array equals to the fastq sequence length. 
According to the above example the first move corresponds with `1 x stride` signal points. 
The second move corresponds with `2 x stride` signale points. The third `4 x stride`, the fourth `2 x stride` and so on.

## Acknowledgement
Code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2).




