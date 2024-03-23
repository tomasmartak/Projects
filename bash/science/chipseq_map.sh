#!/bin/bash/

# This script takes a fastq file of ChIP-seq data, runs FastQC and outputs a BAM file for it that is ready for peak calling. Bowtie2 is the aligner used, and the outputted BAM file is sorted by genomic coordinates and has duplicate reads removed using sambamba. Lastly, samtools and multiqc are used to generate qc reports.

### USAGE: ###
	# sh path/to script/ChIPseq_read_mapping.sh [path to script fq file target]
	# skeleton (presuming an alias is created): sh chipseq_map $1

# initialize a variable to store the path and name of the target fastq (fq) ChIP-seq file
fq=$1

# variables
genome=~/Documents/ChIP-seq/GRCh38_noalt_as/GRCh38_noalt_as # using hg19/GRCh38 as the reference genome; prebuilt by bowtie2 and not changed in any way; change as needed
blacklist=~/Documents/ChIP-seq/GRCh38_noalt_as/blacklist/ENCFF356LFX.bed # official ENCODE blacklist for hg19; change as needed
threads=7

# grab base of filename for naming outputs and define working directory
base=`basename $fq .fastq`
wd=`dirname "$(dirname "$fq")"` # runs 2 levels down from the file

# Sanity check - here because I'm currently still a noob, might remove later.
echo "These are your chosen sample(s):" $base.fastq
sleep 0.75
echo "This is your chosen directory to deposit the results of this analysis into:" $wd
echo "Make sure that this directory contains a raw folder with the needed unzipped fq/fastq files."
sleep 0.75
echo "This is your reference genome:" $genome
sleep 0.75
echo "This is your current selected blacklist:" $blacklist
sleep 0.75
echo "Press CTRL+Z to exit if you wish."
sleep 0.75
echo "The script will wait 5 seconds before automatically continuing."
sleep 5

# Populate the chosen directory.
# The -p option means mkdir will create the whole path if it does not exist already. While inelegant, I do have a batch version I'm working on that has this before the loop that is elegance personified.
cd $wd
mkdir logs results meta
cd results
mkdir -p fastqc bowtie2/intermediate_bams bowtie2/multiqc macs2 R visualisation
cd ..

# set up output filenames, locations, and aliases
fastqc_out=results/fastqc/
sambamba=~/software/sambamba* #The shell script can't call aliases, so one has to be defined here for simplicity; * added to make general for different versions in the unlikely case of an update anytime soon.

## set up file names
align_out=results/bowtie2/${base}_unsorted.sam
align_bam=results/bowtie2/${base}_unsorted.bam
align_sorted=results/bowtie2/${base}_sorted.bam
align_filtered=results/bowtie2/${base}_sorted_filtered.bam
align_blacklisted=results/bowtie2/${base}_aln.bam
samtool_stats=results/bowtie2/${base}_samStats.txt
bowtie_results=results/bowtie2
intermediate_bams=results/bowtie2/intermediate_bams

# Run FastQC -
fastqc -t $threads $fq  # -t threads
mv raw/*_fastqc.* $fastqc_out
echo "##### FastQC done. #####"

### WARNING!!!
echo "Warning! This script is missing Trimmomatic. This is pretty important in many cases for cleaning up data.

# Run bowtie2 -> aligns FASTQ files to reference genome, outputs SAM files
echo "Running bowtie2 with" $threads "threads; this may take a while. Time started:"
date +"%T"
bowtie2 -p $threads -q --local -x $genome -U $fq -S $align_out # p: $threads number of CPUs; q: specifies fq/fastq format; --local: local read alighment
							# -x: basename of index genome; -U: file(s) to be aligned; -S: file to write SAMs to
echo "##### Bowtie2 mapping done. Time finished: #####"
date +"%T"

# Create BAM from bowtie2 SAM
samtools view -h -S -b -@ $threads -o $align_bam $align_out    # h: include header in output; -S: input is SAM; -b: output is bam;
echo "##### samtools sam to bam conversion complete. #####"
						        # -@ $threads: multithreaded - N $threads; -o: output file
# Sort BAM file by genomic coordinates
$sambamba sort -t $threads -o $align_sorted $align_bam  # sambamba is very noisy, we want the output chucked; -t: threads; -o: output
echo "##### sambamba sort complete. #####"

# Filter out multi-mappers and duplicates
$sambamba view -t $threads -f bam -F "[XS] == null and not unmapped and not duplicate" $align_sorted > $align_filtered # -f: format; -F: set custom filter
echo "##### sambamba duplicate and multi-mapper filter complete. #####"

# Create indices for all the bam files for visualization and QC
samtools index $align_filtered
echo "##### samtools index complete. #####"

# Use bedtools to remove potential blacklisted regions
bedtools intersect -v -a $align_filtered -b $blacklist  > $align_blacklisted # -v: only report those entries in -a that have no overlap in -a
echo "##### bedtools blacklist removal complete. #####"

# Use samtools stats to create alignment statistics
samtools stats $align_blacklisted > $samtool_stats
echo "##### samtools alignment stats created. #####"

# Move intermediate files out of the bowtie2 directory. Bear in mind that not all of these are bams. Sorry not sorry.
mv $bowtie_results/${base}*sorted* $intermediate_bams
echo "##### files cleared up to the intermediate bams directory. #####"

# Use multiqc to generate a report based on the samtool stats outputs
multiqc results/bowtie2/. -n $base -o results/bowtie2/multiqc #-n: custom name; -o: output directory
echo "##### multiQC output generated. #####"

### The script ends here because my bash regex and capture skills are pretty inferior. In time I want to enable a batch mode to also take the results and run peak calling (MACS2) and visualisations (R packages) on data using the Input as a base and Control as a control in MACS.
