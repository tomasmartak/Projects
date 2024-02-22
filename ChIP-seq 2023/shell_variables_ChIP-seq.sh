#!/bin/bash

### Global settings
## CPU threads - some software, e.g. samtools, count the number of CPUs used as n+1. Therefore, if n=6, they will actually use 7.
CPU=6
## Base directories
basedir=~/chipseq/
bd_1=~/chipseq/1_experiment/
bd_2=~/chipseq/2_variations/

### Software paths
## QC
trimmomatic=/usr/share/java/trimmomatic-0.39.jar
## Read mapping
minimap2=~/software/minimap2-2.24_x64-linux/minimap2
## Conversions and indexing
sambamba=~/software/sambamba-0.8.2
igvtools=~/software/IGV_2.16.0/igvtools


### Input and output directories
## genome sequence and annotation folder; subfolders
genome_dir_hg38="./genomes/raw/hg38/"
genome_dir_hg38_p14="./genomes/raw/hg38_p14/"
genome_dir_chm13="./genomes/raw/chm13/"
# subdirs: hg38
genome_dir_bt2_hg38="./genomes/bt2_hg38/"
genome_dir_hs2_hg38="./genomes/hs2_hg38/"
genome_dir_mm2_hg38="./genomes/mm2_hg38/"
# subdirs: chm13
genome_dir_bt2_chm13="./genomes/bt2_chm13/"
genome_dir_hs2_chm13="./genomes/hs2_chm13/"
genome_dir_mm2_chm13="./genomes/mm2_chm13/"
# extras - for partial hg38 regions
genome_dir_hg38_extra="./genomes/hg38_extra/"


### Quality control (QC) directories
# FASTQ file directory containing RAW READ FILES
fastq_dir=~/chipseq/raw/ # this is the directory with my raw read files; replace if somewhere else
fastq_file_ext="\.fastq\.gz$" # file extension specification (unzipped .fastq in this case)
# further QC directories
QC_dir="./QC_ChIP-seq/"
QC_trim_dir="./QC_trim_ChIP-seq/" # timommatic-trimmed reads
fastq_trim_dir="./FASTQ_trim_ChIP-seq/" # timommatic-trimmed FASTQ file directory


### Bam directories
# hs2
bams_hs2_hg38="./BAM_ChIP-seq_hs2_hg38/"
bams_hs2_chm13="./BAM_ChIP-seq_hs2_chm13/"
# bt2
bams_bt2_hg38="./BAM_ChIP-seq_bt2_hg38/"
bams_bt2_chm13="./BAM_ChIP-seq_bt2_chm13/"
# mm2
bams_mm2_hg38="./BAM_ChIP-seq_mm2_hg38/"
bams_mm2_chm13="./BAM_ChIP-seq_mm2_chm13/"


### MACS3 directories
# hs2
peaks_hs2_hg38="./MACS3_hs2_hg38/"
peaks_hs2_chm13="./MACS3_hs2_chm13/"
# bt2
peaks_bt2_hg38="./MACS3_bt2_hg38/"
peaks_bt2_chm13="./MACS3_bt2_chm13/"
# mm2
peaks_mm2_hg38="./MACS3_mm2_hg38/"
peaks_mm2_chm13="./MACS3_mm2_chm13/"


### Script paths
diff_fasta="${basedir}scripts/diff.py"
filter_fasta="${basedir}scripts/filter.py"
region_to_bed="${basedir}scripts/region_to_bed.py"
refseq_to_ucsc_refgen="${basedir}scripts/refseq-to-ucsc-refgen.py"
genbank_to_ucsc_refgen="${basedir}scripts/genbank-to-ucsc-refgen.py"


### here-documents
## hs2_map: read mapper for hisat2
hs2_map=$(cat << 'EOA'
## spider read files and do read mapping
#spider
fastq_files=( $(ls -1 "${fastq_trim_dir}" | grep "${fastq_file_ext}") )
#loop through files and map them; also sort/index them for IGV and output their stats
for i in ${fastq_files[@]};
do
	infile="${fastq_trim_dir}${i}"
	outfile="${bam_dir}${i}.bam"
  echo "${outfile}"
  # samtools view: BAM input (-b) to SAM (default)
  hisat2 -p $CPU -x $index -U $infile --no-spliced-alignment --summary-file "${outfile}.log" \
    | samtools view -@ $CPU -b -F 256 - \
    | $sambamba sort -t $CPU -o $outfile /dev/stdin #sambamba is faster than samtools 
	samtools index -@ $CPU $outfile
	echo "reads in ${i} file:"
	samtools view -@ $CPU -c  $outfile
  samtools idxstat $outfile
done

EOA
)

## bt2_map: read mapper for bowtie2
bt2_map=$(cat << 'EOB'
## spider read files and do read mapping
#spider
fastq_files=( $(ls -1 "${fastq_trim_dir}" | grep "${fastq_file_ext}") )
#loop through files and map them; also sort/index them for IGV and output their stats
for i in ${fastq_files[@]};
do
	infile="${fastq_trim_dir}${i}"
	outfile="${bam_dir}${i}.bam"
  echo "${outfile}"
  # samtools view: BAM input (-b) to SAM (default)
  bowtie2 -p $CPU -x $index -U $infile -S /dev/null 2> "${outfile}.log" \
    | samtools view -@ $CPU -b -F 256 - \
    | $sambamba sort -t $CPU -o $outfile /dev/stdin #sambamba is faster than samtools
	samtools index -@ $CPU $outfile
	echo "reads in ${i} file:"
	samtools view -@ $CPU -c  $outfile
  samtools idxstat $outfile
done
EOB
)