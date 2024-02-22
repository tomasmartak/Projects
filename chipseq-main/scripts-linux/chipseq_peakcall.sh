# Usage: sh chipseek_peakcall.sh [name] [dir];   dir --> ChIP/results

# capture info for variables
base=$1
wd=$2

# 


# Call peaks 
cd $wd
macs2 callpeak -t bowtie2/$base*ChIP*aln.bam -c bowtie2/$base*INPUT*aln.bam -f BAM -g hs -n $base --outdir macs2 

