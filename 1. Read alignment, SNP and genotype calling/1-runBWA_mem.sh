#! /bin/bash -l

Cores=8

module load bioinfo-tools
module load samtools/0.1.19
module load bwa



ref="/proj/b2011141/nobackup/reference/nisqV3/Ptrichocarpa_v3.0_210.fa"

fileIn=${1%.trimmomatic.1.fq.gz}
sample=${fileIn##*/}
ind=$sample

BWA_DIR="/proj/b2011141/nobackup/alignments/bwa-mem/nisqV3/davidiana"

if [ ! -d "$BWA_DIR" ]; then
mkdir -p $BWA_DIR
fi

[ -e $ref.sa ] || $bwa index $ref
[ -e $ref.fai ] || samtools faidx $ref

bwa mem -t $Cores -M $ref $fileIn.trimmomatic.1.fq.gz $fileIn.trimmomatic.2.fq.gz | samtools import $ref.fai - - |samtools sort - $BWA_DIR/$sample.PE
bwa mem -t $Cores -M $ref $fileIn.all-trimmomatic.se.fq.gz | samtools import $ref.fai - -| samtools sort - $BWA_DIR/$sample.SE

samtools flagstat $BWA_DIR/$sample.PE.bam > $BWA_DIR/$sample.PE.stats
samtools flagstat $BWA_DIR/$sample.SE.bam > $BWA_DIR/$sample.SE.stats


