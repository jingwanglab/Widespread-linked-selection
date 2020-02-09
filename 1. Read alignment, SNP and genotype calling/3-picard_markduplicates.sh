#! /bin/bash -l


module load bioinfo-tools java
module load samtools/0.1.19

dedup="/proj/b2011141/nobackup/alignments/bwa-mem/nisqV3/davidiana/realignment/deduplication"

if [ ! -d "$dedup" ]; then
mkdir -p $dedup
fi

MarkDuplicates="/proj/b2011141/tools/picard-tools-1.115/MarkDuplicates.jar"

if [ "$1" == "" -o "$2" != "" ] ; then
	echo One BAM file at a time
	exit 1
fi

InputBAM=$1

BAM_index=${1%.bam}.bai
if [ -f $BAM_index ] ; then
  rm -f $BAM_index
fi
(samtools index $1 ) &
echo waiting for indexing of BAM files
wait
echo indexing completed

BAMs=${1##*/}
OutputBAM=$dedup/${BAMs%.bam}.dedup.bam
Metrics=$dedup/${BAMs%.bam}.metrics

java -Xmx2500m -jar $MarkDuplicates INPUT=$InputBAM OUTPUT=$OutputBAM METRICS_FILE=$Metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

echo deduplication completed                                    


