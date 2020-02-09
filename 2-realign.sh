#! /bin/bash -l

# Realign reads in a BAM file.
#
#     ./realign.sh file.fa input1.bam [input2.bam ...]
#
# This script takes 2 or more arguments:
#
# 1. The reference FASTA file.  An index and sequence dictionary are
#    created if they do not already exist under names based off that 
#    of the FASTA file, for example, "file.fa.fai" and "file.dict".
#
# 2. One or more BAM files to realign.  If you provide more than one
#    BAM file, it should be from the same individual (for now).  For
#    example, you have a BAM file for a 300bp insert library and one
#    for a 650bp insert library.
#
# It will produce one BAM per input BAM, ending with ".realign.bam".
# It also may not work... GATK is very memory-hungry and may not
# be able to run in the memory available.  If that is the case, see
# the comment below.
#

set -e
set -x

####################

MAXREADS=50000

TMPID=${SLURM_JOB_ID:-$$}

module load bioinfo-tools
module load samtools/0.1.19
module load BioPerl/1.6.1

GATK="/proj/b2011141/tools/GenomeAnalysisTK-3.8.0/GenomeAnalysisTK.jar"

Ref=$1
if [ ! -f $Ref.fai ] ; then
	samtools faidx $Ref
fi
Dictionary=${Ref%.*}.dict
if [ ! -f $Dictionary ] ; then
	/proj/b2011141/tools/createSequenceDictionary.pl -v -o $Dictionary $Ref
fi

shift

InputBAMList=("$@")
OutputBAMSuffix=".realign.bam"

ArgBAMs=
OutputBAMs=
for bam in ${InputBAMList[@]} ; do
	if [ ! -f ${bam%.bam}.bai -a ! -f $bam.bai ] ; then
		( samtools index $bam ) &
	fi
        ArgBAMs="$ArgBAMs -I $bam"
        Out=${bam##*/}
        Out=${Out%.bam}$OutputBAMSuffix
        OutputBAMs="$OutputBAMs $Out"
	LastBAM=$bam
done
echo waiting for indexing of input BAM files
wait
echo indexing completed or unnecessary
echo ArgBAMs=\"$ArgBAMs\"
echo OutputBAMs=\"$OutputBAMs\"

LastBAM=${LastBAM##*/}
Intervals=${LastBAM%%.*}.realign.intervals

java -jar $GATK $ArgBAMs -R $Ref -T RealignerTargetCreator -o $Intervals --allow_potentially_misencoded_quality_scores

if [ ! -f $Intervals ] ; then
	echo Interval files not created
	exit 1
fi

java -jar $GATK $ArgBAMs -R $Ref -T IndelRealigner --targetIntervals $Intervals --maxReadsForRealignment $MAXREADS --nWayOut $OutputBAMSuffix --allow_potentially_misencoded_quality_scores

for outbam in $OutputBAMs ; do
        out=${outbam##*/}
	GATK_index=${out%.bam}.bai
	if [ -f $GATK_index ] ; then
		rm -f $GATK_index
	fi
	( samtools index $outbam ) &
done
echo waiting for indexing of output BAM files
wait
echo indexing completed

