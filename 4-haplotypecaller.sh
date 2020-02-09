#! /bin/bash -l

#$1 is the species name, $2 is the bed order, $3 is the sample name

samtools="/proj/b2011141/tools/samtools-0.1.19/samtools"
REF="/proj/b2011141/nobackup/reference/nisqV3/Ptrichocarpa_v3.0_210.fa"
AlignmentDir="/proj/b2011141/nobackup/PaperIV-phylogenomics/bam/$1"
OutDir="/proj/b2011141/nobackup/PaperIV-phylogenomics/GATK/HC"
bed="/proj/b2011141/nobackup/PaperIV-phylogenomics/bed/chr/all.filter.$2.bed"
GATK="/proj/b2011141/tools/GATK/3.8.0/GenomeAnalysisTK.jar"

InputBAMs=$AlignmentDir/$3.bam

InputBAMs_index=${InputBAMs}.bai
if [ ! -f $InputBAMs_index ]; then
      	$samtools index $InputBAMs
fi

echo InputBAMs=\"$InputBAMs\"


Outfile="$1.$3.$2.gatk.hap.vcf"

java -jar $GATK -R $REF \
 -T HaplotypeCaller \
 -I $InputBAMs \
 --emitRefConfidence GVCF \
 --variant_index_type LINEAR \
 --variant_index_parameter 128000 \
 --allow_potentially_misencoded_quality_scores \
 -L $bed \
 --heterozygosity 0.015 \
 --indel_heterozygosity 0.0025 \
 -o $OutDir/$Outfile

$Tools/bgzip $OutDir/$Outfile
$Tools/tabix -p vcf $OutDir/$Outfile.gz

echo Haplotype calling completed


