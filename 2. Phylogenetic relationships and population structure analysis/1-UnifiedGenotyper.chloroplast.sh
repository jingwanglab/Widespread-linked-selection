#! /bin/bash -l


module load bioinfo-tools 
module load java
module load samtools
#module load GATK
module load BioPerl/1.6.1

samtools="/proj/b2011141/tools/samtools-0.1.19/samtools"
REF="/proj/b2011141/nobackup/reference/P_trichocarpa/chloroplast/Ptrichocarpa_NC_009143.1.fa"
tremula_AlignmentDir="/proj/b2011141/nobackup/PaperIV-phylogenomics/chloroplast/bam/tremula/deduplication"
tremuloides_AlignmentDir="/proj/b2011141/nobackup/PaperIV-phylogenomics/chloroplast/bam/tremuloides/deduplication"
trichocarpa_AlignmentDir="/proj/b2011141/nobackup/PaperIV-phylogenomics/chloroplast/bam/trichocarpa/deduplication"
davidiana_AlignmentDir="/proj/b2011141/nobackup/PaperIV-phylogenomics/chloroplast/bam/davidiana/deduplication"
OutDir="/proj/b2011141/nobackup/PaperIV-phylogenomics/chloroplast/snp_calling/4species/GATK/UG2"
GATK="/proj/b2011141/tools/GATK/3.8.0/GenomeAnalysisTK.jar"
bed="/proj/b2011141/nobackup/PaperIV-phylogenomics/chloroplast/snp_calling/4species/bed/split/trichocarpa.chloroplast.$1.bed"
Tools="/proj/b2011141/tools"

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

if [ ! -f $REF.fai ] ; then
        samtools faidx $REF
fi

Dictionary=${REF%.*}.dict
if [ ! -f $Dictionary ] ; then
        /proj/b2011141/tools/createSequenceDictionary.pl -v -o $Dictionary $REF
fi

InputBAMs=
for BAMs in $tremula_AlignmentDir/*bam
do

        InputBAMs_index=${BAMs}.bai
        if [ ! -f $InputBAMs_index ]; then
        samtools index $BAMs
        fi
        InputBAMs="$InputBAMs -I $BAMs"
done

for BAMs in $tremuloides_AlignmentDir/*bam
do

        InputBAMs_index=${BAMs}.bai
        if [ ! -f $InputBAMs_index ]; then
        samtools index $BAMs
        fi
        InputBAMs="$InputBAMs -I $BAMs"
done

for BAMs in $trichocarpa_AlignmentDir/*bam
do

        InputBAMs_index=${BAMs}.bai
        if [ ! -f $InputBAMs_index ]; then
        samtools index $BAMs
        fi
        InputBAMs="$InputBAMs -I $BAMs"
done

for BAMs in $davidiana_AlignmentDir/*bam
do

        InputBAMs_index=${BAMs}.bai
        if [ ! -f $InputBAMs_index ]; then
        samtools index $BAMs
        fi
        InputBAMs="$InputBAMs -I $BAMs"
done

echo InputBAMs=\"$InputBAMs\"

Outfile="4species.chloroplast.gatk.ug.$1.vcf"

java -jar /proj/b2011141/tools/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -R $REF \
 -T UnifiedGenotyper \
  $InputBAMs \
 --genotype_likelihoods_model BOTH \
 --allow_potentially_misencoded_quality_scores \
 --sample_ploidy 1 \
 -L $bed \
 -dcov 10000 \
 -o $OutDir/$Outfile \
 --output_mode EMIT_ALL_SITES

bgzip $OutDir/$Outfile
tabix -p vcf $OutDir/$Outfile.gz
echo SNP calling completed


