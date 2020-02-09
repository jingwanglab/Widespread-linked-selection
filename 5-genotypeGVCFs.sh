#! /bin/bash -l

REF="/proj/b2011141/nobackup/reference/nisqV3/Ptrichocarpa_v3.0_210.fa"
VCFdir="/proj/b2011141/nobackup/PaperIV-phylogenomics/GATK/HC"
GATK="/proj/b2011141/tools/GATK/3.8.0/GenomeAnalysisTK.jar"

InputVCFs=
for VCFs in $VCFdir/*vcf.gz
do
        InputVCFs="$InputVCFs --variant $VCFs"
done

echo InputVCFs=\"$InputVCFs\"

java -Xmx60g -jar $GATK \
 -R $REF \
 -T GenotypeGVCFs \
 $InputVCFs \
 --heterozygosity 0.015 \
 --indel_heterozygosity 0.0025 \
 -o $VCFdir/4species.gatk.hap.vcf


