#! /bin/bash -l

module load bioinfo-tools
module load BEDTools
module load vcftools
module load samtools

Inputvcf="/proj/b2011141/nobackup/PaperIV-phylogenomics/GATK/HC/total3/4species.gatk.hap.vcf.gz"
VCFDir=`dirname $Inputvcf`
filter_bed="/proj/b2011141/nobackup/PaperIV-phylogenomics/bed/all.filter.Chr.sort.bed"

################Step1: bgzip the vcf file#################
#bgzip $Inputvcf

################Step2: extract the SNPs and Indels separately from the vcf ######################
Out=${Inputvcf##*/}
echo $Out
OutSNPSuffix=${Out%.vcf.gz}.snp
echo $OutSNPSuffix
OutINDELSuffix=${Out%.vcf.gz}.indel
echo $OutINDELSuffix

vcftools --gzvcf $Inputvcf --remove-indels --out $VCFDir/$OutSNPSuffix --recode --recode-INFO-all
vcftools --gzvcf $Inputvcf --keep-only-indels --out $VCFDir/$OutINDELSuffix --recode --recode-INFO-all
bgzip $VCFDir/${OutSNPSuffix}.recode.vcf
bgzip $VCFDir/${OutINDELSuffix}.recode.vcf

bedtools window -a $Inputvcf -b $VCFDir/${Out%.vcf.gz}.indel.recode.vcf.gz -w 5 -u > $VCFDir/${Out%.vcf.gz}.rm_indel.vcf
###before doing the following, the header should be attached first to the vcf file to run the bedtools
cat $VCFDir/header $VCFDir/${Out%.vcf.gz}.rm_indel.vcf > $VCFDir/temp && mv $VCFDir/temp $VCFDir/${Out%.vcf.gz}.rm_indel.vcf
bgzip $VCFDir/${Out%.vcf.gz}.rm_indel.vcf
bedtools subtract -a $VCFDir/${OutSNPSuffix}.recode.vcf.gz -b $VCFDir/${Out%.vcf.gz}.rm_indel.vcf.gz > $VCFDir/${Out%.vcf.gz}.snp.rm_indel.vcf
cat $VCFDir/header $VCFDir/${Out%.vcf.gz}.snp.rm_indel.vcf > $VCFDir/temp && mv $VCFDir/temp $VCFDir/${Out%.vcf.gz}.snp.rm_indel.vcf

###remove region with extreme coverage and MQ0
bedtools intersect -a $VCFDir/${Out%.vcf.gz}.snp.rm_indel.vcf -b $filter_bed > $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.vcf
cat $VCFDir/header $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.vcf > $VCFDir/temp && mv $VCFDir/temp $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.vcf
bgzip $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.vcf

#Bi-alleleic filtering
vcftools --gzvcf $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.vcf.gz --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --out $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.biallelic
bgzip $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.biallelic.recode.vcf

#Depth filtering
vcftools --gzvcf $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.biallelic.recode.vcf.gz --minDP 5 --recode --recode-INFO-all --out $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.biallelic.DP5
bgzip $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.biallelic.DP5.recode.vcf

#Genotype quality filtering
vcftools --gzvcf $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.biallelic.DP5.recode.vcf.gz --minGQ 10 --recode --recode-INFO-all --out $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.biallelic.DP5.GQ10
bgzip $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.biallelic.DP5.GQ10.recode.vcf

#Max-missing filtering
vcftools --gzvcf $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.biallelic.DP5.GQ10.recode.vcf.gz --max-missing 0.9 --maf 0.00001 --recode --recode-INFO-all --out $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing
bgzip $VCFDir/${Out%.vcf.gz}.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.vcf


