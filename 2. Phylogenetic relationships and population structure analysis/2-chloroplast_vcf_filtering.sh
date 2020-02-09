#! /bin/bash -l


module load bioinfo-tools
module load vcftools
module load samtools

Inputvcf=$1
VCFDir=`dirname $1`
QD_filter="/proj/b2011141/pipeline/vcftools/chloroplast/filter.txt"

Out=${Inputvcf##*/}
echo $Out

Out_snp_Suffix=${Out%.vcf.gz}.snp
snp_vcf=$Out_snp_Suffix.recode.vcf.gz

Out_biallelic_Suffix=${Out%.vcf.gz}.biallelic
biallelic_vcf=$Out_biallelic_Suffix.recode.vcf.gz

Out_DP_Suffix=${Out%.vcf.gz}.DP500
DP_vcf=$Out_DP_Suffix.recode.vcf.gz

Out_GQ_Suffix=${Out%.vcf.gz}.GQ30
GQ_vcf=$Out_GQ_Suffix.recode.vcf.gz

Out_missing_Suffix=${Out%.vcf.gz}.rm_missing
RM_vcf=$Out_missing_Suffix.recode.vcf.gz

Out_QD_Suffix=${Out%.vcf.gz}.QD
QD_vcf=$Out_QD_Suffix.recode.vcf.gz

Out_QD_filtering_Suffix=${Out%.vcf.gz}.QD.filter
QD_filtering_vcf=$Out_QD_filtering_Suffix.recode.vcf.gz

Out_non_ref_Suffix=${Out%.vcf.gz}.non_ref
Non_ref_vcf=$Out_non_ref_Suffix.recode.vcf.gz

#ut_maf_Suffix=${Out%.vcf.gz}.maf

#remove indels
vcftools --gzvcf $Inputvcf --remove-indels --out $VCFDir/$Out_snp_Suffix --recode --recode-INFO-all
bgzip $VCFDir/$Out_snp_Suffix.recode.vcf

#Bi-alleleic filtering
vcftools --gzvcf $VCFDir/$snp_vcf --max-alleles 2 --recode --recode-INFO-all --out $VCFDir/$Out_biallelic_Suffix
bgzip $VCFDir/$Out_biallelic_Suffix.recode.vcf

#Depth filtering
vcftools --gzvcf $VCFDir/$biallelic_vcf --minDP 500 --recode --recode-INFO-all --out $VCFDir/$Out_DP_Suffix
bgzip $VCFDir/$Out_DP_Suffix.recode.vcf

#Genotype quality filtering
vcftools --gzvcf $VCFDir/$DP_vcf --minGQ 30 --recode --recode-INFO-all --out $VCFDir/$Out_GQ_Suffix
bgzip $VCFDir/$Out_GQ_Suffix.recode.vcf

#Max-missing filtering
vcftools --gzvcf $VCFDir/$GQ_vcf --max-missing 0.8 --recode --recode-INFO-all --out $VCFDir/$Out_missing_Suffix
bgzip $VCFDir/$Out_missing_Suffix.recode.vcf

#Quality by depth filtering
zcat $VCFDir/$RM_vcf | vcf-annotate -f $QD_filter > $VCFDir/$Out_QD_Suffix.recode.vcf
bgzip $VCFDir/$Out_QD_Suffix.recode.vcf
vcftools --gzvcf $VCFDir/$QD_vcf --remove-filtered MinQD --recode --recode-INFO-all --out $VCFDir/$Out_QD_filtering_Suffix
bgzip $VCFDir/$Out_QD_filtering_Suffix.recode.vcf

#remove maf <0.001
vcftools --gzvcf $VCFDir/$QD_filtering_vcf --maf 0.001 --recode --recode-INFO-all --out $VCFDir/$Out_non_ref_Suffix
bgzip $VCFDir/$Out_non_ref_Suffix.recode.vcf

rm $VCFDir/$snp_vcf $VCFDir/$biallelic_vcf $VCFDir/$DP_vcf $VCFDir/$GQ_vcf $VCFDir/$QD_vcf $VCFDir/$RM_vcf $VCFDir/$QD_filtering_vcf
