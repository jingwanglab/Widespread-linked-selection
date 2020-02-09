#! /bin/bash -l

module load bioinfo-tools
module load vcftools


####Main aim: Assuming there is 4 species pairs, e.g. (((P.tra,P.dav),P.trs),Ptri), all alleles for each SNP that are the same as the reference (P.tri) are assumed to be ancestral allele A, and the other alleles that were assumed to be derived allele B. Then, for populations samples, we count different site types using the frequency of the derived allele at each sites in each population rather than using the binary counts of fixed sites.

##E.g. C(ABBA)=(1-Pi1)Pi2Pi3(1-Pi4), Pi is the derieved allele frequency in each population
##E.g. C(BABA)=Pi1(1-Pi2)Pi3(1-Pi4), Pi is the derieved allele frequency in each population
# C(BAAA)
# C(ABAA)
# C(AABA)
# C(BBAA)
# C(BBBA)  ###If we conly consider the three aspen species, then we may not need to consider this pattern since no polymorphic sites here

step=$1   ###step to perform analysis, e.g. 1,2,3....
species=$2     ###species to work on, e.g. tremula, tremuloides, davidiana, trichocarpa

Inputvcf="/proj/b2011141/nobackup/PaperIV-phylogenomics/GATK/HC/total3/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.vcf.gz"
Out=${Inputvcf##*/}
echo $Out

vcf_Dir=`dirname $Inputvcf`

#####Step1: calculate the non-reference allele frequence for all the four species
if [ "$step" == "1" ]; then
vcftools --gzvcf $Inputvcf --keep $vcf_Dir/allele_freq/${species}.txt --freq2 --out $vcf_Dir/allele_freq/${species}

#####Step2:summarize the above calculated non-reference alleles from four species into one file
elif [ "$step" == "2" ]; then

for species in {tremula,tremuloides,davidiana,trichocarpa}
do
cut -f 6 $vcf_Dir/allele_freq/${species}.frq |sed '1d' > $vcf_Dir/allele_freq/${species}.cut.frq
done

input=
for species in {tremula,davidiana,tremuloides,trichocarpa}
do
input="$input $vcf_Dir/allele_freq/${species}.cut.frq"
done

echo -e "CHROM\tPOS\ttremula_frq\tdavidiana_frq\ttremuloides_frq\ttrichocarpa_frq" > $vcf_Dir/allele_freq/frq.header
cut -f 1,2 $vcf_Dir/allele_freq/tremula.frq |sed '1d' > $vcf_Dir/allele_freq/chrom
paste $vcf_Dir/allele_freq/chrom $input > $vcf_Dir/allele_freq/temp
cat $vcf_Dir/allele_freq/frq.header $vcf_Dir/allele_freq/temp > $vcf_Dir/allele_freq/species.frq

rm $vcf_Dir/allele_freq/frq.header $vcf_Dir/allele_freq/temp $vcf_Dir/allele_freq/chrom $vcf_Dir/allele_freq/{tremula,tremuloides,davidiana,trichocarpa}.frq $vcf_Dir/allele_freq/{tremula,tremuloides,davidiana,trichocarpa}.cut.frq

fi





