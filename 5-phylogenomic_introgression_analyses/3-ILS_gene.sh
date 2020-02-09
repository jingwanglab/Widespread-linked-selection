#! /bin/bash -l

module load bioinfo-tools
module load BEDOPS

###Main aim: the main aim of this script is to estimate calculate the distance of each SNP in vcf file to the nearest gene
Inputvcf="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/GATK/HC/total3/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.vcf.gz"

gene_bed="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/bed/gene/gene.bed"
exon_bed="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/bed/gene/exon.bed"
cds_bed="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/bed/gene/cds.bed"


vcf_Dir=`dirname $Inputvcf`
Out=${Inputvcf##*/}
echo $Out
Outsuffix=${Out%%.vcf.gz}
echo $Outsuffix
gene_distance_dir="$vcf_Dir/gene_distance"

if [ ! -d "$gene_distance_dir" ]; then
mkdir -p $gene_distance_dir
fi

###Step0. convert gff to bed
#convert2bed --input=fmt [--output=fmt] [options] < input > output


#####step1. gunzip the .gz vcf file
gunzip -c $Inputvcf > $gene_distance_dir/${Outsuffix}.vcf

#####step2. tranfer vcf to bed using the awk script
./vcf2bed.awk $gene_distance_dir/${Outsuffix}.vcf > $gene_distance_diss/${Outsuffix}.bed && rm $gene_distance_dir/${Outsuffix}.vcf

#####step3. find the nearest gene to each SNP and calcuate their distance
closest-features --closest --dist $gene_distance_dir/${Outsuffix}.bed $gene_bed > $gene_distance_dir/${Outsuffix}.gene.bed
closest-features --closest --dist $gene_distance_dir/${Outsuffix}.bed $exon_bed > $gene_distance_dir/${Outsuffix}.exon.bed
closest-features --closest --dist $gene_distance_dir/${Outsuffix}.bed $cds_bed > $gene_distance_dir/${Outsuffix}.cds.bed

#####step4. extract the useful information, and only keep the four columns: Chr, Pos, Gene, Distance
cut -d "|" -f 1 $gene_distance_dir/${Outsuffix}.gene.bed | cut -f 1,3 > $gene_distance_dir/snp.txt
cut -d "|" -f 2 $gene_distance_dir/${Outsuffix}.gene.bed | cut -f 4 > $gene_distance_dir/gene.txt
cut -d "|" -f 3 $gene_distance_dir/${Outsuffix}.gene.bed > $gene_distance_dir/distance.txt
paste $gene_distance_dir/snp.txt $gene_distance_dir/gene.txt $gene_distance_dir/distance.txt > $gene_distance_dir/${Outsuffix}.gene.distance.txt && rm $gene_distance_dir/snp.txt $gene_distance_dir/gene.txt $gene_distance_dir/distance.txt
#

##exon
cut -d "|" -f 1 $gene_distance_dir/${Outsuffix}.exon.bed | cut -f 1,3 > $gene_distance_dir/snp.txt
cut -d "|" -f 3 $gene_distance_dir/${Outsuffix}.exon.bed > $gene_distance_dir/distance.txt
paste $gene_distance_dir/snp.txt $gene_distance_dir/distance.txt > $gene_distance_dir/${Outsuffix}.exon.distance.txt && rm $gene_distance_dir/snp.txt $gene_distance_dir/distance.txt

##cds
cut -d "|" -f 1 $gene_distance_dir/${Outsuffix}.cds.bed | cut -f 1,3 > $gene_distance_dir/snp.txt
cut -d "|" -f 3 $gene_distance_dir/${Outsuffix}.cds.bed > $gene_distance_dir/distance.txt
paste $gene_distance_dir/snp.txt $gene_distance_dir/distance.txt > $gene_distance_dir/${Outsuffix}.cds.distance.txt && rm $gene_distance_dir/snp.txt $gene_distance_dir/distance.txt

