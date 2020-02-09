#! /bin/bash -l

module load bioinfo-tools
module load vcftools
module load samtools
module load bcftools

###This script is used to perform XP-CLR analysis for the three aspen species
#In total, we should have six comparisions
#1. P.tra-P.dav(ref)
#2. P.tra-P.trs(ref)
#3. P.dav-P.tra(ref)
#4. P.dav-P.trs(ref)
#5. P.trs-P.tra(ref)
#6. P.trs-P.dav(ref)


step=$1  ###the step to perfrom the analysis

xpclr="/proj/sllstore2017050/nobackup/milou_files/tools/XPCLR/bin/XPCLR"
vcf="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/GATK/HC/total3/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.vcf.gz"
xpclr_snp_input="/proj/sllstore2017050/pipeline_jingwang/pipeline/PaperIV-phylogenomic/XPCLR/xpclr.snp.input.pl"

txt_dir="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/GATK/HC/total3/txt"

OutDir="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/GATK/HC/total3/XPCLR"
if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

OutDir2=$OutDir/input
if [ ! -d "$OutDir2" ]; then
mkdir -p $OutDir2
fi

OutDir3=$OutDir2/out
if [ ! -d "$OutDir3" ]; then
mkdir -p $OutDir3
fi


#####Step1: extracting the SNPs for the pairs of species
if [ "$step" == "1" ]; then
###create the input file of Betascan from the input file from CalcABS
tra="$txt_dir/tremula.txt"
dav="$txt_dir/davidiana.txt"
trs="$txt_dir/tremuloides.txt"

vcftools --gzvcf $vcf --keep $tra --keep $dav --maf 0.00001 --recode --recode-INFO-all --out $OutDir/tra_dav && bgzip $OutDir/tra_dav.recode.vcf
vcftools --gzvcf $vcf --keep $tra --keep $trs --maf 0.00001 --recode --recode-INFO-all --out $OutDir/tra_trs && bgzip $OutDir/tra_trs.recode.vcf
vcftools --gzvcf $vcf --keep $trs --keep $dav --maf 0.00001 --recode --recode-INFO-all --out $OutDir/trs_dav && bgzip $OutDir/trs_dav.recode.vcf

vcftools --gzvcf $OutDir/tra_dav.recode.vcf.gz --keep $tra --recode --recode-INFO-all --out $OutDir/tra_dav.tra && bgzip $OutDir/tra_dav.tra.recode.vcf && tabix $OutDir/tra_dav.tra.recode.vcf.gz
vcftools --gzvcf $OutDir/tra_dav.recode.vcf.gz --keep $dav --recode --recode-INFO-all --out $OutDir/tra_dav.dav && bgzip $OutDir/tra_dav.dav.recode.vcf && tabix $OutDir/tra_dav.dav.recode.vcf.gz

vcftools --gzvcf $OutDir/tra_trs.recode.vcf.gz --keep $tra --recode --recode-INFO-all --out $OutDir/tra_trs.tra && bgzip $OutDir/tra_trs.tra.recode.vcf && tabix $OutDir/tra_trs.tra.recode.vcf.gz
vcftools --gzvcf $OutDir/tra_trs.recode.vcf.gz --keep $trs --recode --recode-INFO-all --out $OutDir/tra_trs.trs && bgzip $OutDir/tra_trs.trs.recode.vcf && tabix $OutDir/tra_trs.trs.recode.vcf.gz

vcftools --gzvcf $OutDir/trs_dav.recode.vcf.gz --keep $dav --recode --recode-INFO-all --out $OutDir/trs_dav.dav && bgzip $OutDir/trs_dav.dav.recode.vcf && tabix $OutDir/trs_dav.dav.recode.vcf.gz
vcftools --gzvcf $OutDir/trs_dav.recode.vcf.gz --keep $trs --recode --recode-INFO-all --out $OutDir/trs_dav.trs && bgzip $OutDir/trs_dav.trs.recode.vcf && tabix $OutDir/trs_dav.trs.recode.vcf.gz

#####Step2: modify the vcf into the input of the XPCLR
elif [ "$step" == "2" ]; then

for new_vcf in $OutDir/*recode.vcf.gz
do
Out=${new_vcf##*/}
echo $Out
Out_Suffix=${Out%.recode.vcf.gz}.geno.txt

for chr in Chr{01..19}
do
bcftools query -f '[%GT ]\n' -r $chr $new_vcf |sed 's/\// /g' |sed 's/\./9/g' > $OutDir2/$chr.$Out_Suffix
done
done

#####Step3 use the R script to convert the population recombination rate estimated from Ldhelmet to 
elif [ "$step" == "3" ]; then

Rscript xpclr_convert_rho_to_cM.R

#####Step4: use R script to convert the position into genetic map position and to the genetic position
elif [ "$step" == "4" ]; then

snp=$2
Rscript interpolate_map_to_snp.R $snp

#####Step4 use XP-CLR to perform analysis
elif [ "$step" == "5" ]; then
##for chr in Chr{01..19}; do if [[ "$chr" =~ "Chr0" ]]; then chr2="${chr/Chr0/}";  elif [[ "$chr" =~ "Chr" ]]; then chr2="${chr/Chr/}"; fi; echo $chr; echo $chr2; sbatch xpclr.aspen_species.sh 4 $chr $chr2; done
#xpclr2=$OutDir2/XPCLR
chr=$2 #Chr01
chr2=$3 #1

cd $OutDir2
hap1="$chr.trs_dav.trs.geno.txt"
hap2="$chr.trs_dav.dav.geno.txt"
snp="$chr.trs_dav.dav.snp.new"
out1="$chr.trs_dav.trs"
out2="$chr.trs_dav.dav"
/proj/sllstore2017050/nobackup/milou_files/tools/XPCLR/bin/XPCLR -xpclr $hap1 $hap2 $snp out/$out1  -w1 0.001 200 1000 $chr2 -p0 0.7
/proj/sllstore2017050/nobackup/milou_files/tools/XPCLR/bin/XPCLR -xpclr $hap2 $hap1 $snp out/$out2  -w1 0.001 200 1000 $chr2 -p0 0.7

elif [ "$step" == "6" ]; then
chr=$2 #Chr01
chr2=$3 #1

cd $OutDir2
hap1="$chr.tra_trs.tra.geno.txt"
hap2="$chr.tra_trs.trs.geno.txt"
snp="$chr.tra_trs.tra.snp.new"
out1="$chr.tra_trs.tra"
out2="$chr.tra_trs.trs"
/proj/sllstore2017050/nobackup/milou_files/tools/XPCLR/bin/XPCLR -xpclr $hap1 $hap2 $snp out/$out1  -w1 0.001 200 1000 $chr2 -p0 0.7
/proj/sllstore2017050/nobackup/milou_files/tools/XPCLR/bin/XPCLR -xpclr $hap2 $hap1 $snp out/$out2  -w1 0.001 200 1000 $chr2 -p0 0.7

elif [ "$step" == "7" ]; then
chr=$2 #Chr01
chr2=$3 #1

cd $OutDir2
hap1="$chr.tra_dav.tra.geno.txt"
hap2="$chr.tra_dav.dav.geno.txt"
snp="$chr.tra_dav.tra.snp.new"
out1="$chr.tra_dav.tra"
out2="$chr.tra_dav.dav"
/proj/sllstore2017050/nobackup/milou_files/tools/XPCLR/bin/XPCLR -xpclr $hap1 $hap2 $snp out/$out1  -w1 0.001 200 1000 $chr2 -p0 0.7
/proj/sllstore2017050/nobackup/milou_files/tools/XPCLR/bin/XPCLR -xpclr $hap2 $hap1 $snp out/$out2  -w1 0.001 200 1000 $chr2 -p0 0.7

fi

