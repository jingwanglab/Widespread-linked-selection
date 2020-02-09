#! /bin/bash -l

module load bioinfo-tools
module load BEDTools
module load samtools
module load python3

##########################################################################################################
#The main aim of this method is to use the MSMC2 to infer demography history of the four species
#As before, I also use 2,4,8 haplotypes to perform the analysis
#for 2 hapltoypes, we have 8 individuals
#for 4 haplotypes, we have 1-2,1-3,1-4,1-5,1-6,1-7,1-8,2-3,....,7-8, a total of 28 combinations
#for 8 haplotyeps, we have 1-2-3-4,1-2-3-5,1-2-3-6,1-2-3-7,1-2-3-8,1-2-4-5,1-2-4-6,1-2-4-7,1-2-4-8,1-2-5-6,1-2-5-7,1-2-5-8,1-2-6-7,1-2-6-8,1-3-4-5,1-3-4-6,...,5-6-7-8, a total of 70 combinations


##########################################################################################################
#Step 1: extract the vcf file and use msmc-tools to generate the input for msmc-tools
multihetsep="/proj/snic2016-7-89/tools/msmc-tools-master/generate_multihetsep.py"
msmc2="/proj/snic2016-7-89/tools/msmc2/msmc_2.0.0_linux64bit"
vcftools="/crex1/proj/snic2016-7-89/tools/vcftools-0.1.14/bin/vcftools"

Inputvcf="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/vcf/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.gt.beagle.vcf.gz"
tremula_ind="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/samples/tremula.inds"
tremuloides_ind="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/samples/tremuloides.inds"
trichocarpa_ind="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/samples/trichocarpa.inds"
davidiana_ind="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/samples/davidiana.inds"
chr=$1 #Chr01..Chr19
group=$2 ###tremula, tremuloides, trichocarpa or davidiana
hap=$3  ###hap2,hap4,hap8
bed_mask="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/bed/chr/all.filter.$chr.sort.bed"

Out="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/out"


OutDir=$Out/$group/$hap

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

OutDir_analysis=$Out/$group/analysis

if [ ! -d "$OutDir_analysis" ]; then
mkdir -p $OutDir_analysis
fi

group_sample=${group}_ind
nInd=$(cat ${!group_sample} | wc -l)

#########for only 2 haplotypes
if [ "$hap" == "hap2" ]; then

for i in $(seq 1 $nInd); do sample=$(head -n "$i" ${!group_sample} |tail -n 1 );
$vcftools --gzvcf $Inputvcf --chr $chr --indv $sample --maf 0.0001 --out $OutDir/$sample.$chr --recode --recode-INFO-all
bgzip $OutDir/$sample.$chr.recode.vcf
tabix $OutDir/$sample.$chr.recode.vcf.gz
$multihetsep --mask=$bed_mask $OutDir/$sample.$chr.recode.vcf.gz > $OutDir/$sample.$chr.msmc2.input
rm $OutDir/$sample.$chr.recode*
done

fi


####for 4 haplotypes
if [ "$hap" == "hap4" ]; then

###first extract all vcf files for each samples

for i in $(seq 1 $nInd); do 
sample=$(head -n "$i" ${!group_sample} |tail -n 1 );
$vcftools --gzvcf $Inputvcf --chr $chr --indv $sample --maf 0.0001 --out $OutDir/$sample.$chr --recode --recode-INFO-all
bgzip $OutDir/$sample.$chr.recode.vcf
tabix $OutDir/$sample.$chr.recode.vcf.gz
done

###second work on each individual separately
nInd1=$(echo "$nInd-1"| bc )
for i in $(seq 1 $nInd1); do 
sample1=$(head -n "$i" ${!group_sample} |tail -n 1 );
for j in $(seq 2 $nInd); do
if [ $j -le $i ]
then
 continue
else
sample2=$(head -n "$j" ${!group_sample} |tail -n 1 );
echo "${sample1}_${sample2}"
$multihetsep --mask=$bed_mask $OutDir/$sample1.$chr.recode.vcf.gz $OutDir/$sample2.$chr.recode.vcf.gz > $OutDir/${sample1}_${sample2}.$chr.msmc2.input
fi
done
done
fi
###remove all vcf.gz files
rm $OutDir/*.$chr.recode*

####for 8 haplotypes
if [ "$hap" == "hap8" ]; then

###same as hap4, the first step is to generate the vcf.gz files for each sample

for i in $(seq 1 $nInd); do 
sample=$(head -n "$i" ${!group_sample} |tail -n 1 );
$vcftools --gzvcf $Inputvcf --chr $chr --indv $sample --maf 0.0001 --out $OutDir/$sample.$chr --recode --recode-INFO-all
bgzip $OutDir/$sample.$chr.recode.vcf
tabix $OutDir/$sample.$chr.recode.vcf.gz
done

###second, combinations of four samples to make eight haplotypes for inferring recent demography history

nInd1=$(echo "$nInd-3"| bc )
nInd2=$(echo "$nInd-2"| bc )
nInd3=$(echo "$nInd-1"| bc )
for a in $(seq 1 $nInd1); do
sample1=$(head -n "$a" ${!group_sample} |tail -n 1 );
for b in $(seq 2 $nInd2); do
	if [ $b -le $a ]
	then
 	continue
	else
	sample2=$(head -n "$b" ${!group_sample} |tail -n 1 );
	for c in $(seq 3 $nInd3); do
		if [ $c -le $b ]
		then
		continue
		else
		sample3=$(head -n "$c" ${!group_sample} |tail -n 1 );
		for d in $(seq 4 $nInd); do
			if [ $d -le $c ]
			then
			continue
			else
			sample4=$(head -n "$d" ${!group_sample} |tail -n 1 );
			echo "${sample1}_${sample2}_${sample3}_${sample4}"
			$multihetsep --mask=$bed_mask $OutDir/$sample1.$chr.recode.vcf.gz $OutDir/$sample2.$chr.recode.vcf.gz $OutDir/$sample3.$chr.recode.vcf.gz $OutDir/$sample4.$chr.recode.vcf.gz  > $OutDir/${sample1}_${sample2}_${sample3}_${sample4}.$chr.msmc2.input
fi
done
fi
done
fi
done
done

###remove all vcf files
rm $OutDir/*.$chr.recode.*
fi

