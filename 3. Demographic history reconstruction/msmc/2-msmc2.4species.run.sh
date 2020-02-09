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
msmc2="/proj/snic2016-7-89/tools/msmc2/msmc_2.0.0_linux64bit"

##########Inptu file
group=$1 ###tremula, tremuloides, trichocarpa or davidiana
hap=$2  ###hap2,hap4,hap8
inds=$3  ###the specific individual in each species group

Out="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/out"

OutDir=$Out/$group/$hap

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

OutDir_analysis=$OutDir/analysis

if [ ! -d "$OutDir_analysis" ]; then
mkdir -p $OutDir_analysis
fi

#########for only 2 haplotypes
if [ "$hap" == "hap2" ]; then

$msmc2 -o $OutDir_analysis/$inds $OutDir/$inds.*.msmc2.input

fi


####for 4 haplotypes
if [ "$hap" == "hap4" ]; then

$msmc2 -o $OutDir_analysis/$inds $OutDir/$inds.*.msmc2.input -t 6

fi

####for 8 haplotypes
if [ "$hap" == "hap8" ]; then

$msmc2 -o $OutDir_analysis/$inds $OutDir/$inds.*.msmc2.input -t 6

fi


