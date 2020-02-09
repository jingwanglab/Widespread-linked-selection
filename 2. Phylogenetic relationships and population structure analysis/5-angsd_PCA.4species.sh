#! /bin/bash -l


angsd="/proj/b2011141/tools/angsd0.917/angsd"
emOptim2="/proj/b2011141/tools/angsd0.917/misc/emOptim2"
thetaStat="/proj/b2011141/tools/angsd0.917/misc/thetaStat"
bam_list_all="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/PCA/4species.bam.list"
ref="/proj/b2011141/nobackup/reference/nisqV3/Ptrichocarpa_v3.0_210.fa"
OutDir="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/PCA"
region="/proj/b2011141/nobackup/PaperIV-phylogenomics/bed/chr"
chr="$region/all.filter.Chr$1.region"
ngsCovar="/proj/b2011141/tools/ngsTools/ngsPopGen/ngsCovar"

nInd=$(cat $bam_list_all | wc -l)
nChrom=$(echo "2*$nInd" | bc)

echo $nInd
echo $nChrom


if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

#first generate .saf file
$angsd -bam $bam_list_all -GL 1 -doSaf 1 -out $OutDir/4species.Chr$1 -anc $ref -rf $chr -minMapQ 30 -minQ 20

#use emOptim2 to optimization
$emOptim2 $OutDir/4species.Chr$1.saf $nChrom -maxIter 100 > $OutDir/4species.Chr$1.sfs

#calculate posterior probabilities of sample allele frequencies
$angsd -bam $bam_list_all -GL 1 -doSaf 1 -anc $ref -rf $chr -minMapQ 30 -minQ 20 -pest $OutDir/4species.Chr$1.sfs -out $OutDir/4species.Chr$1.rf

#Use ANGSD to compute genotype posterior probabilities
$angsd -bam $bam_list_all -GL 1 -doGeno 32 -doPost 1 -doMaf 1 -rf $chr -out $OutDir/4species.Chr$1 -doMajorMinor 1 -ref $ref -minMapQ 30 -minQ 20

#get covariance matric
input=
for chr in {01..19}
do
input="$input $OutDir/4species.Chr$1.geno.gz"
done
cat $input > $OutDir/4species.all.geno.gz 
# rm $input
gunzip $OutDir/4species.all.geno.gz

input=
for chr in {01..19}
do
input="$input $OutDir/4species.Chr$1.rf.saf"
done
cat $input > $OutDir/4species.all.rf.saf 

$ngsCovar -probfile $OutDir/4species.all.geno -outfile $OutDir/4species.all.covar1 -nind $nInd -call 0 -sfsfile $OutDir/4species.all.rf.saf -norm 0 -nsites 168879808 -block_size 20000 


