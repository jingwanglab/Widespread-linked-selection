#! /bin/bash -l

chr=$1
##-----Input vcf file 
Inputvcf="/proj/b2011141/nobackup/PaperIV-phylogenomics/GATK/HC/total3/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.vcf.gz"
VCFDir=`dirname $Inputvcf`

Out=${Inputvcf##*/}
echo $Out

##------BEAGLE imputation
beagle="/proj/b2011141/tools/beagle4.1/beagle.08Jun17.d8b.jar"

beagle_dir=$VCFDir/beagle

if [ ! -d "$beagle_dir" ]; then
mkdir -p $beagle_dir
fi

module load bioinfo-tools
module load java

Out_beagle=${Out%.vcf.gz}.gt
Out_beagle_phase=${Out%.vcf.gz}.gt.beagle
##the first is to impute the missing genotype
#java -Xmx100g -jar $beagle gtgl=$Inputvcf out=$beagle_dir/$Out_beagle nthreads=16 window=100000 overlap=10000
#java -Xmx100g -jar $beagle gt=$beagle_dir/${Out_beagle}.vcf.gz out=$beagle_dir/$Out_beagle_phase nthreads=16 window=100000 overlap=10000 ibd=true ibdtrim=100 ibdlod=5
###estimate IBD shared haplotypes
java -Xmx8g -jar $beagle gt=$beagle_dir/${Out_beagle_phase}.vcf.gz out=$beagle_dir/${Out_beagle_phase}.ibd.$chr chrom=$chr ibd=true impute=false window=100000 overlap=10000 ibdtrim=100 ibdlod=5 ibdcm=1e-5


