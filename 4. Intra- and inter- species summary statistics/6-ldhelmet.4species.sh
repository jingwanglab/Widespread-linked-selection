#! /bin/bash -l

module load bioinfo-tools
module load LDhelmet
module load BEDTools
module load samtools
module load BioPerl/1.6.1
module load vcftools

step=$1
chr=$2
species=$3


gatk="/proj/b2011141/tools/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"

fa="/proj/b2011141/nobackup/reference/nisqV3/Ptrichocarpa_v3.0_210.fa"
filter_bed="/proj/b2011141/nobackup/PaperIV-phylogenomics/bed/chr/all.filter.Chr$chr.bed"
inputvcf="/proj/b2011141/nobackup/PaperIV-phylogenomics/GATK/HC/total3/beagle/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.vcf.gz"
mask_bed="/proj/b2011141/nobackup/PaperIV-phylogenomics/bed/Potri_bed/mask_bed/trichocarpa.mask.bed"

OutDir="/proj/b2011141/nobackup/PaperIV-phylogenomics/recombination_rates/ldhelmet"
OutDir_species="/proj/b2011141/nobackup/PaperIV-phylogenomics/recombination_rates/ldhelmet/$species"

if [ ! -d "$OutDir_species" ]; then
mkdir -p $OutDir_species
fi

OutDir_fa="$OutDir/fa"
if [ ! -d "$OutDir_fa" ]; then
mkdir -p $OutDir_fa
fi

OutDir_out="$OutDir/$species/out"
if [ ! -d "$OutDir_out" ]; then
mkdir -p $OutDir_out
fi


#####################################################################################
######Step1: create the referece fasta file for each chromosome

if [ $step == "1" ]; then
samtools faidx $fa Chr$chr > $OutDir_fa/ref.Chr$chr.fa

Ref="$OutDir_fa/ref.Chr$chr.fa"
echo "$Ref"
if [ ! -f $Ref.fai ] ; then
        samtools faidx $Ref
fi


elif [ $step == "2" ]; then
#####################################################################################
######Step2: use vcf-consensus to extract the fasta file for each haplotype within each sample
species_sample="/proj/b2011141/nobackup/PaperIV-phylogenomics/recombination_rates/ldhelmet/$species.sample.txt"

for i in {1..8}
	do
	sample=$(cat $species_sample |head -n $i |tail -n 1)
	echo $sample
	####extract the haplotype sequences for chromosome for each sample
	cat $OutDir_fa/ref.Chr$chr.fa | vcf-consensus -H 1 -s $sample $inputvcf > $OutDir_species/Chr$chr.$sample.H1.fa
	cat $OutDir_fa/ref.Chr$chr.fa | vcf-consensus -H 2 -s $sample $inputvcf > $OutDir_species/Chr$chr.$sample.H2.fa
	####mask those low quality regions that were removed in the previous steps
	bedtools maskfasta -fi $OutDir_species/Chr$chr.$sample.H1.fa -bed $mask_bed -fo $OutDir_species/Chr$chr.$sample.H1.mask.fa	
	bedtools maskfasta -fi $OutDir_species/Chr$chr.$sample.H2.fa -bed $mask_bed -fo $OutDir_species/Chr$chr.$sample.H2.mask.fa	
	
	sed 's/'"Chr$chr"'/'"${sample}_H1"'/g' $OutDir_species/Chr$chr.$sample.H1.mask.fa >  $OutDir_species/Chr$chr.$sample.h1.temp && mv $OutDir_species/Chr$chr.$sample.h1.temp $OutDir_species/Chr$chr.$sample.H1.mask.fa
	sed 's/'"Chr$chr"'/'"${sample}_H2"'/g' $OutDir_species/Chr$chr.$sample.H2.mask.fa >  $OutDir_species/Chr$chr.$sample.h2.temp && mv  $OutDir_species/Chr$chr.$sample.h2.temp $OutDir_species/Chr$chr.$sample.H2.mask.fa
done

cat $OutDir_species/Chr$chr.*.mask.fa > $OutDir_species/Chr$chr.$species.mask.fa
#####################################################################################
######Step3:use LDhelmet to estimate the recombination rate

elif [ $step == "3" ]; then

#############################################
###step1: build haplotype configuration list

input=
for chromo in {01..19}
do
input="$input $OutDir_species/Chr$chromo.$species.mask.fa"
done
echo $input

ldhelmet find_confs --num_threads 4 -w 50 -o $OutDir_out/$species.output.conf $input


#############################################
###step2: generating a likelihood table, different population-scaled mutation rates applied to trichcocarpa and other three asepn species

if [ $species == "trichocarpa" ]; then
ldhelmet table_gen --num_threads 1 -c $OutDir_out/$species.output.conf -t 0.005 -r 0.0 0.1 10.0 1.0 100.0 -o $OutDir_out/$species.output.lk
else
ldhelmet table_gen --num_threads 1 -c $OutDir_out/$species.output.conf -t 0.01 -r 0.0 0.1 10.0 1.0 100.0 -o $OutDir_out/$species.output.lk
fi


#############################################
###step3: generating a Pade coefficients table, different population-scaled mutation rates applied to trichcocarpa and other three asepn species

if [ $species == "trichocarpa" ]; then
ldhelmet pade --num_threads 1 -c $OutDir_out/$species.output.conf -t 0.005 -x 11 -o $OutDir_out/$species.output.pade
else
ldhelmet pade --num_threads 1 -c $OutDir_out/$species.output.conf -t 0.01 -x 11 -o $OutDir_out/$species.output.pade
fi


#############################################
elif [ $step == "4" ]; then
###step4: performs rjMCMC to estimate the recombination rate map 

ldhelmet rjmcmc --num_threads 1 -l $OutDir_out/$species.output.lk -p $OutDir_out/$species.output.pade -w 50 -b 50 -s $OutDir_species/Chr$chr.$species.mask.fa --burn_in 100000 -n 1000000 -o $OutDir_out/$species.Chr$chr.post


#############################################
###step5: take the output from rjmcmc and generates a text file with the results of the rjMCMC sampling procedure
ldhelmet post_to_text -m -p 0.025 -p 0.50 -p 0.0975 -o $OutDir_out/$species.Chr$chr.txt $OutDir_out/$species.Chr$chr.post

#####################################################################################
######Step4:summarize the output of the Ldhelmet

elif [ $step == "5" ]; then

input=
for chromo in {01..19}
do
sed '1,3d' $OutDir_out/$species.Chr$chromo.txt |awk '{print "Chr'$chromo'",$0}' > $OutDir_out/$species.Chr$chromo.temp
input="$input $OutDir_out/$species.Chr$chromo.temp"
done
cat $input > $OutDir_out/$species.ldhelmet.txt && rm $input

fi



