#! /bin/bash -l

module load bioinfo-tools
module load python

###Main aim: the main aim of the script is to detect ancient balancing selection for the 3 aspen species using the unfolded site frequency spectrum and also detect ancient balancing selection among the four Populus species using the folded site frequency spectrum
##For the three aspen species, only the polymorphic sites (compared to the ancestral allele) was included.

step=$1  ###the step to perfrom the analysis
species=$2   ###the species to perform the analysis
chr=$3   ###the chromsome to perform the analysis
window=$4  ###window size to perform beta analysis, e.g. 1000, 10000...


betascan="/proj/sllstore2017050/nobackup/milou_files/tools/BetaScan-master/BetaScan.py"

OutDir_unfolded="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/GATK/HC/total3/Betascan/$species/out/unfolded/window$window"
if [ ! -d "$OutDir_unfolded" ]; then
mkdir -p $OutDir_unfolded
fi

OutDir_folded="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/GATK/HC/total3/Betascan/$species/out/folded/window$window"
if [ ! -d "$OutDir_folded" ]; then
mkdir -p $OutDir_folded
fi

InputDir="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/GATK/HC/total3/Betascan/$species/input"
if [ ! -d "$InputDir" ]; then
mkdir -p $InputDir
fi

#####Step1: creating the input files for the betascan
if [ "$step" == "1" ]; then
###create the input file of Betascan from the input file from CalcABS
input="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/GATK/HC/total3/CalcABS/species.$chr.input.txt"

if [ "$species" == "tremula" ];then
cut -f 2,3,4 $input |awk '$2>0' |awk '$2!=$3' | sed '1d' |sed 's/ /\t/g' > $InputDir/$species.$chr.betascan

elif [ "$species" == "davidiana" ]; then
cut -f 2,5,6 $input |awk '$2>0' |awk '$2!=$3' | sed '1d' |sed 's/ /\t/g' > $InputDir/$species.$chr.betascan

elif [ "$species" == "tremuloides" ]; then
cut -f 2,7,8 $input |awk '$2>0' |awk '$2!=$3' | sed '1d' |sed 's/ /\t/g' > $InputDir/$species.$chr.betascan

elif [ "$species" == "trichocarpa" ]; then
cut -f 2,9,10 $input |awk '$2>0' |awk '$2!=$3' | sed '1d' |sed 's/ /\t/g' > $InputDir/$species.$chr.betascan

fi

#####Step2: run betascan use unfolded version
elif [ "$step" == "2" ]; then

###using unfolded version
python $betascan -i $InputDir/$species.$chr.betascan -w $window -m 0.2 
mv Betas_$species.$chr.betascan $OutDir_unfolded/

#####Step2.2: summarize the results from all 19 chromsomes into one file
elif [ "$step" == "2.2" ]; then

input=
for chr in Chr{01..19}
do
sed 's/ /\t/g' $OutDir_unfolded/Betas_$species.$chr.betascan > $OutDir_unfolded/Betas_$species.$chr.betascan.temp
awk '$1="'$chr'"' $OutDir_unfolded/Betas_$species.$chr.betascan |cut -f 1 -d " " > $OutDir_unfolded/Betas_$species.$chr.betascan.chr
paste $OutDir_unfolded/Betas_$species.$chr.betascan.chr $OutDir_unfolded/Betas_$species.$chr.betascan.temp > $OutDir_unfolded/Betas_$species.$chr.betascan.txt
input="$input $OutDir_unfolded/Betas_$species.$chr.betascan.txt"
rm $OutDir_unfolded/Betas_$species.$chr.betascan.temp $OutDir_unfolded/Betas_$species.$chr.betascan.chr $OutDir_unfolded/Betas_$species.$chr.betascan
done

echo -e "Chromo\tPos\tScore" > $OutDir_unfolded/Betas_$species.all.betascan.header
cat $OutDir_unfolded/Betas_$species.all.betascan.header $input > $OutDir_unfolded/Betas_$species.all.betascan.txt && rm $input

#####Step3: run betascan use folded version
elif [ "$step" == "3" ]; then

python $betascan -i $InputDir/$species.$chr.betascan -w $window -m 0.2 -fold
mv Betas_$species.$chr.betascan $OutDir_folded/

elif [ "$step" == "3.2" ]; then

input=
for chr in Chr{01..19}
do
sed 's/ /\t/g' $OutDir_folded/Betas_$species.$chr.betascan > $OutDir_folded/Betas_$species.$chr.betascan.temp
awk '$1="'$chr'"' $OutDir_folded/Betas_$species.$chr.betascan |cut -f 1 -d " " > $OutDir_folded/Betas_$species.$chr.betascan.chr
paste $OutDir_folded/Betas_$species.$chr.betascan.chr $OutDir_folded/Betas_$species.$chr.betascan.temp > $OutDir_folded/Betas_$species.$chr.betascan.txt
input="$input $OutDir_folded/Betas_$species.$chr.betascan.txt"
rm $OutDir_folded/Betas_$species.$chr.betascan.temp $OutDir_folded/Betas_$species.$chr.betascan.chr $OutDir_folded/Betas_$species.$chr.betascan
done

echo -e "Chromo\tPos\tScore" > $OutDir_folded/Betas_$species.all.betascan.header
cat $OutDir_folded/Betas_$species.all.betascan.header $input > $OutDir_folded/Betas_$species.all.betascan.txt && rm $input


fi


