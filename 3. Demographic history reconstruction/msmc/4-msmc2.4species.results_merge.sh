#! /bin/bash -l


####Main aim: the main aim of this script is to merge all .final.txt files from all individual-combinations into a single file for R to make plot 

group=$1 ###tremula, tremuloides, trichocarpa or davidiana
hap=$2  ###hap2,hap4,hap8

##########################################################################################################
#######OutDir of .final results#################################
Out="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/out"

OutDir=$Out/$group/$hap
OutDir_analysis=$OutDir/analysis

OutDir_summary=$OutDir_analysis/summary

if [ ! -d "$OutDir_summary" ]; then
mkdir -p $OutDir_summary
fi

for final in $OutDir_analysis/*.final.txt
do
input=${final##*/}
inds=${input%%\.final.txt}

sed '1d' $final | awk '{print $0,"'$inds'"}' > $OutDir_summary/$inds.final.sample.txt
done

echo -e "time_index\tleft_time_boundary\tright_time_boundary\tlambda\tsample" > $OutDir_summary/$group.$hap.msmc2.summary.txt

input=
for summary in $OutDir_summary/*final.sample.txt
do
input="$input $summary"
done

cat $input >> $OutDir_summary/$group.$hap.msmc2.summary.txt
rm $input

