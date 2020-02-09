#! /bin/bash -l


fsc26="/proj/sllstore2017050/nobackup/milou_files/tools/fsc26_linux64/fsc26"
model="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/demographic_inference/fastsimcoal2/output/ANGSD/model"
directory=`pwd`

#########################Step1##########################################

if [ ! -d "$model" ]; then
mkdir -p $model
fi

for n in {1..50}
do

run=$model/run$n

if [ ! -d "$run" ]; then
mkdir -p $run
fi

$fsc26 -t aspen.3species.model.tpl -n100000 -d -e aspen.3species.model.est -l 10 -L 40 -q -c 4 -M 1e-5 -w 1e-5 -B 32 -N100000

`cp -r $directory/aspen.3species.model/ $run`

done


#############################Step2########################################
##summarize the results into a single summary file

summary_dir=$model2/summary

if [ ! -d "$summary_dir" ]; then
mkdir -p $summary_dir
fi

cat $model/run1/aspen.3species.model/aspen.3species.model.bestlhoods > $summary_dir/model.summary.txt

for n in {2..50}
do
sed '1d' $model/run$n/aspen.3species.model/aspen.3species.model.bestlhoods >> $summary_dir/model.summary.txt
done




