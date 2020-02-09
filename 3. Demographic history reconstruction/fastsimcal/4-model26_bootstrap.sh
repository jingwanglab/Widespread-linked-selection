#! /bin/bash -l


step=$1

fsc26="/proj/sllstore2017050/nobackup/milou_files/tools/fsc26_linux64/fsc26"

max_par="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/demographic_inference/fastsimcoal2/output/ANGSD/best_model26/aspen.3species.model26_maxL.par"


if [ $step == "1" ];then
#########Step1:generating the 170Mb DNA sequence data set and SFS, each run include 1.7Mb and add 100 runs of data into one single 2d-SFS file

###total number of sites: 168,950,389  ##approx. 170Mb
#for rep in {1..100}
#do
rep=$2
bootstrap="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/demographic_inference/fastsimcoal2/output/ANGSD/best_model26/bootstrap/rep$rep"

if [ ! -d "$bootstrap" ]; then
mkdir -p $bootstrap
fi

$fsc26 -i $max_par -n 100 -d -s0 -x -k10000000 -q 

sed '1,2d' aspen.3species.model26_maxL/aspen.3species.model26_maxL_jointDAFpop1_0.obs |cut -f 2- >  $bootstrap/pop1_0.obs
sed '1,2d' aspen.3species.model26_maxL/aspen.3species.model26_maxL_jointDAFpop2_0.obs |cut -f 2- >  $bootstrap/pop2_0.obs
sed '1,2d' aspen.3species.model26_maxL/aspen.3species.model26_maxL_jointDAFpop2_1.obs |cut -f 2- >  $bootstrap/pop2_1.obs
Rscript model26.bootstrap.jointSFS.R $bootstrap 

rm $bootstrap/pop1_0.obs $bootstrap/pop2_0.obs $bootstrap/pop2_1.obs
#done


#########Step2: run fastsimcoal on simulated datasets for each replication (100 in total), with 50 independent runs in each replication
elif [ $step == "2" ];then
rep=$2
bootstrap="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/demographic_inference/fastsimcoal2/output/ANGSD/best_model26/bootstrap/rep$rep"

if [ ! -d "$bootstrap" ]; then
mkdir -p $bootstrap
fi
 

cp aspen.3species.model26.tpl $bootstrap/aspen.3species.model26_maxL.tpl
cp aspen.3species.model26.est $bootstrap/aspen.3species.model26_maxL.est
cp aspen.3species.model26.pv $bootstrap/aspen.3species.model26_maxL.pv

cd $bootstrap 

for n in {1..50}
do
run=$bootstrap/run$n

if [ ! -d "$run" ]; then
mkdir -p $run
fi

$fsc26 -t aspen.3species.model26_maxL.tpl -n100000 -d -e aspen.3species.model26_maxL.est -l 10 -L 40 -q -c 4 -M 1e-5 -w 1e-5 -B 32 -N100000  --initValues  aspen.3species.model26_maxL.pv
`cp -r $bootstrap/aspen.3species.model26_maxL $run`

done


###Step3: for each replication (100 in total) for bootstrap estimates, only choose the one from the 50 independent runs with lowest MaxEstLhood
elif [ $step == "3" ];then
rep=$2
bootstrap="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/demographic_inference/fastsimcoal2/output/ANGSD/best_model26/bootstrap/rep$rep"

cat $bootstrap/run1/aspen.3species.model26_maxL/aspen.3species.model26_maxL.bestlhoods > $bootstrap/rep$rep.run.compare.bestlhoods


for n in {2..50}
do
run=$bootstrap/run$n
sed '1d' $run/aspen.3species.model26_maxL/aspen.3species.model26_maxL.bestlhoods >> $bootstrap/rep$rep.run.compare.bestlhoods
done

Rscript model.compare.R $bootstrap $rep

###Step4. summarize the 100 bootstrap best runs into a table for calculating 95% confidence internvals for the parameters of the best model
elif [ $step == "4" ];then

bootstrap="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/demographic_inference/fastsimcoal2/output/ANGSD/best_model26/bootstrap"

cat $bootstrap/rep1/rep1.best.run.bestlhoods > $bootstrap/model26.bootstrap.run1_100.bestlhoods

for rep in {2..100}
do
sed '1d' $bootstrap/rep$rep/rep$rep.best.run.bestlhoods >>  $bootstrap/model26.bootstrap.run1_100.bestlhoods
done


fi


