#! /bin/bash -l

angsd="/proj/b2011141/tools/angsd_Github_dev/angsd/angsd"
estAvg="/proj/b2011141/tools/angsd_Github_dev/angsd/R/estAvgError.R"

bam="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/abbababa/4species.bam.list"
OutDir_1="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/abbababa"
region="/proj/b2011141/nobackup/PaperIV-phylogenomics/bed/chr"
window=$1  ###the window size, e.g. 10000,100000 
step=$2 ##which step to take the performance
chromosome=$3   ###the chromosomes to perform the analysis, e.g. 01,02,03...

 
OutDir=$OutDir_1/$window

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

if [ "$step" == "1" ]; then

#file=$chromosome
#chr=$region/all.filter.Chr$file.region
#$angsd -doAbbababa2 1 -bam $bam -doCounts 1 -blockSize $window -sample 0 -enhance 0 -out $OutDir/4species.abbababa.chr$file.$window.out -sizeH1 8 -sizeH2 8 -sizeH3 8 -sizeH4 8 -combFile 1 -rf $chr -minMapQ 30 -minQ 20

for file in {02..19}
do
abbababa="$OutDir/4species.abbababa.chr$file.$window.out.abbababa2"
sed '1d' $abbababa > $OutDir/temp && mv $OutDir/temp $abbababa
done

all=
for file in $OutDir/*abbababa2; do
all="$all $file"; done
cat $all > $OutDir/4species.abbababa.all.$window.out.abbababa2

Rscript $estAvg angsdFile=$OutDir/4species.abbababa.all.$window.out out=$OutDir/4species.outDstat


elif [ "$step" == "2" ]; then

#if [ "$chromosome" == "01"]; then
#Rscript $estAvg angsdFile=$OutDir/4species.abbababa.chr$chromosome.$window.out out=$OutDir/4species.$chromosome.$window.outDstat
#else
#cat $OutDir/header $OutDir/4species.abbababa.chr$chromosome.$window.out.abbababa2 > $OutDir/4species.abbababa.chr$chromosome.$window.out.temp && mv $OutDir/4species.abbababa.chr$chromosome.$window.out.temp $OutDir/4species.abbababa.chr$chromosome.$window.out.abbababa2
#Rscript $estAvg angsdFile=$OutDir/4species.abbababa.chr$chromosome.$window.out out=$OutDir/4species.$chromosome.$window.outDstat
#fi

Rscript $estAvg angsdFile=$OutDir/4species.abbababa.chr$chromosome.$window.out out=$OutDir/4species.$chromosome.$window.outDstat

elif [ "$step" == "3" ]; then

#####################################################
###Make D statistics and ABBA-BABA summary table for each window
nInd=$(cat $OutDir/4species.abbababa.all.$window.out.abbababa2 | wc -l)

all_obs=$OutDir/4species.all.outDstat.Observed.txt
all_rem=$OutDir/4species.all.outDstat.RemTrans.txt

for (( i=2; i<=$nInd; i++ ))
do
sed -n ${i}p $OutDir/4species.abbababa.all.$window.out.abbababa2 > $OutDir/temp
cat $OutDir/header $OutDir/temp > $OutDir/temp.abbababa2

Rscript $estAvg angsdFile=$OutDir/temp out=$OutDir/4species.$i.outDstat
sed '1d' $OutDir/4species.$i.outDstat.Observed.txt >> $all_obs
sed '1d' $OutDir/4species.$i.outDstat.RemTrans.txt >> $all_rem
rm $OutDir/temp $OutDir/temp.abbababa2 $OutDir/4species.$i.outDstat.Observed.txt $OutDir/4species.$i.outDstat.RemTrans.txt
done

#####################################################
###Summarize the results
head -n 1 $OutDir/4species.outDstat.Observed.txt > $OutDir/out.header
cat $OutDir/out.header $all_obs > $OutDir/temp && mv $OutDir/temp $all_obs
cat $OutDir/out.header $all_rem > $OutDir/temp && mv $OutDir/temp $all_rem

cut -f 1-3 $OutDir/4species.abbababa.all.$window.out.abbababa2 > $OutDir/chro
paste $OutDir/chro $all_obs > $OutDir/temp && mv $OutDir/temp $all_obs
paste $OutDir/chro $all_rem > $OutDir/temp && mv $OutDir/temp $all_rem

rm $OutDir/out.header $OutDir/chro
fi

