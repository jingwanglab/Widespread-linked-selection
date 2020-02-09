#! /bin/bash -l

angsd="/proj/b2011141/tools/angsd0.917/angsd"
emOptim2="/proj/b2011141/tools/angsd0.917/misc/emOptim2"
thetaStat="/proj/b2011141/tools/angsd0.917/misc/thetaStat"
bam_list_tremuloides="/proj/b2011141/nobackup/PaperIV-phylogenomics/bam/tremuloides/tremuloides.bam.list"
ref="/proj/b2011141/nobackup/reference/nisqV3/Ptrichocarpa_v3.0_210.fa"
OutDir="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/tremuloides/tremuloides_$1"
region="/proj/b2011141/nobackup/PaperIV-phylogenomics/bed/chr"
chr=$region/all.filter.Chr$1.region

nInd=$(cat $bam_list_tremuloides | wc -l)
nChrom=$(echo "2*$nInd" | bc)
#nChrom=$nInd

echo $nInd
echo $nChrom

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

###first generate .saf file
#$angsd -bam $bam_list_tremuloides -minMapQ 30 -minQ 20 -GL 1 -doSaf 1 -out $OutDir/tremuloides_$1 -anc $ref -rf $chr

###use emOptim2 to optimization
#$emOptim2 $OutDir/tremuloides_$1.saf $nChrom -maxIter 100 -P 4 > $OutDir/tremuloides_$1.sfs

###calculate thetas
#$angsd -bam $bam_list_tremuloides -out $OutDir/tremuloides_$1 -doThetas 1 -GL 1 -doSaf 1 -anc $ref -rf $chr -pest $OutDir/tremuloides_$1.sfs -minMapQ 30 -minQ 20

###calculate Tajimas
#$thetaStat make_bed $OutDir/tremuloides_$1.thetas.gz
#$thetaStat do_stat $OutDir/tremuloides_$1.thetas.gz -nChr $nChrom
#$thetaStat do_stat $OutDir/tremuloides_$1.thetas.gz -nChr $nChrom -win 100000 -step 100000 -outnames $OutDir/tremuloides_$1.thetas100000.window.gz
#$thetaStat do_stat $OutDir/tremuloides_$1.thetas.gz -nChr $nChrom -win 10000 -step 10000 -outnames $OutDir/tremuloides_$1.thetas10000.window.gz
#$thetaStat do_stat $OutDir/tremuloides_$1.thetas.gz -nChr $nChrom -win 50000 -step 50000 -outnames $OutDir/tremuloides_$1.thetas50000.window.gz

#calculate posterior probabilities of sample allele frequencies
#$angsd -bam $bam_list_tremuloides -GL 1 -doSaf 1 -anc $ref -rf $chr -minMapQ 30 -minQ 20 -pest $OutDir/tremuloides_$1.sfs -out $OutDir/tremuloides_$1.rf

#####use r script to summarize the thetas table

thetas_r="/proj/b2011141/pipeline/R/davidiana_paper/angsd/thetas_table.R"

#Rscript $thetas_r tremuloides $1 100000 100000
#Rscript $thetas_r tremuloides $1 10000 10000
#Rscript $thetas_r tremuloides $1 50000 50000

#####Summarize all chromsome-level output files to a single file

OutDir_new="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/tremuloides"
for chr in {02..19}
do
sed '1d' $OutDir_new/tremuloides_$chr/tremuloides_$chr.10000bp10000bp.thetas.txt > $OutDir_new/tremuloides_all/tremuloides_$chr.10000bp10000bp.thetas.txt
done
cp $OutDir_new/tremuloides_01/tremuloides_01.10000bp10000bp.thetas.txt  $OutDir_new/tremuloides_all/tremuloides_01.10000bp10000bp.thetas.txt

input=
for chr in {01..19}
do
input="$input $OutDir_new/tremuloides_all/tremuloides_$chr.10000bp10000bp.thetas.txt"
done
cat $input > $OutDir_new/tremuloides_all/tremuloides_all.10000bp10000bp.thetas.txt && rm $input




