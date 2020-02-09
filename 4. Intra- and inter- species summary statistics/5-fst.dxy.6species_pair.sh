#! /bin/bash -l

###Step1
####After estimation of fst and dxy from ANGSD and NGStools, add the Chr and Pos for all fst and dxy values in order for downstream sliding windows

##Rscript /proj/b2011141/pipeline/R/davidiana_paper/angsd/4species.fix.pos.R 01

###Step2
####Combine dxy and fst together 
species_pair=$1   ###tremula_tremuloides,tremula_davidiana,tremula_trichocarpa,tremuloides_davidiana,tremuloides_trichocarpa,davidiana_trichocarpa

all=
for chr in {01..19}
do

fst="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/ngstools/fst/$species_pair/$species_pair.chr$chr.fst"
dxy="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/ngstools/dxy/$species_pair/$species_pair.chr$chr.stat"
chr_pos="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/davidiana/davidiana_$chr/4species.intersect.$chr.rf.saf.fix.pos"

OutDir="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/ngstools/fst_dxy/$species_pair"

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

cut -f 1,2,3,4,5 $fst > $OutDir/chr$chr.fst
cut -f 7,8 $dxy > $OutDir/chr$chr.stat

paste $chr_pos $OutDir/chr$chr.fst $OutDir/chr$chr.stat > $OutDir/$species_pair.chr$chr.fst.stat && rm $OutDir/chr$chr.fst $OutDir/chr$chr.stat

all="$all $OutDir/$species_pair.chr$chr.fst.stat"
done

cat $all > $OutDir/$species_pair.all.fst_stat && rm $all

script_r="/proj/b2011141/pipeline/R/davidiana_paper/angsd/angsd_fst_stat_win.noChr.R"

Rscript $script_r $OutDir $species_pair 100000 100000
Rscript $script_r $OutDir $species_pair 10000 10000
Rscript $script_r $OutDir $species_pair 50000 50000

for file in {02..19}
do
sed '1d' $OutDir/$species_pair.win.fst.stat.Chr$file.window1e+05_step1e+05.fst.dxy.txt > $OutDir/temp && mv $OutDir/temp $OutDir/$species_pair.win.fst.stat.Chr$file.window1e+05_step1e+05.fst.dxy.txt
sed '1d' $OutDir/$species_pair.win.fst.stat.Chr$file.window10000_step10000.fst.dxy.txt > $OutDir/temp && mv $OutDir/temp $OutDir/$species_pair.win.fst.stat.Chr$file.window10000_step10000.fst.dxy.txt
sed '1d' $OutDir/$species_pair.win.fst.stat.Chr$file.window50000_step50000.fst.dxy.txt > $OutDir/temp && mv $OutDir/temp $OutDir/$species_pair.win.fst.stat.Chr$file.window50000_step50000.fst.dxy.txt
done

all=
all_1=
all_2=
for file in {01..19}
do
all="$all $OutDir/$species_pair.win.fst.stat.Chr$file.window1e+05_step1e+05.fst.dxy.txt"
all_1="$all_1 $OutDir/$species_pair.win.fst.stat.Chr$file.window10000_step10000.fst.dxy.txt"
all_2="$all_2 $OutDir/$species_pair.win.fst.stat.Chr$file.window50000_step50000.fst.dxy.txt"
done
cat $all > $OutDir/$species_pair.all.window1e+05_step1e+05.fst.dxy.txt && rm $all
cat $all_1 > $OutDir/$species_pair.all.window10000_step10000.fst.dxy.txt && rm $all_1
cat $all_2 > $OutDir/$species_pair.all.window50000_step50000.fst.dxy.txt && rm $all_2


