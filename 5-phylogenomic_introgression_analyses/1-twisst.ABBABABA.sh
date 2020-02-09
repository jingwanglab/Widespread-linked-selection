#! /bin/bash -l

##Main aim: This script is mainly used to caculate topology weightying by iterative sampling of sub-trees (Twisst) developed by Simon Martin. There are several steps needed to be performed in order to finally calculate Twisst weightings, e.g. BEAGLE imputation, python script transfer .vcf to .geno, make trees for sliding windows using Phyml and RAxML......

module load bioinfo-tools
module load python/2.7.11
module load java
module load raxml/8.2.4-gcc

step=$1   ####the step to perform the analysis
window=$2   ###the windows to perform the analysis: 10000,100000

###Step1: Use BEAGLE 4.1 to phase the SNP data for the four aspen species: P. tremula, P. davidana, P. tremuloides, P. trichocarpa

Inputvcf="/proj/b2011141/nobackup/PaperIV-phylogenomics/GATK/HC/total3/beagle/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.vcf.gz"
Out=${Inputvcf##*/}
echo $Out

vcf_Dir=`dirname $Inputvcf`

twisst_dir=$vcf_Dir/twisst

if [ ! -d "$twisst_dir" ]; then
mkdir -p $twisst_dir
fi
###Step2:Parsing VCF files, use python script to transfer vcf file into .geno file, with both un-phased data and phased data

if [ "$step" == "1" ]; then
parse_vcf="/proj/b2011141/tools/twisst/genomics_general-master/VCF_processing/parseVCF.py"

python $parse_vcf -i $Inputvcf |gzip > $twisst_dir/${Out%.vcf.gz}.geno.gz

elif [ "$step" == "2" ]; then
###Step3: Diversity and divergence analyses in sliding windows
popgen_window="/proj/b2011141/tools/twisst/genomics_general-master/popgenWindows.py"

##10kb
python $popgen_window -w 10000 -m 50 -g $twisst_dir/${Out%.vcf.gz}.geno.gz -o $twisst_dir/${Out%.vcf.gz}.10kb.popgen.csv.gz -f phased -p tremula SwAsp001,SwAsp014,SwAsp021,SwAsp033,SwAsp087,SwAsp096,SwAsp110,SwAsp114 -p davidiana TN1410D3187,TN1410D3189,TN1410D3190,TN1410D3191,TN1410D3192,TN1410D3194,TN1410D3590,TN1410D3591 -p tremuloides Alb10-3,Alb13-1,Alb18-3,Alb27-1,Alb31-1,Alb33-2,Alb6-3,Albb15-3 -p trichocarpa SRR1569629,SRR1569814,SRR1570762,SRR1571152,SRR1571263,SRR1571362,SRR1571416,SRR1571500 -T 2
##100kb
python $popgen_window -w 100000 -m 200 -g $twisst_dir/${Out%.vcf.gz}.geno.gz -o $twisst_dir/${Out%.vcf.gz}.100kb.popgen.csv.gz -f phased -p tremula SwAsp001,SwAsp014,SwAsp021,SwAsp033,SwAsp087,SwAsp096,SwAsp110,SwAsp114 -p davidiana TN1410D3187,TN1410D3189,TN1410D3190,TN1410D3191,TN1410D3192,TN1410D3194,TN1410D3590,TN1410D3591 -p tremuloides Alb10-3,Alb13-1,Alb18-3,Alb27-1,Alb31-1,Alb33-2,Alb6-3,Albb15-3 -p trichocarpa SRR1569629,SRR1569814,SRR1570762,SRR1571152,SRR1571263,SRR1571362,SRR1571416,SRR1571500 -T 2

elif [ "$step" == "3" ];then
###Step4: Compute ABBA-BABA statistics in sliding windows
abbababa="/proj/b2011141/tools/twisst/genomics_general-master/ABBABABAwindows.py"

if [ "$window" == "10000" ];then
#10kb
python $abbababa -w 10000 -m 50 -g $twisst_dir/${Out%.vcf.gz}.geno.gz -o $twisst_dir/${Out%.vcf.gz}.10kb.abbababa.csv -f phased -P1 tremula SwAsp001,SwAsp014,SwAsp021,SwAsp033,SwAsp087,SwAsp096,SwAsp110,SwAsp114 -P2 davidiana TN1410D3187,TN1410D3189,TN1410D3190,TN1410D3191,TN1410D3192,TN1410D3194,TN1410D3590,TN1410D3591 -P3 tremuloides Alb10-3,Alb13-1,Alb18-3,Alb27-1,Alb31-1,Alb33-2,Alb6-3,Albb15-3 -O trichocarpa SRR1569629,SRR1569814,SRR1570762,SRR1571152,SRR1571263,SRR1571362,SRR1571416,SRR1571500 -T 6
elif [ "$window" == "100000" ]; then
##100kb
python $abbababa -w 100000 -m 200 -g $twisst_dir/${Out%.vcf.gz}.geno.gz -o $twisst_dir/${Out%.vcf.gz}.100kb.abbababa.csv -f phased -P1 tremula SwAsp001,SwAsp014,SwAsp021,SwAsp033,SwAsp087,SwAsp096,SwAsp110,SwAsp114 -P2 davidiana TN1410D3187,TN1410D3189,TN1410D3190,TN1410D3191,TN1410D3192,TN1410D3194,TN1410D3590,TN1410D3591 -P3 tremuloides Alb10-3,Alb13-1,Alb18-3,Alb27-1,Alb31-1,Alb33-2,Alb6-3,Albb15-3 -O trichocarpa SRR1569629,SRR1569814,SRR1570762,SRR1571152,SRR1571263,SRR1571362,SRR1571416,SRR1571500 -T 6
fi

elif [ "$step" == "4" ];then
###Step5: Using raxml to make phylogenetic trees

raxml_tool="/proj/b2011141/tools/twisst/genomics_general-master/phylo/raxml_sliding_windows.py"
raxml="/proj/b2011141/tools/twisst/genomics_general-master/phylo/raxml"

###reconstruct the tree in 10kb or 100kb window size
if [ "$window" == "10000" ];then
##10kb 
python $raxml_tool -g $twisst_dir/${Out%.vcf.gz}.geno.gz --prefix $twisst_dir/${Out%%.vcf.gz}.raxml.w10kb -w 10000 --windType coordinate --model GTRCATI --genoFormat phased --raxml $raxml --outgroup SRR1569629,SRR1569814,SRR1570762,SRR1571152,SRR1571263,SRR1571362,SRR1571416,SRR1571500 -T 16 -M 50

elif [ "$window" == "100000" ]; then
##100kb
python $raxml_tool -g $twisst_dir/${Out%.vcf.gz}.geno.gz --prefix $twisst_dir/${Out%%.vcf.gz}.raxml.w100kb -w 100000 --windType coordinate --model GTRCATI --genoFormat phased --raxml $raxml --outgroup SRR1569629,SRR1569814,SRR1570762,SRR1571152,SRR1571263,SRR1571362,SRR1571416,SRR1571500 -T 16 -M 200
fi


elif [ "$step" == "5" ];then
###Step6: Using twisst to calculate the phylogenetic tree weightings
twisst_tool="/proj/b2011141/tools/twisst/twisst-master/twisst.py"
group="/proj/b2011141/nobackup/PaperIV-phylogenomics/GATK/HC/total3/groups.hap.tsv"
export PATH=~/anaconda_ete/bin:$PATH

###NOTICE!  I just realized that some of the windows do not contain SNPs, so those corresponding lines are shown by "NA", I therefore use the following command to delete those lines with SNPs number <50 for 10kb windows and <200 for 100kb windows

gunzip $twisst_dir/${Out%%.vcf.gz}.raxml.w10kb.trees.gz
sed '1d' $twisst_dir/${Out%%.vcf.gz}.raxml.w10kb.data.tsv > $twisst_dir/raxml.w10kb.data
paste $twisst_dir/raxml.w10kb.data $twisst_dir/${Out%%.vcf.gz}.raxml.w10kb.trees |awk '$5>50' |cut -f 6- |gzip - > $twisst_dir/${Out%%.vcf.gz}.raxml.w10kb.shrink.trees.gz
gzip $twisst_dir/${Out%%.vcf.gz}.raxml.w10kb.trees 
rm $twisst_dir/raxml.w10kb.data

weighting_dir=$vcf_Dir/twisst/weightings

if [ ! -d "$weighting_dir" ]; then
mkdir -p $weighting_dir
fi

cp $twisst_dir/${Out%%.vcf.gz}.raxml.w10kb.shrink.trees.gz $weighting_dir
#10kb
python $twisst_tool -t $weighting_dir/${Out%%.vcf.gz}.raxml.w10kb.shrink.trees.gz -w $weighting_dir/${Out%%.vcf.gz}.twisst.weights.10kb.csv.gz -o $weighting_dir/${Out%%.vcf.gz}.twisst.topologies.10kb.trees -g P.tremula -g P.tremuloides -g P.davidiana -g P.trichocarpa --method complete --groupsFile $group -D $weighting_dir/${Out%%.vcf.gz}.twisst.distance.10kb.csv.gz

#100kb

gunzip $twisst_dir/${Out%%.vcf.gz}.raxml.w100kb.trees.gz
sed '1d' $twisst_dir/${Out%%.vcf.gz}.raxml.w100kb.data.tsv > $twisst_dir/raxml.w100kb.data
paste $twisst_dir/raxml.w100kb.data $twisst_dir/${Out%%.vcf.gz}.raxml.w100kb.trees |awk '$5>200' |cut -f 6- |gzip - > $twisst_dir/${Out%%.vcf.gz}.raxml.w100kb.shrink.trees.gz
gzip $twisst_dir/${Out%%.vcf.gz}.raxml.w100kb.trees
rm $twisst_dir/raxml.w100kb.data

cp $twisst_dir/${Out%%.vcf.gz}.raxml.w100kb.shrink.trees.gz $weighting_dir
python $twisst_tool -t $weighting_dir/${Out%%.vcf.gz}.raxml.w100kb.shrink.trees.gz -w $weighting_dir/${Out%%.vcf.gz}.twisst.weights.100kb.csv.gz -o $weighting_dir/${Out%%.vcf.gz}.twisst.topologies.100kb.trees -g P.tremula -g P.tremuloides -g P.davidiana -g P.trichocarpa --method complete --groupsFile $group -D $weighting_dir/${Out%%.vcf.gz}.twisst.distance.100kb.csv.gz

fi


