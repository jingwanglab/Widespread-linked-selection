#! /bin/bash -l

module load bioinfo-tools
module load vcftools

Inputvcf=$1
VCFDir=`dirname $1`
MLtree=$VCFDir/ml_tree

if [ ! -d "$MLtree" ]; then
mkdir -p $MLtree
fi

vcf_to_fasta="/proj/b2011141/tools/vcf_to_fasta/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl"
muscle="/proj/b2011141/tools/snphylo/bin/muscle"
DNAML="/proj/b2011141/tools/snphylo/phylip/exe/dnaml"
SCRIPTS_DIR="/proj/b2011141/tools/snphylo/scripts"

Out=${Inputvcf##*/}
echo $Out
OutSuffix=${Out%.recode.vcf.gz}

out_sample_id=$2

##transfer vcf file to fasta file
zcat $Inputvcf | vcf-to-tab > $MLtree/$OutSuffix.tab
if [ -z "${out_sample_id}" ]
then
perl $vcf_to_fasta --output_ref -i $MLtree/$OutSuffix.tab > $MLtree/$OutSuffix.fasta
else
perl $vcf_to_fasta -i $MLtree/$OutSuffix.tab > $MLtree/$OutSuffix.fasta
fi

##muscle alignment and transfer fasta file to phylip file
$muscle -phyi -in $MLtree/$OutSuffix.fasta -out $MLtree/$OutSuffix.phylip.txt
ln -sf $MLtree/$OutSuffix.phylip.txt $MLtree/infile

##add outgroup sample number
cd $VCFDir

if [ ! -z "${out_sample_id}" ]
then
    echo -e "y\n" | "${DNAML}"; [ $? != 0 ] && exit 1
else
    out_sample_no=$[$(grep -ne "\<${out_sample_id}\>" infile | cut -f1 -d':') - 1]
    if [ ${out_sample_no} -eq -1 ]
    then
        echo "Error!!! There is no sample name (${out_sample_id}) to use as a outgroup in a tree input file ($OutSuffix.phylip.txt)." 1>&2
        echo "Please check the name and restart this script." 1>&2
        rm -f infile
        exit 1
    else
        echo -e "o\n${out_sample_no}\ny\n" | "${DNAML}"; [ $? != 0 ] && exit 1
    fi
fi
mv outfile $OutSuffix.ml.txt
mv outtree $OutSuffix.ml.tree

rm -f infile

if [ -e $OutSuffix.ml.tree ]
then
    Rscript --slave --vanilla --file="${SCRIPTS_DIR}/draw_unrooted_tree.R" --args -i "${OutSuffix}.ml.tree" -o "${OutSuffix}"
    [ $? != 0 ] && exit 1
fi

