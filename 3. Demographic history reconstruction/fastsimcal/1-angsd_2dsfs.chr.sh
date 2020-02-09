#! /bin/bash -l

module load bioinfo-tools
module load ANGSD

chr=$1
species=$2  ###species to perform the analysis

tremula_list="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/bam/tremula/tremula.bam.list"
tremuloides_list="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/bam/tremuloides/tremuloides.bam.list"
davidiana_list="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/bam/davidiana/davidiana.bam.list"
trichocarpa_list="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/bam/trichocarpa/trichocarpa.bam.list"

ref="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/reference/Ptrichocarpa_v3.0_210.fa"
chr_region="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/bed/chr/all.filter.Chr$chr.region"

OutDir="/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/demographic_inference/fastsimcoal2/input/angsd/chr"


if [ "$species" == "tremula" ]; then
angsd -GL 1 -b $tremula_list -anc $ref -rf $chr_region -P 2 -out $OutDir/tremula.Chr$chr -doSaf 1

elif [ "$species" == "tremuloides" ]; then
angsd -GL 1 -b $tremuloides_list -anc $ref -rf $chr_region -P 2 -out $OutDir/tremuloides.Chr$chr -doSaf 1

elif [ "$species" == "davidiana" ]; then
angsd -GL 1 -b $davidiana_list -anc $ref -rf $chr_region -P 2 -out $OutDir/davidiana.Chr$chr -doSaf 1
fi

#angsd -GL 1 -b $trichocarpa_list -anc $ref -rf $chr_region -P 10 -out $OutDir/trichocarpa -doSaf 1



##########Step2:create 2dsfs for each chromosome dataset
##2.1 sfs for each species
#realSFS $OutDir/$species.Chr$chr.saf.idx -P 4 > $OutDir/$species.Chr$chr.saf.sfs

#2.2 2dsfs for pairs of species
realSFS $OutDir/tremula.Chr$chr.saf.idx $OutDir/davidiana.Chr$chr.saf.idx -P 4 > $OutDir/tremula.davidiana.Chr$chr.2dsfs.sfs
realSFS $OutDir/tremula.Chr$chr.saf.idx $OutDir/tremuloides.Chr$chr.saf.idx -P 4 > $OutDir/tremula.tremuloides.Chr$chr.2dsfs.sfs
realSFS $OutDir/davidiana.Chr$chr.saf.idx $OutDir/tremuloides.Chr$chr.saf.idx -P 4 > $OutDir/davidiana.tremuloides.Chr$chr.2dsfs.sfs


#########Step3: summarize the output of ANGSD into the input file of fastsimcoal
Rscript angsd_2dsfs.matrix.R

