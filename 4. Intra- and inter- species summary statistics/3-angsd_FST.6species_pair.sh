#! /bin/bash -l

angsd="/proj/b2011141/tools/angsd0.917/angsd"
emOptim2="/proj/b2011141/tools/angsd0.917/misc/emOptim2"
thetaStat="/proj/b2011141/tools/angsd0.917/misc/thetaStat"
ngsFst="/proj/b2011141/tools/ngsTools/ngsPopGen/ngsFST"


davidiana_out_saf="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/davidiana/davidiana_$1/davidiana_$1.rf.fix.saf"
tremula_out_saf="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/tremula/tremula_$1/tremula_$1.rf.fix.saf"
tremuloides_out_saf="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/tremuloides/tremuloides_$1/tremuloides_$1.rf.fix.saf"
trichocarpa_out_saf="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/trichocarpa/trichocarpa_$1/trichocarpa_$1.rf.fix.saf"

davidiana_sfs="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/davidiana/davidiana_$1/davidiana_$1.sfs"
tremula_sfs="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/tremula/tremula_$1/tremula_$1.sfs"
tremuloides_sfs="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/tremuloides/tremuloides_$1/tremuloides_$1.sfs"
trichocarpa_sfs="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/trichocarpa/trichocarpa_$1/trichocarpa_$1.sfs"

final_fix_pos="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/davidiana/davidiana_$1/4species.intersect.$1.rf.pos"

tremula_tremuloides="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/ngstools/fst/tremula_tremuloides"
davidiana_trichocarpa="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/ngstools/fst/davidiana_trichocarpa"
tremula_davidiana="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/ngstools/fst/tremula_davidiana"
tremula_trichocarpa="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/ngstools/fst/tremula_trichocarpa"
tremuloides_davidiana="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/ngstools/fst/tremuloides_davidiana"
tremuloides_trichocarpa="/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/ngstools/fst/tremuloides_trichocarpa"


ngsFst="/proj/b2011141/tools/ngsTools/ngsPopGen/ngsFST"


nsites=$(cat $final_fix_pos | wc -l)

#using marginal spectra as priors to calculate Fst between populations, estimated using optimSFS (ANGSD)
$ngsFst -postfiles $tremula_out_saf $tremuloides_out_saf -priorfiles $tremula_sfs $tremuloides_sfs -nind 8 8 -outfile $tremula_tremuloides/tremula_tremuloides.chr$1.fst -nsites $nsites -block_size 20000 -islog 1
$ngsFst -postfiles $davidiana_out_saf $trichocarpa_out_saf -priorfiles $davidiana_sfs $trichocarpa_sfs -nind 8 8 -outfile $davidiana_trichocarpa/davidiana_trichocarpa.chr$1.fst -nsites $nsites -block_size 20000 -islog 1
$ngsFst -postfiles $tremula_out_saf $davidiana_out_saf -priorfiles $tremula_sfs $davidiana_sfs -nind 8 8 -outfile $tremula_davidiana/tremula_davidiana.chr$1.fst -nsites $nsites -block_size 20000 -islog 1
$ngsFst -postfiles $tremula_out_saf $trichocarpa_out_saf -priorfiles $tremula_sfs $trichocarpa_sfs -nind 8 8 -outfile $tremula_trichocarpa/tremula_trichocarpa.chr$1.fst -nsites $nsites -block_size 20000 -islog 1
$ngsFst -postfiles $tremuloides_out_saf $davidiana_out_saf -priorfiles $tremuloides_sfs $davidiana_sfs -nind 8 8 -outfile $tremuloides_davidiana/tremuloides_davidiana.chr$1.fst -nsites $nsites -block_size 20000 -islog 1
$ngsFst -postfiles $tremuloides_out_saf $trichocarpa_out_saf -priorfiles $tremuloides_sfs $trichocarpa_sfs -nind 8 8 -outfile $tremuloides_trichocarpa/tremuloides_trichocarpa.chr$1.fst -nsites $nsites -block_size 20000 -islog 1


