#! /bin/bash -l

module load bioinfo-tools
module load BEDTools
module load samtools
module load python3

##########################################################################################################
#The main aim of this method is to use the MSMC2 to infer demography history of the four species
#As before, I also use 2,4,8 haplotypes to perform the analysis
#for 2 hapltoypes, we have 8 individuals
#for 4 haplotypes, we have 1-2,1-3,1-4,1-5,1-6,1-7,1-8,2-3,....,7-8, a total of 28 combinations
#for 8 haplotyeps, we have 1-2-3-4,1-2-3-5,1-2-3-6,1-2-3-7,1-2-3-8,1-2-4-5,1-2-4-6,1-2-4-7,1-2-4-8,1-2-5-6,1-2-5-7,1-2-5-8,1-2-6-7,1-2-6-8,1-3-4-5,1-3-4-6,...,5-6-7-8, a total of 70 combinations


tremula_ind="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/samples/tremula.inds"
tremuloides_ind="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/samples/tremuloides.inds"
trichocarpa_ind="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/samples/trichocarpa.inds"
davidiana_ind="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/samples/davidiana.inds"
group=$1 ###tremula, tremuloides, trichocarpa or davidiana
hap=$2  ###hap2,hap4,hap8

##########################################################################################################
#######OutDir#################################
Out="/proj/snic2016-7-89/nobackup/phylogenomic_paper/msmc2/out"

OutDir=$Out/$group/$hap

if [ ! -d "$OutDir" ]; then
mkdir -p $OutDir
fi

OutDir_analysis=$OutDir/analysis

if [ ! -d "$OutDir_analysis" ]; then
mkdir -p $OutDir_analysis
fi



group_sample=${group}_ind
nInd=$(cat ${!group_sample} | wc -l)

#########for only 2 haplotypes
if [ "$hap" == "hap2" ]; then

for i in $(seq 1 $nInd); do sample=$(head -n "$i" ${!group_sample} |tail -n 1 );
sbatch msmc2.4species.run.sh $group $hap $sample
done

fi

####for 4 haplotypes
###############because some files failed because of memory issue, so here I first tested whetehr ".final.txt" already exits or not 

if [ "$hap" == "hap4" ]; then

nInd1=$(echo "$nInd-1"| bc )
for i in $(seq 1 $nInd1); do
sample1=$(head -n "$i" ${!group_sample} |tail -n 1 );
for j in $(seq 2 $nInd); do
if [ $j -le $i ]
then
 continue
else
sample2=$(head -n "$j" ${!group_sample} |tail -n 1 );
echo "${sample1}_${sample2}"
###########test whether the .final.txt already exits or not
final="$OutDir_analysis/${sample1}_${sample2}.final.txt"

if [ -f $final ] 
then
break
else 
sbatch msmc2.4species.run.sh $group $hap ${sample1}_${sample2}
fi
fi
done
done
fi

####for 8 haplotypes

if [ "$hap" == "hap8" ]; then

nInd1=$(echo "$nInd-3"| bc )
nInd2=$(echo "$nInd-2"| bc )
nInd3=$(echo "$nInd-1"| bc )
for a in $(seq 1 $nInd1); do
sample1=$(head -n "$a" ${!group_sample} |tail -n 1 );
for b in $(seq 2 $nInd2); do
        if [ $b -le $a ]
        then
        continue
        else
        sample2=$(head -n "$b" ${!group_sample} |tail -n 1 );
        for c in $(seq 3 $nInd3); do
                if [ $c -le $b ]
                then
                continue
                else
                sample3=$(head -n "$c" ${!group_sample} |tail -n 1 );
                for d in $(seq 4 $nInd); do
                        if [ $d -le $c ]
                        then
                        continue
                        else
                        sample4=$(head -n "$d" ${!group_sample} |tail -n 1 );
                        echo "${sample1}_${sample2}_${sample3}_${sample4}"
#			test whether the .final.txt already exits or not
			final="$OutDir_analysis/${sample1}_${sample2}_${sample3}_${sample4}.final.txt"
			if [ -f $final ]
			then
			echo "already exits"
			continue
			else
			echo "not exits"
			sbatch msmc2.4species.run.sh $group $hap ${sample1}_${sample2}_${sample3}_${sample4}
			fi
fi
done
fi
done
fi
done
done
fi

