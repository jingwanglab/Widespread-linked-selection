# Widespread-linked-selection
Scripts for Wang et al (2020) Evidence for widespread selection in shaping the genomic landscape during speciation of Populus.

<h2>Documentation of Scripts</h2>

<b> 1. Read alignment, SNP and genotype calling </b>
        1-runBWA_mem.sh - Use bwa-mem to do read alignment
        2-realign.sh - Use GATK to do read realignment near indels

3-picard_markduplicates.sh
Use Picard to correct for artifacts of PCR duplication

4-haplotypecaller.sh
Use HaplotypeCaller of GATK to do SNP calling

5-genotypeGVCFs.sh
Use GenotypeGVCFs of GATK to do genotype calling among multiple individuals
 
6-vcf_geno_dep.filtering.4species.sh
Use multiple commands of BEDtools and VCFtools to do SNP filtering



