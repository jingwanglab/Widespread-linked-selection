# Widespread-linked-selection
Scripts for Wang et al (2020) Evidence for widespread selection in shaping the genomic landscape during speciation of Populus. Molecular Ecology, accepted

<h2>Documentation of Scripts</h2>

<b> 1. Read alignment, SNP and genotype calling </b>

        1-runBWA_mem.sh - Use bwa-mem to do read alignment
        2-realign.sh - Use GATK to do read realignment near indels
        3-picard_markduplicates.sh - Use Picard to correct for artifacts of PCR duplication
        4-haplotypecaller.sh - Use HaplotypeCaller of GATK to do SNP calling
        5-genotypeGVCFs.sh - Use GenotypeGVCFs of GATK to do genotype calling among multiple individuals
        6-vcf_geno_dep.filtering.4species.sh - Use multiple commands of BEDtools and VCFtools to do SNP filtering

<b> 2. Phylogenetic relationships and population structure analysis </b>

        1-UnifiedGenotyper.chloroplast.sh - Use UnifiedGenotyper of GATK to call chloroplast SNPs
        2-chloroplast_vcf_filtering.sh - Use a set of VCFtools commands to filter chloroplast SNPs
        3-SNPhylo-chloroplast.sh - Use SNPhylo to construct a consensus tree based on chloroplast SNPs
        4-SNPhylo-chloroplast.tree.R - Use R to plot the chloroplast phylogenetic tree
        5-angsd_PCA.4species.sh - Use ANGSD to perform principal component analysis
        6-pca.plot.R - Use R to plot PCA results
        7-beagle_phasing.sh - Use BEAGLE to impute and phase SNPs and perform identity-by-descent block analysis
        8-beagle_ibd.plot.R - Use R to plot the IBD and haplotype sharing between individuals of species

<b> 3. Demographic history reconstruction </b>

     fastsimcal
        1-angsd_2dsfs.chr.sh - Use ANGSD to construct two-dimensional joint site frequency spectrum (2d-sfs) 
        2-angsd_2dsfs.matrix.R - Use R to transfer the output from angsd to the input format of fastsimcoal2
        3-aspen.3species.model26.sh - One example showing the script running the fastsimcoal for one model
        4-model26_bootstrap.sh - Script to running bootstrap analysis for the best model (model26) inferred
     msmc
        1-msmc2.4species.input.sh- Shell script to generate the input file format for msmc 
        2-msmc2.4species.run.sh - Script for running msmc on different combinations of haplotypes
        3-msmc2.4species.submit.sh - Because there are many individual configurations and we use this script to submit these jobs respectively
        4-msmc2.4species.results_merge.sh - Script to merging the msmc results for different haplotype configuration
        5-msmc2.plot.R  - Use R to plot the msmc results for the four Populus species

<b> 4. Intra- and inter- species summary statistics </b>
        
        1-angsd_SFS_davidiana.sh - Use ANGSD to estimate the genetic diversity at each site of P. davidiana
        1-angsd_SFS_tremula.sh - Use ANGSD to estimate the genetic diversity at each site of P. tremula
        1-angsd_SFS_tremuloides.sh - Use ANGSD to estimate the genetic diversity at each site of P. tremuloides
        1-angsd_SFS_trichocarpa.sh - Use ANGSD to estimate the genetic diversity at each site of P. trichocarpa
        2-thetas_table.R - Use R to summarize the genetic diversity along the windows
        3-angsd_FST.6species_pair.sh - Use ANGSD to estimate genetic divergence, Fst, at each site between pairs of species
        4-angsd_dxy.6species_pair.sh - Use ANGSD to estimate genetic divergence, dxy, at each site between pairs of species
        5-angsd_fst_stat_win.noChr.R - Use R to summarize the genetic divergence (Fst and dxy) along the windows
        5-fst.dxy.6species_pair.sh - Script to summarize the fst and dxy estimates from ANGSD to 10Kbp and 100Kbp non-overlapping windows
        6-ldhelmet.4species.sh - Use Ldhelmet to estimate the population-scaled recombination rate of the four species
        6-ldhelmet.window.summary.R - Use R to summarize the population-scaled recombination rate from Ldhelmet into windows
        7-diversity_divergence_rho.4species.10kb.R - Summary analysis of genetic diversity, divergence and population-scaled recombination rate in 10kbp windows
        7-Figure2_diversity_divergence_rho.4species.100kb.R - Summary analysis of genetic diversity, divergence and population-scaled recombination rate in 100kbp windows and also include scripts for plotting Figure 2
        8-linked_selection.correlation.R - Use R to integrate the analysis of genetic diversity, divergence, recombination rate, coding density, incomplete lineage sorting, phylogenetic topologies and introgression proportion (fd)
        9-Figure4.R - Use R script to make analysis and plot the figure4

<b> 5-phylogenomic_introgression_analyses </b>

        1-twisst.ABBABABA.sh - Use twisst and ABBABABA to perform topology weighting and gene flow analysis (D statistics and fd)
        2-ILS.allele_freq.sh - Calculate incomplete lineage sorting (ILS) from allele frequencies of the four species
        3-ILS_gene.sh - Use script to calculate the distance to genes, exons and coding sequences and make the output files
        3-vcf2bed.awk - Awk script to transfer the vcf file into bed file format
        4-ILS.R - Use R script to estimate the level of incomplete linage sorting and plot the level of ILS with increasing distance of coding density
        5-Figure3_FigureS10S11.R - Use R script to plot Figure3 and Figure S10 and Figure S11
        6-angsd_abbababa2.sh - Use the abbababa2 option in ANGSD to perform the D statistics

<b> 6-positive_balancing_selection </b>

        1-xpclr_convert_rho_to_cM.R - Use R script to convert recombination rate estimates from Ldhelmet to centi-Morgan (cM) for P.tremula
        1-interpolate_map_to_snp.R - Use R script to interpolate genetic distance between adjacent SNPs according to the centimorgan values
        1-xpclr.aspen_species.rec_cM.sh - pipeline to run XPCLR
        1-Figure5.R - R script to plot the Figure5
        1-manhanttan_plot.R - Function to make the manhanttan plot
        2-betascan.sh - Script to run BetaScan analysis
        2-Figure6.R - R script to plot the Figure6
