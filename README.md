# FischleLab_KarthikEswara_ribo0seq_repeats

## Preparing gtf

Nonredunant genome annotations for repeats were obtained fom [Dfam database version 3.5](https://www.dfam.org/releases/Dfam_3.5/annotations/ce10/) and lifted over to the ce11 version of the genome. Repeat loci were given unique names by appending "_rpt#" to the repeat family name, and a gtf was constructed with the same annotation given gene, transcript and exon types. This was then concatenated with canonical gene annotation gtf downloaded from [Wormbase (WS298 version)](https://downloads.wormbase.org/releases/WS298/species/c_elegans/PRJNA13758/). This was performed with the `getRepeatData.R` script.

Repeat classification was downloaded with [FamDB tools](https://github.com/Dfam-consortium/FamDB) for Dfam version 4.0, and then merged with Dfam 3.5 annotations using the `getRepeatAnnotation.R`script.

## Downloading RNAseq from ENA

RNAseq samples were downloaded from ENA usig a list of ids recovered from the database

## Alignment with nf-core/rnaseq


## Final plots

**upreguledGenesOnArms.R** - boxplots, violin plots and heatmaps of expression of genes on autosomal arms. 

**integrationWithChIP.R** - euler plot of autosomal arm genes and their overlap with ChIPseq peaks of HPL-2, LIN-61 and H3K9me2. Boxplot and violin plot of gene expression of those autosomal arm genes that overlap chipseq peaks of all three types. Also produces lists of gene names of significantly upregulated genes on autosomal arms (in general, or those overlapping all three chipseq peak types) for running GO analysis; and bed files of significantly upregulated genes on autosomal arms with triple peaks for deeptools.

**deeptools_sigGenes_v_ChIP.sh** - creates heatmaps with deeptools of HPL-2 LIN-61 and H3K9me2 ChIPseq signal (from GSE271919) over upregulated gene on autosomal arms with triple peaks.

**functions_finalFigures.R** - some functions used by plotting scripts.

**getCountTablesForGEO.R** - script to extract raw counts of gene expression for the different datasets to upload to GEO.

