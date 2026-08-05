#!/bin/bash
#SBATCH --time=1-17:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

genomeVer=WS298
genomeDir=/mnt/meister.data/publicData/genomes
genomeFile=${genomeDir}/${genomeVer}/c_elegans.PRJNA13758.${genomeVer}.genomic.fa.gz
gtfFile=${genomeDir}/dfam35/${genomeVer}_canonicalgenes_Dfam3.5nr_repeats.uniq_gene_id.gtf


WORK_DIR=/mnt/meister.data/FischleLab_KarthikEswara/ribo0seq
CONFIG_FILE=/mnt/meister.data/nf-core/unibe_izb.config


nextflow run nf-core/rnaseq -profile singularity -r 3.26.0 --input ${WORK_DIR}/samplesheet_all.csv --multiqc_title multiqc_rnaseq --outdir $WORK_DIR -c $CONFIG_FILE --fasta $genomeFile --gtf $gtfFile --extra_star_align_args '--outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --alignIntronMax 100000 --outFilterMismatchNoverLmax 0.04' --pseudo_aligner kallisto -resume


# Mapping parameter recommendations from Schwarz, Robert, Philipp Koch, Jeanne Wilbrandt, and Steve Hoffmann. 2022. “Locus-Specific Expression Analysis of Transposable Elements.” Briefings in Bioinformatics 23 (1): bbab417.
