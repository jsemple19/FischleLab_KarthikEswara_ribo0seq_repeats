#!/bin/bash
#SBATCH --time=0-14:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"
# nf-core/differentialabundance 2.0.0 uses workflow outputs; pin to a stable
# Nextflow version where this is no longer a preview feature.
export NXF_VER=25.10.4

genomeVer=WS298
genomeDir=/mnt/meister.data/publicData/genomes/${genomeVer}
#gtfFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.canonical_geneset.gtf
gtfFile=${genomeDir}/dfam35/${genomeVer}_canonicalgenes_Dfam3.5nr_repeats.uniq_gene_id.gtf

WORK_DIR=/mnt/meister.data/FischleLab_KarthikEswara/ribo0seq
CONFIG_FILE=/mnt/meister.data/nf-core/unibe_izb.config
#SRR_FILE=${WORK_DIR}/SRR_file.csv
lfcShrink=true
minAbund=5
minSamples=3

suffix="noShrink"
$lfcShrink && suffix="lfcShrink"

ALIGNER="kallisto" # should be one of kallisto, salmon or star
case "$ALIGNER" in
    kallisto)
        TABLE_PATH="kallisto/kallisto"
        ;;
    salmon)
        TABLE_PATH="salmon/salmon"
        ;;
    star)
        TABLE_PATH="star_salmon/salmon"
        ;;
    *)
        echo "Error: Unknown aligner '$ALIGNER'" >&2
        exit 1
        ;;
esac

echo "Using aligner: $ALIGNER with $TABLE_PATH"


runName=da2_${ALIGNER}_rpts_minAbund${minAbund}_minSamples${minSamples}_${suffix}


samplesheet=${WORK_DIR}/samplesheet_all.csv
contrasts=${WORK_DIR}/contrasts.csv
matrix=${WORK_DIR}/${TABLE_PATH}.merged.gene_counts.tsv
txptlength=${WORK_DIR}/${TABLE_PATH}.merged.gene_lengths.tsv

nextflow run nf-core/differentialabundance -profile rnaseq,singularity --input ${samplesheet} --outdir ${WORK_DIR}/${runName} \
        -r 2.0.0 -c $CONFIG_FILE \
	--gtf $gtfFile   \
	--matrix ${matrix} \
	--transcript_length_matrix ${txptlength} \
	--contrasts ${contrasts} \
	--deseq2_shrink_lfc ${lfcShrink} --filtering_min_abundance ${minAbund} --filtering_min_samples ${minSamples}
