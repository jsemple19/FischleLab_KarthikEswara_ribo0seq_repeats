#!/bin/bash
#SBATCH --time=0-05:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

WORK_DIR=/mnt/meister.data/jsemple/ribo0seq
CONFIG_FILE=/mnt/meister.data/nf-core/unibe_izb_apptainer.config


nextflow run nf-core/fetchngs -r dev  -profile apptainer --input ids.txt --outdir $WORK_DIR -c $CONFIG_FILE --nf_core_pipeline rnaseq --download-method ftp
