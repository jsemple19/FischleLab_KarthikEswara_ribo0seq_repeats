#!/bin/bash
#SBATCH --time=1-12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --array=9
#SBATCH --mem-per-cpu=96G

TMPDIR=/scratch/${USER}/tmp_tele_$SLURM_ARRAY_TASK_ID
mkdir -p $TMPDIR
rerunNSORT=false


genomeVer=WS298
genomeDir=/mnt/meister.data/publicData/genomes
genomeFile=${genomeDir}/${genomeVer}/c_elegans.PRJNA13758.${genomeVer}.genomic.fa

gtfFile=${genomeDir}/dfam35/${genomeVer}_canonicalgenes_Dfam3.5nr_repeats.uniq_gene_id.gtf

WORK_DIR=/mnt/meister.data/FischleLab_KarthikEswara/ribo0seq
CONFIG_FILE=/mnt/meister.data/nf-core/unibe_izb.config


#singularity pull docker://quay.io/biocontainers/telescope:1.0.3_fix--py36h87e0c26_0
TELESCOPE_SIF=/mnt/meister.data/containers/telescope_1.0.3_fix--py36h87e0c26_0.sif

SAMPLE_SHEET=$WORK_DIR/bamlist.csv

# Print the task index.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

line_num=$((SLURM_ARRAY_TASK_ID + 1))
line=$(sed -n "${line_num}p" $SAMPLE_SHEET)

IFS=',' read -r -a sample_data <<< "$line"
sample_name=${sample_data[0]}
bamFile=${sample_data[1]}

mkdir -p $WORK_DIR/bam
filename=`basename $bamFile`
namesorted="${filename/.sorted.bam/.nsorted.bam}"
if [  "$rerunNSORT" ]; then
  rm $WORK_DIR/bam/$namesorted
  rm $TMPDIR/$filename*
  echo "name sorting " $filename
  samtools collate -o $WORK_DIR/bam/$namesorted -T $TMPDIR/$filename --output-fmt BAM -@  $SLURM_CPUS_PER_TASK $bamFile
fi

rm -rf $TMPDIR/*

mkdir -p $WORK_DIR/telescope/$sample_name
echo "running telescope for " $namesorted
# need --cleanenv to use container R rather than system R
singularity exec --cleanenv --bind ${genomeDir} --bind $TMPDIR $TELESCOPE_SIF telescope assign $WORK_DIR/bam/$namesorted  $gtfFile --outdir $WORK_DIR/telescope/$sample_name --exp_tag $sample_name --tempdir $TMPDIR --attribute gene_id --updated_sam

echo "telescope run finished" && date

rm -rf $TMPDIR
echo "tmpdir cleanup finished" && date
