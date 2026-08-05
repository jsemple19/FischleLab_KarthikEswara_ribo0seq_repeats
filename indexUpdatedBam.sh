#!/bin/bash
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --array=7-27%9
#SBATCH --mem-per-cpu=32G

TMPDIR=/scratch/jsemple/tmp_tele_$SLURM_ARRAY_TASK_ID
mkdir -p $TMPDIR

WORK_DIR=/mnt/meister.data/FischleLab_KarthikEswara/ribo0seq_telescope
CONFIG_FILE=/mnt/meister.data/nf-core/unibe_izb.config

SAMPLE_SHEET=$WORK_DIR/bamlist.csv

# Print the task index.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

line_num=$((SLURM_ARRAY_TASK_ID + 1))
line=$(sed -n "${line_num}p" $SAMPLE_SHEET)

IFS=',' read -r -a sample_data <<< "$line"
sample_name=${sample_data[0]}
#bamFile=${sample_data[1]}

mkdir -p $WORK_DIR/bam
namesorted="$WORK_DIR/telescope/${sample_name}/${sample_name}-updated.bam"
possorted="$WORK_DIR/telescope/${sample_name}/${sample_name}-updated.sort.bam"

if [ ! -f "$possorted" ]; then
  echo "position sorting " $filename
  samtools sort -o $possorted -T $TMPDIR -@ $SLURM_CPUS_PER_TASK $namesorted
  samtools index -@ $SLURM_CPUS_PER_TASK $possorted
fi


rm -rf $TMPDIR
echo "tmpdir cleanup finished" && date
