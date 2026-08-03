#! /usr/bin/bash
#SBATCH --time=0-08:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --ntasks=2

source $CONDA_ACTIVATE deeptools

serverPath=/mnt/meister.data
workDir=${serverPath}/FischleLab_KarthikEswara/ribo0seq
runName=da1_star_canonical_minAbund5_minSamples3_lfcShrink
regionData="autosomalArms_lfcSorted"
chipData="publicChIPseq"
sortRegions="descend" #descend, ascend or no
redoMatrix=false
outDir=${workDir}/${runName}/custom/deeptools_upGenes
mkdir -p $outDir

echo "sort option? " ${sortRegions}

bigwigPath=${serverPath}/FischleLab_KarthikEswara/GSE271919_ChIPseq_processed/bigwig

sampleNames=( "HPL-2" "LIN-61" "H3K9me2" )

fileNames=( ${bigwigPath}/hpl2_N2_merged.bwa_aln.rmdup.bamCompare.bw \
   ${bigwigPath}/lin61_N2_merged.bwa_aln.rmdup.bamCompare.bw \
   ${bigwigPath}/H3K9me2_N2_merged.bwa_aln.rmdup.bamCompare.bw )


echo ${sampleNames[@]}
echo ${fileNames[@]}


## bigwig files
bwPaths=${fileNames[@]}
bwNames=${sampleNames[@]}

## bed files
upRegN2=( `ls ${workDir}/${runName}/custom/bed/*_vs_N2__sigUp_${regionData}.bed` )
upRegNamesN2=( `basename -a ${upRegN2[@]} | awk -F"__" '{print $1}'` )

upRegEM88=( `ls ${workDir}/${runName}/custom/bed/*_vs_EM88__sigUp_${regionData}.bed` )
upRegNamesEM88=( `basename -a ${upRegEM88[@]} | awk -F"__" '{print $1}'` )


echo "number of tasks: " $SLURM_NTASKS

## compute matrix
if $redoMatrix; then
      blacklist=https://github.com/Boyle-Lab/Blacklist/raw/master/lists/ce11-blacklist.v2.bed.gz
      blacklistFile=`basename ${blacklist%.gz}`
      if [ ! -f "$blacklistFile" ]; then
            wget $blacklist
            gunzip `basename $blacklist`
      fi

      echo "computing matrix"
      computeMatrix scale-regions -S ${bwPaths[@]} \
            -R ${upRegN2[@]} \
            --beforeRegionStartLength 1000 \
            --regionBodyLength 2000 \
            --afterRegionStartLength 1000 \
            --skipZeros \
            -o ${outDir}/matrix_${regionData}_upN2_${chipData}.mat.gz \
            --outFileNameMatrix ${outDir}/matrix_${regionData}_upN2_${chipData}.tab \
            --blackListFileName $blacklistFile \
            --numberOfProcessors $SLURM_NTASKS \
            --verbose

      echo "computing matrix"
      computeMatrix scale-regions -S ${bwPaths[@]} \
            -R ${upRegEM88[@]} \
            --beforeRegionStartLength 1000 \
            --regionBodyLength 2000 \
            --afterRegionStartLength 1000 \
            --skipZeros \
            -o ${outDir}/matrix_${regionData}_upEM88_${chipData}.mat.gz \
            --outFileNameMatrix ${outDir}/matrix_${regionData}_upEM88_${chipData}.tab \
            --blackListFileName $blacklistFile \
            --numberOfProcessors $SLURM_NTASKS \
            --verbose
fi


maxs=( 20 20 20 )
mins=( `printf ' 0 %.0s' {1..3}` )

## make plots
echo "plotting heatmap"
plotHeatmap -m ${outDir}/matrix_${regionData}_upN2_${chipData}.mat.gz \
      -out ${outDir}/${regionData}_upN2_${chipData}_heatmap_${sortRegions}sort.pdf \
      --colorMap Blues  \
      --startLabel "TSS" --endLabel "TES" \
      -y "" -x "Distance" \
      --regionsLabel ${upRegNamesN2[@]}  \
      --plotTitle "Autosomal arm up regulated genes" \
      --samplesLabel ${bwNames[@]} \
      --zMin ${mins[@]} \
      --zMax ${maxs[@]} \
      --yMin ${mins[@]} \
      --yMax ${maxs[@]} \
      --sortRegions ${sortRegions} \
      --sortUsingSamples 1

plotHeatmap -m ${outDir}/matrix_${regionData}_upEM88_${chipData}.mat.gz \
      -out ${outDir}/${regionData}_upEM88_${chipData}_heatmap_${sortRegions}sort.pdf \
      --colorMap Blues  \
      --startLabel "TSS" --endLabel "TES" \
      -y "" -x "Distance" \
      --regionsLabel ${upRegNamesEM88[@]}  \
      --plotTitle "Autosomal arm up regulated genes" \
      --samplesLabel ${bwNames[@]} \
      --zMin ${mins[@]} \
      --zMax ${maxs[@]} \
      --yMin ${mins[@]} \
      --yMax ${maxs[@]} \
      --sortRegions ${sortRegions} \
      --sortUsingSamples 1

