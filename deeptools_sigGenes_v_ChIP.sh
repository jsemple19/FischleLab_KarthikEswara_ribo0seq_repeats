#! /usr/bin/bash
#SBATCH --time=0-08:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --ntasks=2

source $CONDA_ACTIVATE deeptools

workDir=/mnt/meister.data/FischleLab_KarthikEswara/ribo0seq
runName=da1_star_canonical_minAbund5_minSamples3_lfcShrink
regionData="autosomalArms"
chipData="publicChIPseq"
outDir=${workDir}/${runName}/custom/deeptools_upGenes
mkdir -p $outDir

bigwigPath=/mnt/meister.data/FischleLab_KarthikEswara/GSE271919_ChIPseq_processed/bigwig

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
upReg=( `ls ${workDir}/${runName}/custom/bed/*_sigUp_autosomalArms_triplePeaks.bed` )
upRegNames=( `basename -a ${upReg[@]} | awk -F"__" '{print $1}'` )
#downReg=( `ls ${workDir}/${runName}/custom/exploreChIP/bed/*autosomalArmsDown*.bed` )
#downRegNames=( `basename -a ${downReg[@]} | awk -F"__" '{print $1}'` )


echo "number of tasks: " $SLURM_NTASKS

## compute matrix
redoMatrix=true
if $redoMatrix; then
      blacklist=https://github.com/Boyle-Lab/Blacklist/raw/master/lists/ce11-blacklist.v2.bed.gz
      blacklistFile=`basename ${blacklist%.gz}`
      if [ ! -f "$blacklistFile" ]; then
            wget $blacklist
            gunzip `basename $blacklist`
      fi

      echo "computing matrix"
      computeMatrix scale-regions -S ${bwPaths[@]} \
            -R ${upReg[@]} \
            --beforeRegionStartLength 1000 \
            --regionBodyLength 2000 \
            --afterRegionStartLength 1000 \
            --skipZeros \
            -o ${outDir}/matrix_${regionData}_up_${chipData}.mat.gz \
	      --outFileNameMatrix ${outDir}/matrix_${regionData}_up_${chipData}.tab \
            --blackListFileName $blacklistFile \
            --numberOfProcessors $SLURM_NTASKS \
            --verbose

#       computeMatrix scale-regions -S ${bwPaths[@]} \
#             -R ${downReg[@]} \
#             --beforeRegionStartLength 1000 \
#             --regionBodyLength 2000 \
#             --afterRegionStartLength 1000 \
#             --skipZeros \
#             -o ${outDir}/matrix_${regionData}_down_${chipData}.mat.gz \
# 	      --outFileNameMatrix ${outDir}/matrix_${regionData}_down_${chipData}.tab \
#             --blackListFileName $blacklistFile \
#             --numberOfProcessors $SLURM_NTASKS \
#             --verbose
fi

# multiBigwigSummary BED-file --bwfiles ${bwPaths[@]} --BED $activeEnhancers_fountain \
#     --outRawCounts matrix_enhVhistoneModEncode_avrRaw.tab -o matrix_enhVhistoneModEncode_avrRaw.npz \
#     --numberOfProcessors 4

# maxs=( `awk 'NR > 1 { for (i = 4; i <= NF; i++) { if ($i > max[i]) max[i] = $i+0 } } END { for (i = 4; i <= NF; i++) print max[i]/2+0 }' matrix_enhVhistoneModEncode_avrRaw.tab` )
#maxs=( 30 `printf ' 2 %.0s' {1..21}` )
#mins=( `printf ' -2 %.0s' {1..22}` )
maxs=( 50 150 700 )
mins=( `printf ' 0 %.0s' {1..3}` )


## make plots
echo "plotting heatmap"
plotHeatmap -m ${outDir}/matrix_${regionData}_up_${chipData}.mat.gz \
      -out ${outDir}/${regionData}_up_${chipData}_heatmap.pdf \
      --colorMap Blues  \
      --startLabel "TSS" --endLabel "TES" \
      -y "" -x "Distance" \
      --regionsLabel ${upRegNames[@]}  \
      --plotTitle "Autosomal arm up regulated genes" \
      --samplesLabel ${bwNames[@]} \
      --sortRegions no \
      --zMin ${mins[@]} \
      --zMax ${maxs[@]} \
      --yMin ${mins[@]} \
      --yMax ${maxs[@]}
# --sortUsingSamples 1 2 3 \

# ## make plots
# plotHeatmap -m ${outDir}/matrix_${regionData}_down_${chipData}.mat.gz \
#       -out ${outDir}/${regionData}_down_${chipData}_heatmap.pdf \
#       --colorMap Blues  \
#       --startLabel "TSS" --endLabel "TES" \
#       -y "" -x "Distance" \
#       --regionsLabel ${downRegNames[@]} \
#       --plotTitle "Autosomal arm down regulated genes" \
#       --samplesLabel ${bwNames[@]} \
#       --sortRegions no \
#       --zMin ${mins[@]} \
#       --zMax ${maxs[@]} \
#       --yMin ${mins[@]} \
#       --yMax ${maxs[@]}
#--sortUsingSamples 1 2 3 \

# plotProfile -m ${outDir}/matrix_${regionData}_${chipData}.mat.gz  \
#               -out ${outDir}/${regionData}_${chipData}_profile.png \
#               --numPlotsPerRow 5 \
#               --regionsLabel "upReg" "downReg" \
#               --startLabel "prom" --endLabel "" \
#               --colors "cyan" "magenta" "brown" "grey" \
#               --yMax ${maxs[@]} \
#               --yMin ${mins[@]} \
#               --samplesLabel ${bwNames[@]} \
# 	      --outFileNameData ${outDir}/${regionData}_${chipData}_profile.tab \
#               --plotTitle "Promoters of COH-1 regulated genes"
