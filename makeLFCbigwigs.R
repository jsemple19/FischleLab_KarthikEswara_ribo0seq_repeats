library(rtracklayer)
library(GenomicRanges)
library(DESeq2)
library(BSgenome.Celegans.UCSC.ce11)

serverPath="/Volumes/meister.data"
#serverPath="Z:/MeisterLab"

workDir=paste0(serverPath,"/FischleLab_KarthikEswara/ribo0seq")
runName="/da1_star_canonical_minAbund5_minSamples3_lfcShrink_noRR_noSP_noRpts"

source(paste0(workDir,"/functions_finalFigures.R"))

batch="EM88"
#batch="EM88"

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)
contrasts<-contrasts[contrasts$id!="EM38_vs_N2",]

prefix=""

setwd(workDir)

dir.create(paste0(workDir,runName,"/custom/lfcBigwigs"), showWarnings = FALSE, recursive = TRUE)

## data -------
results<-readRDS(paste0(workDir,runName,"/custom/rds/.results_annotated.RDS"))
results$padj[is.na(results$padj)]<-1
results<-results[results$group %in% contrasts$id,]

gr<-tableToGranges(results,sort=FALSE)
seqlevelsStyle(gr)<-"UCSC"
seqinfo(gr)<-seqinfo(Celegans)
strand(gr)<-"*"

gr$score<-gr$log2FoldChange

g=results$group[1]
for(g in levels(results$group)){
  subgr<-gr[gr$group==g]
  # Create disjoint intervals from all your ranges
  disjoint_gr <- disjoin(subgr)

  # For each disjoint interval, find which original genes overlap it
  hits <- findOverlaps(disjoint_gr, subgr)

  # Aggregate score per disjoint interval — choose your logic:
  agg_score <- tapply(subgr$score[subjectHits(hits)], queryHits(hits),
                      FUN = function(x) x[which.max(abs(x))])   # or max, or which.max by baseMean, etc.

  disjoint_gr$score <- NA
  disjoint_gr$score[as.integer(names(agg_score))] <- agg_score
  disjoint_gr <- disjoint_gr[!is.na(disjoint_gr$score)]

  export.bw(sort(disjoint_gr), paste0(workDir,runName,"/custom/lfcBigwigs/",g,"_log2FoldChange.bw"))
}

