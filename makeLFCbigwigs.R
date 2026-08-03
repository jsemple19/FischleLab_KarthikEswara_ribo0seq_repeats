library(rtracklayer)
library(GenomicRanges)
library(DESeq2)
library(BSgenome.Celegans.UCSC.ce11)

serverPath="/Volumes/meister.data"
#serverPath="Z:/MeisterLab"

workDir=paste0(serverPath,"/FischleLab_KarthikEswara/ribo0seq")
runName="/da1_star_canonical_minAbund5_minSamples3_lfcShrink"

source(paste0(workDir,"/functions_finalFigures.R"))

batch="N2"
#batch="EM88"

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)

prefix=""

setwd(workDir)

dir.create(paste0(workDir,runName,"/custom/lfcBigwigs"), showWarnings = FALSE, recursive = TRUE)

## data -------
results<-readRDS(paste0(workDir,runName,"/custom/rds/.results_annotated.RDS"))
results$padj[is.na(results$padj)]<-1

gr<-tableToGranges(results,sort=FALSE)
seqlevelsStyle(gr)<-"UCSC"
seqinfo(gr)<-seqinfo(Celegans)
strand(gr)<-"*"

gr$score<-gr$log2FoldChange

for(g in levels(results$group)){
  # Create disjoint intervals from all your ranges
  disjoint_gr <- disjoin(gr[gr$group==g])

  # For each disjoint interval, find which original genes overlap it
  hits <- findOverlaps(disjoint_gr, gr)

  # Aggregate score per disjoint interval — choose your logic:
  agg_score <- tapply(gr$score[subjectHits(hits)], queryHits(hits),
                      FUN = function(x) x[which.max(abs(x))])   # or max, or which.max by baseMean, etc.

  disjoint_gr$score <- NA
  disjoint_gr$score[as.integer(names(agg_score))] <- agg_score
  disjoint_gr <- disjoint_gr[!is.na(disjoint_gr$score)]

  export.bw(sort(disjoint_gr), paste0(workDir,runName,"/custom/lfcBigwigs/",g,"_log2FoldChange.bw"))
}
