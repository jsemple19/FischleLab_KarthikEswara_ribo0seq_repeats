library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(BSgenome.Celegans.UCSC.ce11)

theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9),
          title=ggtext::element_markdown(size=9)
    )
)

### find chip seq peaks that overlap between HPL-2 LIN-61 and H3K9me2
### In the end these peaks were not used, but rather genes that overlap all
### three types of peaks (but any given peak does not necessarily contain all three,
### though in reality there is a lot of overlap)


serverPath="/Volumes/meister.data"
#serverPath="Z:/MeisterLab"

workDir=paste0(serverPath,"/FischleLab_KarthikEswara/ribo0seq")
runName="/da1_star_canonical_minAbund5_minSamples3_lfcShrink"

# peaks
bedFiles<-list.files(paste0(workDir,"/../GSE271919_ChIPseq/bed/"),pattern="_N2")
targets<-sapply(strsplit(bedFiles,"_"),"[[",1)
chosen<-c(3,5,7)
beddf<-data.frame(bedFiles=paste0(workDir,"/../GSE271919_ChIPseq/bed/",bedFiles[chosen]),targets=targets[chosen])
beddf

# signal
bwFiles<-list.files(paste0(workDir,"/../GSE271919_ChIPseq/bigwig/"),pattern="_N2_merged")
targets<-sapply(strsplit(bwFiles,"_"),"[[",1)
chosen<-c(3,5,7)
bwdf<-data.frame(bwFiles=paste0(workDir,"/../GSE271919_ChIPseq/bigwig/",bwFiles[chosen]),targets=targets[chosen])
bwdf


i=1
triple<-import(beddf$bedFiles[i],format="broadPeak")
seqlevelsStyle(triple)<-"UCSC"


for(i in 2:nrow(beddf)){
  bed<-import(beddf$bedFiles[i],format="broadPeak")
  seqlevelsStyle(bed)<-"UCSC"
  cname<-paste0(beddf$target[i],"_peaks")
  mcols(triple)[cname]<-countOverlaps(triple,bed)
}

peakCols<-grep("_peaks",colnames(mcols(triple)))
triple$peakTypeNumber<-rowSums(data.frame(mcols(triple)[,peakCols])>0)
#triple<-triple[triple$peakTypeNumber==2]


bwdf$totalScore<-0
for(i in 1:nrow(bwdf)){
  bw<-import(bwdf$bwFiles[i])
  bwdf$totalScore[i]<-sum(bw$score)
  seqlevelsStyle(bw)<-"UCSC"
  seqlevels(bw)<-seqlevels(Celegans)
  bw<-sort(bw)
  bw<-bw[seqnames(bw) %in% seqlevels(triple)]
  seqlevels(bw)<-seqlevels(triple)
  cov<-coverage(bw,weight="score")
  cname<-paste0(bwdf$target[i],"_signal")
  print(cname)
  triple<-binnedAverage(triple,cov,varname=cname,na.rm=T)
}


triple<-triple[triple$peakTypeNumber==2]
#2116
saveRDS(triple,paste0(workDir,runName,"/custom/rds/triple_ChIPpeaks.RDS"))

# df<-triple |> data.frame()  |> pivot_longer(cols=c("H3K9me2_signal","hpl2_signal","lin61_signal"),names_to="target",values_to="signal")
#
# ggplot(df,aes(x=signal)) + geom_histogram() +
#   facet_grid(peakTypeNumber~target)
