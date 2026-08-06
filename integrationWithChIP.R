# Plotting significant gene number and LFC by chromosome, chromosome region (arm/center),
# volcano plots and heatmaps

library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggpubr)
library(plotly)
library(ggrepel)
library(rstatix)
library(htmlwidgets)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyr)
library(grid)
library(gridExtra)
library(ggplotify)
library(RColorBrewer)
library(eulerr)
library(ComplexUpset)
library(ggVennDiagram)

options(width=100)

theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9),
          plot.title=ggtext::element_markdown(size=9)
    )
)

serverPath="/Volumes/meister.data"
#serverPath="Z:/MeisterLab"

workDir=paste0(serverPath,"/FischleLab_KarthikEswara/ribo0seq")
runName="/da1_star_canonical_minAbund5_minSamples3_lfcShrink_noRR_noSP_noRpts"

source(paste0(workDir,"/functions_finalFigures.R"))

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)
contrasts<-contrasts[contrasts$id!="EM38_vs_N2",]

prefix="noEM38_"

setwd(workDir)

dir.create(paste0(workDir,runName,"/custom/integrationWithChIP"), showWarnings = FALSE, recursive = TRUE)


# proper triple peaks
#triple<-readRDS(paste0(workDir,runName,"/custom/rds/triple_ChIPpeaks.RDS"))

# annotate with chipseq data for all samples -----
# genes
results<-readRDS(paste0(workDir,runName,"/custom/rds/.results_annotated.RDS"))
results$padj[is.na(results$padj)]<-1
results<-results[results$group %in% contrasts$id,]
results$group<-factor(results$group,levels=c("hpl2tm_vs_N2","lin61tm_vs_N2","Btm_vs_N2",
                                             "EM90_vs_EM88","EM92_vs_EM88","EM91_vs_EM88"))
#results$group<-relevel(results$group,ref="EM38_vs_N2")
results<-results[results$seqnames!="MtDNA" & results$seqnames!="chrM",]

res<-tableToGranges(results)
seqlevelsStyle(res)<-"UCSC"

#arms vs center
domains<-read.delim(paste0(serverPath,"/publicData/various/Rockman_Kruglyak_2009_armsVcenter/celegans_rockman_kruglyak_2009_table1_domains_ce11_bp.tsv"))
domaingr<-GRanges(domains)
seqlevelsStyle(domaingr)<-"UCSC"

res$chrRegionType<-NA
ol<-findOverlaps(resize(res,width=1,fix="start"),domaingr)
res$chrRegionType[queryHits(ol)]<-domaingr$DomainType[subjectHits(ol)]

table(res$chrRegionType,res$group)

#res$triplePeaks<-countOverlaps(res,triple)
#table(res$chrRegionType[res$triplePeaks>0])


# peaks overlapping individually
bedFiles<-list.files(paste0(workDir,"/../GSE271919_ChIPseq_processed/bed/"),pattern="_N2")
targets<-sapply(strsplit(bedFiles,"_"),"[[",1)
chosen<-c(3,5,7)
beddf<-data.frame(bedFiles=paste0(workDir,"/../GSE271919_ChIPseq_processed/bed/",bedFiles[chosen]),targets=targets[chosen])
beddf


for(i in 1:nrow(beddf)){
  bed<-import(beddf$bedFiles[i],format="broadPeak")
  seqlevelsStyle(bed)<-"UCSC"
  cname<-paste0(beddf$target[i],"_peaks")
  mcols(res)[cname]<-countOverlaps(res,bed)
}

peakCols<-grep("_peaks",colnames(mcols(res)))
res$peakTypeNumber<-rowSums(data.frame(mcols(res)[,peakCols])>0)

# signal
bwFiles<-list.files(paste0(workDir,"/../GSE271919_ChIPseq_processed/bigwig/"),pattern="_N2_merged")
targets<-sapply(strsplit(bwFiles,"_"),"[[",1)
chosen<-c(3,5,7)
bwdf<-data.frame(bwFiles=paste0(workDir,"/../GSE271919_ChIPseq_processed/bigwig/",bwFiles[chosen]),targets=targets[chosen])
bwdf

bwdf$totalScore<-0
i=1
for(i in 1:nrow(bwdf)){
  bw<-import(bwdf$bwFiles[i])
  bwdf$totalScore[i]<-sum(bw$score)
  seqlevelsStyle(bw)<-"UCSC"
  seqlevels(bw)<-seqlevels(Celegans)
  bw<-sort(bw)
  bw<-bw[seqnames(bw)!="chrM"]
  bw<-dropSeqlevels(bw,value="chrM")
  cov<-coverage(bw,weight="score")
  cname<-paste0(bwdf$target[i],"_signal")
  print(cname)
  res<-binnedAverage(res,cov,varname=cname,na.rm=T)
}

saveRDS(data.frame(res),paste0(workDir,runName,"/custom/rds/_ChIPpeaks.results_annotated.RDS"))



# start plotting ------

res<-readRDS(paste0(workDir,runName,"/custom/rds/_ChIPpeaks.results_annotated.RDS"))

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)
contrastsToKeep<-c(3,4,1,5,7,6)
contrasts<-contrasts[contrastsToKeep,]

levels(res$group)<-contrasts$id
## euler diagram -----

ss<-res |> data.frame() |>  filter(chrRegionType %in% c("arm","tip"),
                   seqnames %in% seqnames(Celegans)[1:5],
                   group=="EM90_vs_EM88")

table(ss$group)

peakCols<-grep("_peaks",colnames(ss))
binMat<-ss[,peakCols]>0
row.names(binMat)<-ss$gene_id
binMat
#colors <- brewer.pal(3, "Accent")
#named_colors <- setNames(colors,colnames(binMat))
named_colors<-getColorCombos(binMat)
fit<-euler(binMat)
p1<-plot(fit,quantities=T, fills = list(fill = named_colors, alpha = 0.7))
p1

ggsave(filename = paste0(workDir, runName, "/custom/integrationWithChIP/euler_autosomalArmTripleChIPpeaks.pdf"),
       plot = p1, width = 15, height = 15,units="cm")


## boxplots -------
if(length(levels(res$group))>3){
  colorSet<-c(brewer.pal(nrow(contrasts)-3, "Dark2"),"#6accdd","#f5ed20","#bf60a5")
  colorSet3<-c("#6accdd","#f5ed20","#bf60a5")
} else {
  colorSet<-c("#6accdd","#f5ed20","#bf60a5")
}
names(colorSet)<-contrasts$id

ss<-res |> data.frame() |>  filter(chrRegionType %in% c("arm","tip"),
                   seqnames %in% seqnames(Celegans)[1:5],
                   peakTypeNumber==3)

table(ss$peakTypeNumber,ss$group)

obsCounts<-data.frame(ss) %>% group_by(group) %>%
  summarize(count = n())
ss$batch<-"EM88"
ss$batch[grepl("_vs_N2",ss$group)]<-"N2"
ss1<-ss[grepl("_vs_N2",ss$group),]
ss1$group<-droplevels(ss1$group)
tt1<-data.frame(ss1) %>%
  wilcox_test(log2FoldChange~group, p.adjust.method="fdr") %>%
  p_format(p.adj, new.col=T,accuracy=1e-32)
ss2<-ss[grepl("_vs_EM88",ss$group),]
ss2$group<-droplevels(ss2$group)
tt2<-data.frame(ss2) %>%
  wilcox_test(log2FoldChange~group,p.adjust.method="fdr") %>%
  p_format(p.adj, new.col=T,accuracy=1e-32)

ylimits<-c(-1.5,2)
p2<-ggplot(ss,aes(x=group,y=log2FoldChange)) +
  geom_boxplot(aes(fill=group),outlier.shape=NA,notch=T)  +
  scale_fill_manual(values = colorSet)+
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey40")+
  #geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt1,label="p.adj.format",y.position=ylimits[2]*0.6,step.increase=0.023,
                     angle=0, hide.ns=T,tip.length=0.01,label.size=3) +
  stat_pvalue_manual(data=tt2,label="p.adj.format",y.position=ylimits[2]*0.6,step.increase=0.023,
                     angle=0, hide.ns=T,tip.length=0.01,label.size=3) +
  labs(title=paste0(length(unique(ss$gene_id)),
                    " genes on autosomal arms<br />with HPL-2/LIN-61/H3K9me2 ChIP peaks"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank()) +
  geom_vline(xintercept=3.5)
p2

ggsave(filename = paste0(workDir, runName, "/custom/integrationWithChIP/boxplot_lfcByGroup_autosomalArmTripleChIPpeaks.pdf"),
       plot = p2, width = 12, height = 12,units="cm")




## violinplots -----
ylimits<-c(-3,6)
p2aa<-ggplot(ss1,aes(x=group,y=log2FoldChange,fill=group)) +
  geom_violin(width=0.9)+
  geom_boxplot(width=0.12,color="white",alpha=1,outlier.shape=NA)  +
  scale_fill_manual(values = colorSet)+
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey40")+
  #geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt1,label="p.adj.format",y.position=ylimits[2]*0.7,step.increase=0.08,
                     angle=0, hide.ns=T,tip.length=0.01,label.size=3) +
  labs(title=paste0(length(unique(ss$gene_id)),
                    " genes on autosomal arms<br>with HPL-2/LIN-61/H3K9me2 ChIP peaks"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p2aa

## violinplot
ylimits<-c(-3,4)
p2ab<-ggplot(ss2,aes(x=group,y=log2FoldChange,fill=group)) +
  geom_violin(width=0.9)+
  geom_boxplot(width=0.12,color="white",alpha=1,outlier.shape=NA)  +
  scale_fill_manual(values = c("#6accdd","#f5ed20","#bf60a5"))+
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey40")+
  #geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt2,label="p.adj.format",y.position=ylimits[2]*0.8,step.increase=0.08,
                     angle=0, hide.ns=T,tip.length=0.01,label.size=3) +
  labs(title=paste0(length(unique(ss$gene_id)),
                    " genes on autosomal arms<br>with HPL-2/LIN-61/H3K9me2 ChIP peaks"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p2ab

p2a<-ggarrange(p2aa,p2ab,ncol=2)

ggsave(filename = paste0(workDir, runName, "/custom/integrationWithChIP/violin_lfcByGroup_autosomalArmTripleChIPpeaks_white.pdf"),
       plot = p2a, width = 12, height = 12,units="cm")



## boxplots by quantile ------
ss<-ss[ss$batch=="EM88",]

for(i in 1:nrow(bwdf)){
  print(bwdf$target[i])
  print(bwdf$totalScore[i])
  chipValues<-as.matrix(ss[paste0(bwdf$target[i],"_signal")])
  ss[,paste0(bwdf$target[i],"_quantile")]<-factor(cut(chipValues,breaks=quantile(chipValues,probs=c(0,0.33,0.66,1)),labels=F,include.lowest=T))
}

ylimits<-c(-1.5,2.5)
obsCounts<-data.frame(ss) %>% group_by(group,hpl2_quantile) %>%
  summarize(count = n())
tt<-data.frame(ss) %>% group_by(hpl2_quantile) %>%
  wilcox_test(log2FoldChange~group, p.adjust.method="fdr") %>%
  p_format(p.adj, new.col=T,accuracy=1e-32) %>%
  mutate(y.position=ylimits[2]*0.75)
#rep(seq(ylimits[2]*0.5,ylimits[2]*0.95,ylimits[2]*0.45/5),3))


p3<-ggplot(ss,aes(x=group,y=log2FoldChange)) +
  geom_boxplot(aes(fill=group),outlier.shape=NA)  +
  scale_fill_manual(values = colorSet)+
  coord_cartesian(ylim=ylimits) +
  facet_wrap(.~hpl2_quantile)+
  geom_hline(yintercept=0,linetype="dashed",color="grey40")+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt,label="p.adj.format",
                     angle=0, hide.ns=T,tip.length=0.01,label.size=3) +
  labs(title=paste0(length(unique(ss$gene_id)),
                    " genes on autosomal arms with HPL-2/LIN-61/H3K9me2 ChIP peaks<br>grouped by HPL-2 ChIP signal quantile"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "right",
        axis.title.x=element_blank(),axis.text.x=element_blank())
p3


ylimits<-c(-1.5,2.5)
obsCounts<-data.frame(ss) %>% group_by(group,lin61_quantile) %>%
  summarize(count = n())
tt<-data.frame(ss) %>% group_by(lin61_quantile) %>%
  wilcox_test(log2FoldChange~group, p.adjust.method="fdr") %>%
  p_format(p.adj, new.col=T,accuracy=1e-32) %>%
  mutate(y.position=ylimits[2]*0.75)
#rep(seq(ylimits[2]*0.5,ylimits[2]*0.95,ylimits[2]*0.45/5),3))

p4<-ggplot(ss,aes(x=group,y=log2FoldChange)) +
  geom_boxplot(aes(fill=group),outlier.shape=NA) +
  scale_fill_manual(values = colorSet)+
  coord_cartesian(ylim=ylimits) +
  facet_wrap(.~lin61_quantile)+
  geom_hline(yintercept=0,linetype="dashed",color="grey40")+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt,label="p.adj.format",
                     angle=0, hide.ns=T,tip.length=0.01,label.size=3) +
  labs(title=paste0(length(unique(ss$gene_id)),
                    " genes on autosomal arms with HPL-2/LIN-61/H3K9me2 ChIP peaks<br>grouped by LIN-61 ChIP signal quantile"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "right",
        axis.title.x=element_blank(),axis.text.x=element_blank())
p4


ylimits<-c(-1.5,2.5)
obsCounts<-data.frame(ss) %>% group_by(group,H3K9me2_quantile) %>%
  summarize(count = n())
tt<-data.frame(ss) %>% group_by(H3K9me2_quantile) %>%
  wilcox_test(log2FoldChange~group, p.adjust.method="fdr") %>%
  p_format(p.adj, new.col=T,accuracy=1e-32) %>%
  mutate(y.position=ylimits[2]*0.75)
#rep(seq(ylimits[2]*0.5,ylimits[2]*0.95,ylimits[2]*0.45/5),3))

p5<-ggplot(ss,aes(x=group,y=log2FoldChange)) +
  geom_boxplot(aes(fill=group),outlier.shape=NA)  +
  scale_fill_manual(values = colorSet)+
  coord_cartesian(ylim=ylimits) +
  facet_wrap(.~H3K9me2_quantile)+
  geom_hline(yintercept=0,linetype="dashed",color="grey40")+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt,label="p.adj.format",
                     angle=0, hide.ns=T,tip.length=0.01,label.size=3) +
  labs(title=paste0(length(unique(ss$gene_id)),
                    " genes on autosomal arms with HPL-2/LIN-61/H3K9me2 ChIP peaks<br>grouped by H3K9me2 ChIP signal quantile"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p5

p<-ggarrange(p3,p4,p5,nrow=3,heights=c(1.3,1.3,2))
p
ggsave(filename = paste0(workDir, runName, "/custom/integrationWithChIP/lfcByGroup_autosomalArmTripleChIPpeaks_quantile_wilcoxTest.pdf"),
       plot = p, width = 15, height = 30,units="cm")


## output file lists of upregulated genes-------

res<-readRDS(paste0(workDir,runName,"/custom/rds/_ChIPpeaks.results_annotated.RDS"))

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)
contrastsToKeep<-c(3,4,1,5,7,6)
contrasts<-contrasts[contrastsToKeep,]

padjVal=0.05
lfcVal=0.5

dir.create(paste0(workDir,runName,"/custom/GO"),showWarnings=F)
dir.create(paste0(workDir,runName,"/custom/bed"),showWarnings=F)

### all significant lfc>0.5 autosomal arm up genes
ss<-res |> data.frame() |>  filter(res$group %in% contrasts$id,
                                    chrRegionType %in% c("arm","tip"),
                                   seqnames %in% seqnames(Celegans)[1:5],
                                   padj<padjVal,log2FoldChange>lfcVal)

ss$group<-droplevels(ss$group)

table(ss$group)

forBed<-GRanges(ss)
forBed$name<-forBed$gene_id
forBed$score<- forBed$log2FoldChange #-log10(forBed$padj)

for(g in levels(ss$group)){
  write.table(ss$gene_id[ss$group==g],paste0(workDir,runName,"/custom/GO/",g,"_upGenes_padj0.05_lfc0.5_autosomal_arms.txt"),row.names=F, quote=F,col.names=F)
  print(length(forBed[forBed$group==g]))
  forBed_sorted <- forBed[order(mcols(forBed)$score, decreasing = TRUE)]
  export.bed(forBed_sorted[forBed_sorted$group==g],paste0(workDir,runName,"/custom/bed/",g,"__sigUp_autosomalArms_lfcSorted.bed"))
}

### all significant lfc>0.5 autosomal arm up genes with triple peaks
ss<-res |> data.frame() |>  filter(res$group %in% contrasts$id,
                                   chrRegionType %in% c("arm","tip"),
                                   seqnames %in% seqnames(Celegans)[1:5],
                                   padj<padjVal,log2FoldChange>lfcVal,
                                   peakTypeNumber==3)

ss$group<-droplevels(ss$group)

table(ss$group)

# forBed<-GRanges(ss)
# forBed$name<-forBed$gene_id
# forBed$score<- -log10(forBed$padj)

dir.create(paste0(workDir,runName,"/custom/bed"),showWarnings=F)

for(g in levels(ss$group)){
  write.table(ss$gene_id[ss$group==g],paste0(workDir,runName,"/custom/GO/",g,"_upGenes_padj0.05_lfc0.5_autosomal_arms_triplePeaks.txt"),row.names=F, quote=F,col.names=F)
  # print(length(forBed[forBed$group==g]))
  # forBed_sorted <- forBed[order(mcols(forBed)$score, decreasing = TRUE)]
  # export.bed(forBed_sorted[forBed_sorted$group==g],paste0(workDir,runName,"/custom/bed/",g,"__sigUp_autosomalArms_triplePeaks_mlogpadjSorted.bed"))
}

