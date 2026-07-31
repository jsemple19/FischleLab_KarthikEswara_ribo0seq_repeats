# Plotting significant gene number and LFC by chromosome, chromosome region (arm/center),
# volcano plots and heatmaps

library(GenomicRanges)
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
runName="/da1_star_canonical_minAbund5_minSamples3_lfcShrink"

source(paste0(workDir,"/functions_finalFigures.R"))

batch="N2"
#batch="EM88"

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)


if(batch == "N2"){
  contrastsToKeep<-c(3,4,1)
  contrasts<-contrasts[contrastsToKeep,]
} else if(batch =="EM88"){
  contrastsToKeep<-c(5,7,6)
  contrasts<-contrasts[contrastsToKeep,]
} else {
  print("unrecognised batch name")
}

prefix=paste0(batch,"_")


setwd(workDir)

dir.create(paste0(workDir,runName,"/custom/upregulatedOnArms"), showWarnings = FALSE, recursive = TRUE)

#arms vs center Rockman (2009)
domains<-read.delim(paste0(serverPath,"/publicData/various/Rockman_Kruglyak_2009_armsVcenter/celegans_rockman_kruglyak_2009_table1_domains_ce11_bp.tsv"))
rockman<-GRanges(domains)
seqlevelsStyle(rockman)<-"UCSC"

## data -------
results<-readRDS(paste0(workDir,runName,"/custom/rds/.results_annotated.RDS"))
results$padj[is.na(results$padj)]<-1

gr<-tableToGranges(results,sort=FALSE)
seqlevelsStyle(gr)<-"UCSC"
ol<-findOverlaps(resize(gr,width=1,fix="start"),domaingr,ignore.strand=T)
gr$DomainType[queryHits(ol)] <- paste0(domaingr$DomainType[subjectHits(ol)])
gr$chrRegion[queryHits(ol)] <- paste0(seqnames(domaingr)[subjectHits(ol)],"_",domaingr$DomainType[subjectHits(ol)])
res<-as.data.frame(gr)
res<-res[res$seqnames != "chrM",]
saveRDS(res,paste0(workDir,runName,"/custom/rds/",prefix,"_regionType.results_annotated.RDS"))
contrasts

lfcVal=0.5
padjVal=0.05

### significant genes on arms -----
ss<-res %>% filter(group %in% contrasts$id,
               DomainType %in% c("arm", "tip"),
               seqnames %in% seqnames(Celegans)[1:5],
               padj<padjVal,
               log2FoldChange>lfcVal)

ss$group<-factor(ss$group,levels=contrasts$id)
table(ss$seqnames,ss$group)
table(ss$group)


# keep values in all samples for genes significant in at least one
genesToKeep<-unique(ss$gene_id)
ss<-res[res$gene_id %in% genesToKeep & res$group %in% contrasts$id,]
ss$group<-droplevels(ss$group)
ss$group<-factor(ss$group,levels=contrasts$id)
table(ss$seqnames,ss$group)
saveRDS(ss, paste0(workDir,runName,"/custom/rds/",prefix,"_subsetByGene.results_annotated.RDS"))


rm(results)
rm(ss)
rm(genesToKeep)

## Heatmap of significantly upregulated autosomal genes -----
chromosomes<-"autosomal"
direction<-"up"
numSamplesSignificant=1

res<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,"_subsetByGene.results_annotated.RDS"))

mat_padj<-gatherResults(res,valueColumn="padj")
mat_padj[is.na(mat_padj)]<-1
mat_lfc<-gatherResults(res,valueColumn="log2FoldChange")
sig<-(mat_padj[,2:ncol(mat_padj)]<padjVal & mat_lfc[,2:ncol(mat_lfc)]>lfcVal)
sigGenes<-mat_padj$gene_id[rowSums(sig)>=numSamplesSignificant]
mat_lfc<-mat_lfc[mat_lfc$gene_id %in% sigGenes,]
rownames(mat_lfc)<-mat_lfc$gene_id
mat_lfc<-mat_lfc[,2:ncol(mat_lfc)]
mat_lfc<-as.matrix(mat_lfc)
colnames(mat_lfc)<-gsub("log2FoldChange_","",colnames(mat_lfc))

clusters<-data.frame(name=apply(sig,1,function(x) paste0(gsub("padj_","",names(x)[x]),collapse=" & ")),
                     rank=rowSums(sig))

groups<-clusters[!duplicated(clusters),]
groups<-groups[order(-groups$rank),]
#groups<-groups[c(1:4,7,5,6),]
groups
clusters$name<-factor(clusters$name,levels=groups$name)

mat_list <- lapply(split(seq_len(nrow(mat_lfc)), clusters$name), function(i) mat_lfc[i, , drop = FALSE])

plotTitle<-paste0(nrow(mat_lfc)," ", chromosomes," arm genes sig. in \n>=",
                  numSamplesSignificant," sample (padj<",padjVal," ",
                  ifelse(direction=="both","|LFC|","LFC"),
                  ifelse(direction=="down","< -",">"),lfcVal,")")
# Shared color mapping
col_fun <- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("blue", "white", "red"))

pdf(paste0(workDir, runName,"/custom/upregulatedOnArms/",prefix,"hclust_heatmap_sigSamples",
           numSamplesSignificant,"_padj",padjVal,"_lfc",lfcVal,"_",direction,"_",chromosomes,"Chr.pdf"),
    width=4,height=9)

ht<-Heatmap(mat_lfc,show_row_names=F,
            cluster_columns=F, cluster_rows=T,  show_row_dend = F,
            column_names_rot=90, column_title=plotTitle,
            row_title_rot=0,
            column_names_gp=gpar(fontsize=10),
            col=col_fun,
            column_title_gp = gpar(fontsize = 12, fontface = "bold", just = 0),
            heatmap_legend_param=list(title=expression("log"[2]~"FC")),
            )
ht<-draw(ht)
dev.off()


pdf(paste0(workDir, runName,"/custom/upregulatedOnArms/",prefix,"hclust_heatmap_sigSamples",
           numSamplesSignificant,"_padj",padjVal,"_lfc",lfcVal,"_",direction,"_",chromosomes,"Chr_grouped.pdf"),
    width=8,height=11)

# Create a heatmap for each matrix
ht_list <- lapply(names(mat_list), function(name) {
  Heatmap(mat_list[[name]],row_title=name, show_row_names=F,
          cluster_columns=F, cluster_rows=T,  show_row_dend = F,
          column_names_rot=45, column_title=plotTitle,
          column_names_gp = gpar(fontsize = 10),
          column_title_gp = gpar(fontsize = 10, fontface = "bold", just = -5),
          row_title_rot=0, row_title_gp=gpar(fontsize=10),
          col=col_fun,
          heatmap_legend_param=list(title=expression("log"[2]~"FC"))
  )
})

# Combine them into one plot
ht_combined <- Reduce(`%v%`, ht_list)
draw(ht_combined,merge_legend=T, heatmap_legend_side="right")
dev.off()

sapply(mat_list,nrow)



## significant genes boxplots  & violin plots -----
res<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,"_subsetByGene.results_annotated.RDS"))

if(batch=="N2"){
  colorSet<-brewer.pal(nrow(contrasts), "Dark2")
} else if(batch=="EM88"){
  colorSet<-c("#6accdd","#f5ed20","#bf60a5")
}
inside<-rep("white",nrow(contrasts))
names(colorSet)<-contrasts$id
names(inside)<-contrasts$id

obsCounts<-data.frame(res) %>% group_by(group) %>%
  summarize(count = n())
tt<-data.frame(res) %>%
  wilcox_test(log2FoldChange~group, p.adjust.method="fdr") %>%
  p_format(p.adj, new.col=T,accuracy=1e-300)
ylimits<-c(-1.5,2)
p1<-ggplot(res,aes(x=group,y=log2FoldChange)) +
  geom_boxplot(aes(fill=group),outlier.shape=NA) +
  scale_fill_manual(values = colorSet)+
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey40") +
  stat_pvalue_manual(data=tt,label="p.adj.format",y.position=ylimits[2]*0.65,step.increase=0.005,
                     angle=0, hide.ns=T,tip.length=0.001,label.size=3) +
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  labs(title=paste0(length(unique(res$gene_id)),
                    " genes on autosomal<br>arms significant in >=",
                    numSamplesSignificant," sample"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p1
ggsave(filename = paste0(workDir, runName, "/custom/upregulatedOnArms/",prefix,"boxplot_lfcByGroup_autosomalArmUp",
                         "_padj",padjVal,"_lfc",lfcVal,".pdf"),
       plot = p1, width = 7, height = 11,units="cm")
#y.position=seq(ylimits[2]*0.95,ylimits[2]-0.83,-ylimits[2]*0.24/nrow(wilcoxt))

if(batch=="N2"){
  ylimits<-c(-3,6)
} else if(batch=="EM88"){
  ylimits<-c(-2,4)
}

p1a<-ggplot(res,aes(x=group,y=log2FoldChange,fill=group)) +
  geom_violin(width=0.9)+
  geom_boxplot(width=0.12,color="white",alpha=1,outlier.shape=NA)  +
  scale_fill_manual(values = colorSet)+
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey40") +
  stat_pvalue_manual(data=tt,label="p.adj.format",y.position=ylimits[2]*0.75,step.increase=0.01,
                     angle=0, hide.ns=T,tip.length=0.001,label.size=3) +
  labs(title=paste0(length(unique(res$gene_id)),
                    " genes on autosomal<br>arms significant in >=",
                    numSamplesSignificant," sample"),
       y="log<sub>2</sub>FC")+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p1a
ggsave(filename = paste0(workDir, runName, "/custom/upregulatedOnArms/",prefix,"violin_lfcByGroup_autosomalArmUp",
                         "_padj",padjVal,"_lfc",lfcVal,".pdf"),
       plot = p1a, width = 7, height = 11,units="cm")


## all arm genes - boxplots and violin plots -----
results<-readRDS(paste0(workDir,runName,"/custom/rds/.results_annotated.RDS"))

gr<-tableToGranges(results,sort=FALSE)
seqlevelsStyle(gr)<-"UCSC"
ol<-findOverlaps(resize(gr,width=1,fix="start"),rockman,ignore.strand=T)
gr$DomainType[queryHits(ol)] <- paste0(rockman$DomainType[subjectHits(ol)])
gr$chrRegion[queryHits(ol)] <- paste0(seqnames(rockman)[subjectHits(ol)],"_",rockman$DomainType[subjectHits(ol)])

res<-as.data.frame(gr)
ss<-res %>% filter(group %in% contrasts$id,
                   DomainType %in% c("arm", "tip"),
                   seqnames %in% seqnames(Celegans)[1:5])
ss$group<-droplevels(ss$group)
ss$group<-factor(ss$group,levels=contrasts$id)
saveRDS(ss, paste0(workDir,runName,"/custom/rds/",prefix,"_allArmGenes.results_annotated.RDS"))

obsCounts<-data.frame(ss) %>% group_by(group) %>%
  summarize(count = n())
tt<-data.frame(ss) %>%
  wilcox_test(log2FoldChange~group, p.adjust.method="fdr") %>%
  p_format(p.adj,new.col=T,accuracy=1e-100)
ylimits<-c(-1.5,2)
p2<-ggplot(ss,aes(x=group,y=log2FoldChange)) +
  geom_boxplot(aes(fill=group),outlier.shape=NA) +
  scale_fill_manual(values = colorSet)+
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey40") +
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt,label="p.adj.format",y.position=ylimits[2]*0.8,step.increase=0.003,
                     angle=0, hide.ns=T,tip.length=0.001,label.size=3) +
  labs(title=paste0(length(unique(ss$gene_id)),
                    " (all expressed) genes<br>on autosomal arms"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p2
ggsave(filename = paste0(workDir, runName, "/custom/upregulatedOnArms/",prefix,"boxplot_lfcByGroup_autosomalArm",
                         "_allGenes.pdf"),
       plot = p2, width = 7, height = 11,units="cm")


if(batch=="N2"){
  ylimits<-c(-3,6)
} else if(batch=="EM88"){
  ylimits<-c(-2,4)
}
p2a<-ggplot(ss,aes(x=group,y=log2FoldChange,fill=group)) +
  geom_violin(width=0.9)+
  geom_boxplot(width=0.12,color="white",alpha=1,outlier.shape=NA)  +
  scale_fill_manual(values = colorSet)+
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey40") +
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt,label="p.adj.format",y.position=ylimits[2]*0.8,step.increase=0.003,
                     angle=0, hide.ns=T,tip.length=0.001,label.size=3) +
  labs(title=paste0(length(unique(ss$gene_id)),
                    " (all expressed) genes<br>on autosomal arms"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p2a
ggsave(filename = paste0(workDir, runName, "/custom/upregulatedOnArms/",prefix,"violin_lfcByGroup_autosomalArm",
                         "_allGenes.pdf"),
       plot = p2a, width = 7, height = 11,units="cm")


