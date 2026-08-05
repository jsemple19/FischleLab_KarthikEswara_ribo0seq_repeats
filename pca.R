rm(list = ls())

library(RAPToR)
library(wormRef)
library(rtracklayer)
library(DESeq2)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggpubr)
#library(plotly)
#library(ggrepel)
#library(rstatix)
library(PCAtools)

options(tible.width=10000)

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

workDir=paste0(serverPath,"/FischleLab_KarthikEswara/ribo0seq")
runName="/da1_star_canonical_minAbund5_minSamples3_lfcShrink_noRR_noSP_noRpts"
prefix=""
setwd(workDir)

dir.create(paste0(workDir,runName,"/custom/pca"), showWarnings = FALSE, recursive = TRUE)

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)
contrasts<-contrasts[contrasts$id!="EM38_vs_N2",]

samplesheet<-read.csv(paste0(workDir,"/samplesheet_all.csv"), header=T, stringsAsFactors = F)

samplesheet<-samplesheet[samplesheet$strain %in% c(contrasts$target, contrasts$reference),]

## pca ------
# from precomputed vsd
prefix="all"
vsd<-read.delim(paste0(workDir,runName,"/tables/processed_abundance/rnaseq/all.vst.tsv"))

strains<-c("N2","hpl2tm","lin61tm","Btm","EM88","EM90","EM92","EM91")
colourSet<-c("grey60",brewer.pal(4, "Dark2"),"#6accdd","#f5ed20","#bf60a5")
names(colourSet)<-strains

minipdat<-samplesheet[,c("sample","strain")]
minipdat<-minipdat[!duplicated(minipdat),]
rownames(minipdat)<-minipdat$sample
minipdat$strain<-factor(minipdat$strain,levels=strains)

mat<-vsd[,colnames(vsd) %in% samplesheet$sample]
row.names(mat)<-vsd$gene_id

removeVar=NULL
vsdpca<-pca(mat,metadata=minipdat,center=T,removeVar=removeVar)
p1<-biplot(vsdpca,colby="strain",lab=NULL,
           gridlines.major=F,gridlines.minor=F,
           legendPosition='right',
           legendIconSize = 2,
           pointSize=2)

# remove the existing point layer(s)
p1$layers <- p1$layers[!sapply(p1$layers, function(l) inherits(l$geom, "GeomPoint"))]
# add your own styled points
p1 <- p1 +
  geom_point(aes(fill = col), shape = 21, color = "black",
             size = 3, alpha = 0.7, stroke = 0.8) +
  scale_fill_manual(values=colourSet)

p1

ggsave(paste0(workDir,runName,"/custom/pca/",prefix,"_PCA",
              ifelse(is.null(removeVar),"_all",paste0("_top",(1-removeVar)*100,"pc")),".pdf"),
       p1,height=15, width=20, units="cm")




## pca ------
# from precomputed vsd
prefix="EM88"
vsd<-read.delim(paste0(workDir,runName,"/tables/processed_abundance/rnaseq/all.vst.tsv"))

strains<-c("EM88","EM90","EM92","EM91")
colourSet<-c("grey60","#6accdd","#f5ed20","#bf60a5")
names(colourSet)<-strains

samplesheet

minipdat<-samplesheet[samplesheet$strain %in% strains,c("sample","strain")]
minipdat<-minipdat[!duplicated(minipdat),]
rownames(minipdat)<-minipdat$sample
minipdat$strain<-factor(minipdat$strain,levels=strains)

mat<-vsd[,colnames(vsd) %in% minipdat$sample]
row.names(mat)<-vsd$gene_id

vsdpca<-pca(mat,metadata=minipdat)
p1<-biplot(vsdpca,colby="strain",lab=NULL,
           gridlines.major=F,gridlines.minor=F,
           legendPosition='right',
           legendIconSize = 2,
           pointSize=2)

# remove the existing point layer(s)
p1$layers <- p1$layers[!sapply(p1$layers, function(l) inherits(l$geom, "GeomPoint"))]
# add your own styled points
p1 <- p1 +
  geom_point(aes(fill = col), shape = 21, color = "black",
             size = 3, alpha = 0.7, stroke = 0.8) +
  scale_fill_manual(values=colourSet)

p1

ggsave(paste0(workDir,runName,"/custom/pca/",prefix,"_PCA.pdf"),p1,height=15, width=20, units="cm")

# from dds object (rnaseq pipeline):
# load(paste0(workDir,"/star_salmon/deseq2_qc/deseq2.dds.RData"))
# dds<-DESeq(dds)
#
# vsd <- vst(dds, blind=TRUE)
# orderedcols<-match(rownames(minipdat),colnames(vsd))
# vsd <- vsd[, orderedcols]
#
# minipdat<-samplesheet[,c("sample","replicate","strain")]
# minipdat<-minipdat[!duplicated(minipdat),]
# rownames(minipdat)<-minipdat$sample
# vsdpca<-pca(assay(vsd),metadata=minipdat,center=T,scale=F)
# p1<-biplot(vsdpca,colby="strain",gridlines.major=F,gridlines.minor=F)
# p1
# result qualitatively similar


