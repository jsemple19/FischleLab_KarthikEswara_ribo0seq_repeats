## add annotation to results tables with gene names and locations.
## Combine all results tables into a single .rds object
## extract table of Number of significant up/down regulated genes by different thresholds

## some special considerations for repeat families:
## If you are going to just look at family level alignment you can use star and randomly
## assign a multimapping read to a single locus
## If you want to get per locus information, it is better to use salmon pseudo alignment
## as it assigns multimappers based on EM from other reads.

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
library(DESeq2)
library(stringr)


theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9),
          title=ggtext::element_markdown(size=9)
    )
)

options(width=10000)

serverPath="/Volumes/meister.data"
#serverPath="Z:/MeisterLab"
workDir=paste0(serverPath,"/FischleLab_KarthikEswara/ribo0seq")
runName="/results_telescope"
prefix=""

setwd(workDir)
dir.create(paste0(workDir,runName,"/rds"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(workDir,runName,"/txt"), showWarnings = FALSE, recursive = TRUE)
genomeVer<-"WS298"

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)
samplesheet<-read.csv(paste0(workDir,"/samplesheet_all.csv"),sep=",",header=T)

dfamFam<-read.delim(paste0(serverPath,"/publicData/genomes/dfam35/repeatFamilies.tsv"),header=T,sep="\t")
dfamFam$classification<-NULL
#rmsk<-import("/Volumes/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq_squire/ce11_rmsk_TE.gtf")

## collect counts
countsTable<-NULL
#i=1
for(i in 1:length(samplesheet$sample)){
  sample<-samplesheet$sample[i]
  if(sample %in% colnames(countsTable)){
    next()
  }
  tmpCnts<-read.delim(paste0(workDir,"/telescope/",sample,"/",sample,"-telescope_report.tsv"),skip=1)
  tmpCnts_sample<-tmpCnts[,c("transcript","final_count")]
  colnames(tmpCnts_sample)<-c("gene_id",sample)
  if(is.null(countsTable)){
    countsTable<-tmpCnts_sample
  } else {
    countsTable<-left_join(countsTable,tmpCnts_sample,by="gene_id")
  }
}

dir.create(paste0(workDir,runName,"/txt"),showWarnings = F)
write.table(countsTable,paste0(workDir,runName,"/txt/allCounts.tsv"))

# sum counts per rptName
rpts<-countsTable[!grepl("WBGene",countsTable$gene_id),]
rpts$gene_id<-gsub("_rpt.*$","",rpts$gene_id)
#rpts$gene_id<-gsub("-I_CE","_CE",rpts$gene_id)
#rpts$gene_id<-gsub("-LTR_CE","_CE",rpts$gene_id)
rptFam<-tibble(rpts) %>% group_by(gene_id) %>%
  summarise(across(everything(),\(x) sum(x,na.rm=T)),countRpt=n()) %>% as.data.frame()

dfamFam<-read.delim(paste0(serverPath,"/publicData/genomes/dfam35/repeatFamilies.tsv"),header=T)
dfamFam$classification<-NULL

rptFam<-left_join(rptFam,dfamFam,by=join_by("gene_id"=="name"))
row.names(rptFam)<-rptFam$gene_id

#rptFam[is.na(rptFam)]<-0


# do DESeq2 on repeat families
results<-NULL
#i=2
for(i in 1:nrow(contrasts)){
  contrast<-contrasts$id[i]
  # if(contrast %in% c("lin61_vs_N2","lin61_vs_HPL2GFP__lin61")){
  #   next()
  # }
  print(contrast)
  ctrl<-contrasts$reference[i]
  treat<-contrasts$target[i]
  ctrls<-samplesheet$sample[grepl(paste0("^",ctrl,"$"),samplesheet$strain)]
  treats<-unique(samplesheet$sample[grepl(paste0("^",treat,"$"),samplesheet$strain)])
  tmpss<-samplesheet[samplesheet$strain %in% c(ctrl,treat),]
  tmpss<-tmpss[!duplicated(tmpss$sample),]
  tmpss<-tmpss[match(colnames(rptFam)[colnames(rptFam) %in% tmpss$sample], tmpss$sample),]
  #rownames(tmpss)<-tmpss$sample
  tmpss$replicate<-factor(tmpss$replicate)
  tmpss$strain<-factor(tmpss$strain)

  dds<-DESeqDataSetFromMatrix(countData=round(rptFam[,colnames(rptFam) %in% tmpss$sample],0),
                              colData=tmpss,
                              design=~replicate+strain)
  # Repeat-family matrices are sparse; poscounts handles zeros during size-factor estimation.
  dds <- DESeq(dds, sfType="poscounts")
  resultsNames(dds)
  res <- lfcShrink(dds, contrast=c("strain",treat,ctrl),type="ashr")
  res$strain<-contrast
  res$gene_id<-rownames(res)
  res<-left_join(data.frame(res),dfamFam,by=join_by("gene_id"=="name"))
  #res$famSize<-rptFam$count.y
  #res$avgLength<-rptFam$avgLength
  if(is.null(results)){
    results<-res
  } else {
    results<-rbind(results,res)
  }
}

dir.create(paste0(workDir,runName,"/rds"),showWarnings=F,recursive=T)
saveRDS(results,paste0(workDir,runName,"/rds/allResults.rds"))
results
#head(rptFam)
