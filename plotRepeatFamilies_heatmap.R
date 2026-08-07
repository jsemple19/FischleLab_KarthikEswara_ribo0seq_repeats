library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)



theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9),
          title=ggtext::element_markdown(size=9)
    )
)


serverPath="/Volumes/meister.data"
#serverPath="Z:/MeisterLab"

workDir=paste0(serverPath,"/FischleLab_KarthikEswara/ribo0seq")
runName="/results_telescope"
#runName="/results_repeats_kallisto"

dir.create(paste0(workDir,runName,"/plots"), showWarnings = FALSE, recursive = TRUE)

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)
batch="lin61" #lin61, N2 or all
prefix=paste0(batch,"_")
genomeVer="WS298"
setwd(workDir)
minFamSize=1

# raw tpm

results<-readRDS(paste0(workDir,runName,"/rds/allResults.rds"))
results$strain<-factor(results$strain,levels=contrasts$id)
results<-results[!grepl("rRNA",results$gene_id),]
results<-results[results$gene_id!="__no_feature",]
results<-results[results$rptAnnCount>minFamSize,]
length(unique(results$gene_id)) #177

# log2fc ------
results$sig<-F
results$sig[results$padj<0.05]<-T
results$avrExpr<-data.frame(results) |>  group_by(gene_id) |>  mutate(avrExpr=mean(baseMean)) %>%
  ungroup() %>% pull(avrExpr)
#results<-results[results$avrExpr>20,]

if(batch=="lin61"){
  lin61contrasts<-grep("vs_EM88$",levels(results$strain),value=T)
  res<-results[results$strain %in% lin61contrasts,]
  colOrder<-c("EM90_vs_EM88",  "EM92_vs_EM88","EM91_vs_EM88")
}else if(batch=="N2"){
  n2contrasts<-grep("vs_N2$",levels(results$strain),value=T)
  res<-results[results$strain %in% n2contrasts,]
  colOrder<-c("Btm_vs_N2", "hpl2tm_vs_N2","lin61tm_vs_N2", "EM38_vs_N2")
}else{
  res<-results
  colOrder<-c("Btm_vs_N2", "hpl2tm_vs_N2","lin61tm_vs_N2", "EM38_vs_N2", "EM90_vs_EM88",  "EM92_vs_EM88","EM91_vs_EM88")
}

sigGenes<-unique(res$gene_id[res$sig])
resSig<-res[res$gene_id %in% sigGenes,]

resWide<-pivot_wider(resSig,id_cols=c("gene_id","class","type","rptAnnCount"),names_from="strain",values_from="log2FoldChange")

resWide<-resWide[,c("gene_id","class","type","rptAnnCount",colOrder)]

resWide$class[is.na(resWide$class)]<-"IR"
resWide$class<-factor(resWide$class)

levels(resWide$class)<-c("ClassI\n(RNA)","ClassII\n(DNA)","IR","Satellite")

table(resWide$class)
table(resWide$type)

rowOrder<-order(resWide$class,resWide$type,resWide$gene_id)
resWide<-resWide[rowOrder,]
#resWide$type<-factor(resWide$type,levels=unique(resWide$type))

mat<-as.matrix(resWide[,c(colOrder)])
row.names(mat)<-gsub("_CE$","",resWide$gene_id)


# add type annotation
type_colors <- structure(
  rep(c("darkgreen", "orange"),length(unique(resWide$type))/2),
  names = unique(resWide$type)
)

row_ha <- rowAnnotation(
  Type = anno_text(
    resWide$type,
    gp = gpar(col = type_colors[resWide$type], fontsize = 8),
    location = 0.5,
    just = "center"
  )
)


pdf(paste0(workDir, runName,"/plots/",prefix,"repeats_heatmap.pdf"),
    width = 5, height = ifelse(batch=="lin61",5.5,7))

ht<-Heatmap(mat,cluster_columns = F, cluster_rows=F, name="Log2FC",
        row_split=resWide$class, row_labels=rownames(mat),
        column_split=sapply(strsplit(colnames(mat),"_"),"[[",3),
        right_annotation = rowAnnotation("log10(n)" = anno_barplot(log10(resWide$rptAnnCount)),
        annotation_name_gp=gpar(fontsize=10)),
        left_annotation=row_ha,
        column_gap = unit(3, "mm"), row_gap = unit(2,"mm"),
        row_names_gp = gpar(fontsize = 8),row_title_rot = 0,
        row_title_gp= gpar(fontsize= 10),column_title_gp=gpar(fontsize=10),
        column_names_gp=gpar(fontsize=10))
draw(ht)
dev.off()





