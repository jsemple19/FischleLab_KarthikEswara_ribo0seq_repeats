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
library(cowplot)


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

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)
prefix=""
genomeVer="WS298"
setwd(workDir)
minFamSize=1

# raw tpm

results<-readRDS(paste0(workDir,runName,"/rds/allResults.rds"))
results$strain<-factor(results$strain,levels=contrasts$id)
results<-results[!grepl("rRNA",results$gene_id),]
results<-results[results$gene_id!="__no_feature",]
results<-results[results$rptAnnCount>minFamSize,]

# log2fc ------
results$sig<-F
results$sig[results$padj<0.05]<-T
results$avrExpr<-data.frame(results) |>  group_by(gene_id) |>  mutate(avrExpr=mean(baseMean)) %>%
  ungroup() %>% pull(avrExpr)
#results<-results[results$avrExpr>20,]


n2contrasts<-grep("vs_N2$",levels(results$strain),value=T)
lin61contrasts<-grep("vs_EM88$",levels(results$strain),value=T)

res_n2<-results[results$strain %in% n2contrasts,]
res_lin61<-results[results$strain %in% lin61contrasts,]

# n2 subset
pa<-ggplot(res_n2,aes(x=reorder(gene_id,-log2FoldChange),y=log2FoldChange,fill=sig)) +
  facet_grid(strain~.,switch = "y")+
  geom_col() +
  geom_errorbar(
    aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE),
    width = 0.2
  ) +
  geom_hline(yintercept=0,color="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none") +
  scale_fill_manual(values=c("grey40","blue"))
pa
pb<-ggplot(res_n2,aes(x=reorder(gene_id,-log2FoldChange),y=rptAnnCount)) +
  geom_col(fill="darkgreen")+
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
  )
pb
pc<-ggplot(res_n2,aes(x=reorder(gene_id,-log2FoldChange),y=avrExpr)) +
  geom_col(fill="purple")+
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
  )
pc
pd<-ggplot(res_n2,aes(x=reorder(gene_id,-log2FoldChange),y=avgLength)) +
  geom_col(fill="orange")+
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
  )
pd
p<-plot_grid(pd,pb,pc,pa,ncol=1,align="v",axis = "l",rel_heights=c(1,1,1,11))
ggsave(paste0(workDir, runName,"/plots/",prefix,"results_n2_lfc_rptFamilies",
              ifelse(minFamSize>0,paste0("_minSize",minFamSize),""),".png"),p,width=50,height=40,units="cm")


# lin-61 subset
pa<-ggplot(res_lin61,aes(x=reorder(gene_id,-log2FoldChange),y=log2FoldChange,fill=sig)) +
  facet_grid(strain~.,switch = "y")+
  geom_col() +
  geom_errorbar(
    aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE),
    width = 0.2
  ) +
  geom_hline(yintercept=0,color="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none") +
  scale_fill_manual(values=c("grey40","blue"))
pa
pb<-ggplot(res_lin61,aes(x=reorder(gene_id,-log2FoldChange),y=rptAnnCount)) +
  geom_col(fill="darkgreen")+
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
  )
pb
pc<-ggplot(res_lin61,aes(x=reorder(gene_id,-log2FoldChange),y=avrExpr)) +
  geom_col(fill="purple")+
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
  )
pc
pd<-ggplot(res_lin61,aes(x=reorder(gene_id,-log2FoldChange),y=avgLength)) +
  geom_col(fill="orange")+
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
  )
pd
p<-plot_grid(pd,pb,pc,pa,ncol=1,align="v",axis = "l",rel_heights=c(1,1,1,11))
p
ggsave(paste0(workDir, runName,"/plots/",prefix,"results_lin61_lfc_rptFamilies",
              ifelse(minFamSize>0,paste0("_minSize",minFamSize),""),".png"),p,width=50,height=35,units="cm")

# lin61 subset reordered by hpl2
order_vec<-data.frame(res_lin61) %>% filter(strain=="Btm_vs_EM88") %>%
  arrange(-log2FoldChange) %>% pull(gene_id)
df<-res_lin61[res_lin61$strain!="lin61tm_vs_EM88",]
df$gene_id<-factor(df$gene_id,levels=order_vec)
df<-data.frame(df) %>% arrange(gene_id)

pa<-ggplot(df,aes(x=gene_id,y=log2FoldChange,fill=sig)) +
  facet_grid(strain~.,switch = "y")+
  geom_col() +
  geom_errorbar(
    aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE),
    width = 0.2
  ) +
  geom_hline(yintercept=0,color="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none") +
  scale_fill_manual(values=c("grey40","blue"))
pa
pb<-ggplot(df,aes(x=gene_id,y=rptAnnCount)) +
  geom_col(fill="darkgreen")+
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
  )
pb
pc<-ggplot(df,aes(x=gene_id,y=avrExpr)) +
  geom_col(fill="purple")+
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
  )
pc
pd<-ggplot(df,aes(x=reorder(gene_id,-log2FoldChange),y=avgLength)) +
  geom_col(fill="orange")+
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
  )
pd
p<-plot_grid(pd,pb,pc,pa,ncol=1,align="v",axis = "l",rel_heights=c(1,1,1,8))
ggsave(paste0(workDir, runName,"/plots/",prefix,"results_lin61_lfc_rptFamilies_hpl2mutOrder",
              ifelse(minFamSize>0,paste0("_minSize",minFamSize),""),".png"),p,width=40,height=35,units="cm")

