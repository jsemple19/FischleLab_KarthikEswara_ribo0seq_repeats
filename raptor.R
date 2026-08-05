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
runName="/_2_raptorAgePrior_star_canonical_minAbund5_minSamples3_noRR_noSP_noRpts"
prefix=""
setwd(workDir)

dir.create(paste0(workDir,runName,"/custom/raptor"), showWarnings = FALSE, recursive = TRUE)

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)
contrasts<-contrasts[contrasts$id!="EM38_vs_N2",]

samplesheet<-read.csv(paste0(workDir,"/samplesheet_all.csv"), header=T, stringsAsFactors = F)

samplesheet<-samplesheet[samplesheet$strain %in% c(contrasts$target, contrasts$reference),]

## raptor ------
minAbund=5
minSamples=3
runNum=2
lfcShrink=T
pc=F
prior=c(220,200) # either mean and range of prior or NULL

runName=paste0("/_",runNum,"_raptorAge",
               ifelse(is.null(prior),"","Prior"),
               "_star_canonical_minAbund",minAbund,
               "_minSamples",minSamples,
               ifelse(pc,"_pc",""),
               "_noRR_noSP_noRpts")

prefix=""
genomeVer<-"WS298"

source(paste0(workDir,"/raptor_functions.R"))
col.palette1 <- c('grey20', 'firebrick', 'royalblue', 'forestgreen')
#mapMethod="salmon"
#prefix=paste0(mapMethod,"_raptor_age_deseq2_pc_nomito_lfcShrink_")

counts<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts_noRR_noSP_noRpts.tsv"))
tpm<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_tpm.tsv"))
tpm<-tpm[tpm$gene_id %in% counts$gene_id,]


gtf<-import(paste0(serverPath,"/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.gtf"))
gtf <- gtf[gtf$type == "gene"]
mcols(gtf)<-mcols(gtf)[c("source","type","gene_id","gene_biotype","gene_name")]
gtf$source<-genomeVer
gtf<-sort(gtf)

if(pc){
  gtf<-gtf[gtf$gene_biotype=="protein_coding"]
  counts<-counts[counts$gene_id %in% gtf$gene_id,]
  tpm<-tpm[tpm$gene_id %in% gtf$gene_id,]
}


idx<-which(colnames(counts) %in% samplesheet$sample)
counts<-counts[,c(1,2,idx)]
tpm<-tpm[,c(1,2,idx)]


dir.create(paste0(workDir,runName,"/plots"), recursive=T, showWarnings = FALSE)
# use same directory structure as for nf-core pipelines
dir.create(paste0(workDir,runName,"/tables/differential/rnaseq"), recursive=T, showWarnings = FALSE)



# load reference
list_refs(datapkg="wormRef")
ref <- prepare_refdata("Cel_embryo", "wormRef", 600)


# create list object with rnaseq data
rnaseqData=list()

# get counts and tpm
counts<-counts[rowSums(counts[, c(-1,-2)]>minAbund)>minSamples,]
rnaseqData$count<-round(counts[, c(-1,-2)],0)
row.names(rnaseqData$count)<-counts$gene_id
colnames(rnaseqData$count)<-gsub("^X","",colnames(rnaseqData$count))

tpm<-tpm[tpm$gene_id %in% row.names(rnaseqData$count),]
rnaseqData$tpm<- tpm[, c(-1,-2)]
row.names(rnaseqData$tpm)<-tpm$gene_id
colnames(rnaseqData$tpm)<-gsub("^X","",colnames(rnaseqData$tpm))

# get phenotype data  (metadata) for samples
pdatOrder<-match(colnames(rnaseqData$tpm),samplesheet$sample)
rnaseqData$pdat<-samplesheet[pdatOrder,]

str(rnaseqData$pdat)
rnaseqData$pdat$strain<-factor(rnaseqData$pdat$strain)
#rnaseqData$pdat$strain<-factor(rnaseqData$pdat$strain)
rnaseqData$pdat$replicate<-factor(rnaseqData$pdat$replicate)

summary(rnaseqData$pdat)

# estimate ages
# estimate ages
if(is.null(prior)){
  ageEstimate <- ae(rnaseqData$tpm, ref)
} else {
  ageEstimate <- ae(rnaseqData$tpm, ref,prior=prior[1],prior.params=prior[2])
}


print(ageEstimate)
par(mfrow=c(1,1))

# plot ages
pdf(paste0(workDir,runName,"/plots/",prefix,"ages.pdf"),width=9,height=7)
plot(ageEstimate,group=rnaseqData$pdat$strain, main="Ages of samples by strain",
     xlab="Age (min)", cex=0.5, col=col.palette1[rnaseqData$pdat$replicate])
#plot(ageEstimate,group=rnaseqData$pdat$group,main="Ages of samples by group",
#     xlab="Age (h)", cex=0.5, col=col.palette1[rnaseqData$pdat$replicate])
plot(ageEstimate,group=rnaseqData$pdat$replicate, main="Ages of samples by replicate",
     xlab="Age (min)", cex=0.5, col=col.palette1[rnaseqData$pdat$group])
dev.off()

# plot full correlation
pdf(paste0(workDir,runName,"/plots/",prefix,"_full_age_correlations.pdf"),width=9,height=6)
par(mfrow=c(2,3))
plot_cor(ageEstimate,subset=1:6)
plot_cor(ageEstimate,subset=7:length(unique(samplesheet$sample)))
#plot_cor(ageEstimate,subset=13:17)
dev.off()
par(mfrow=c(1,1))


## look at contribution development has to expression
rnaseqData_rc <- ref_compare(
  X = log1p(rnaseqData$tpm), # sample data, log(TPM+1)
  ref = ref, # ref object
  ae_obj = ageEstimate, # ae object
  group = rnaseqData$pdat$strain # factor defining compared groups (wt/mut)
)

print(rnaseqData)
plotList=list()
c=2
for(c in 1:nrow(contrasts)){
  rnaseqData_lfc<-get_logFC(rnaseqData_rc,
                            l=contrasts$target[c],
                            l0=contrasts$reference[c],
                            verbose=T)

  plotList[[contrasts$id[c]]]<-gg_logFC(rnaseqData_lfc,
                                        main = contrasts$id[c],
                                        xlab = contrasts$id[c],
                                        ylab = "Developmental log2FC")
}

p<-ggarrange(plotlist=plotList,nrow=2, ncol=3)
ggsave(paste0(workDir,runName,"/plots/",prefix,"_development_lfc.pdf"),p,height=20, width=29, units="cm")


## add Z normalised ages to metadata as age covariate
idx<-match(rnaseqData$pdat$sample,row.names(ageEstimate$age.estimates))
rnaseqData$pdat$age<-ageEstimate$age.estimates[idx,1]
rnaseqData$pdat$ageZnorm<- (rnaseqData$pdat$age-mean(rnaseqData$pdat$age))/sd(rnaseqData$pdat$age)

# make sure at least one sample has > 5 reads
#rnaseqData$count <- rnaseqData$count[apply(rnaseqData$count, 1, max)>5, ]

# use age instead of batch
dd0 <- DESeqDataSetFromMatrix(countData = rnaseqData$count,
                              colData = rnaseqData$pdat,
                              design = ~ageZnorm+replicate+strain)
dd0<-DESeq(dd0)


for(c in 1:nrow(contrasts)){
  res <- results(dd0,contrast=c("strain",contrasts$target[c],contrasts$reference[c]),alpha=0.05)
  if(lfcShrink){
    res <- lfcShrink(dd0,contrast=c("strain",contrasts$target[c],contrasts$reference[c]),alpha=0.05,type="ashr")
  }
  summary(res)
  write.table(res, file=paste0(workDir,runName,"/tables/differential/rnaseq/",prefix,contrasts$id[c],"_study.deseq2.results.tsv"), sep="\t", quote=F)
}



vsd <- vst(dd0, blind=TRUE)
minipdat<-rnaseqData$pdat[,c("sample","ageZnorm","replicate","strain")]
rownames(minipdat)<-minipdat$sample
vsdpca<-pca(assay(vsd),metadata=minipdat)
p1<-biplot(vsdpca,colby="strain",gridlines.major=F,gridlines.minor=F)
p1a<-screeplot(vsdpca,gridlines.major=F,gridlines.minor=F)
p2<-biplot(vsdpca,x="PC3",y="PC4",colby="strain",gridlines.major=F,gridlines.minor=F)
p3<-biplot(vsdpca,x="PC5",y="PC6",colby="strain",gridlines.major=F,gridlines.minor=F)

p<-ggarrange(p1,p1a,p2,p3,nrow=2,ncol=2)
p
ggsave(paste0(workDir,runName,"/plots/",prefix,"_PCA.pdf"),p,height=25, width=25, units="cm")


file.copy(from=paste0(workDir,"/raptor_deseq2.R"),
          to=paste0(workDir,runName,"/raptor_deseq2.R"),
          overwrite=TRUE)

print(runName)


