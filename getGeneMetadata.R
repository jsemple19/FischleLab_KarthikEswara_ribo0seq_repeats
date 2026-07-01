library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(dplyr)
library(tidyr)

gtf<-rtracklayer::import("/Volumes/meister.data/publicData/genomes/WS298/c_elegans.PRJNA13758.WS298.canonical_geneset.gtf")
gtf<-gtf[gtf$type=="gene"]

table(gtf$gene_biotype)
geneids<-read.csv("/Volumes/meister.data/publicData/genomes/WS298/c_elegans.PRJNA13758.WS298.geneIDs",header=F)

genes<-left_join(data.frame(gtf),geneids,by=join_by("gene_id"=="V2"))

toKeep<-c("seqnames","start","end","width","strand","source","gene_id","gene_biotype","gene_name","V4")
genes<-genes[,toKeep]
colnames(genes)[ncol(genes)]<-"sequence_id"
genegr<-GRanges(genes)

seqnames(genegr)
seqlevels(genegr)<-c("I","II","III","IV","V","X","MtDNA")
genegr<-sort(genegr)
saveRDS(genegr,"/Volumes/meister.data/publicData/genomes/WS298/c_elegans.PRJNA13758.WS298.metadata.rds")
