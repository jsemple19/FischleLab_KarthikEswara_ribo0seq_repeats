library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)

serverDir="/Volumes/meister.data"
workDir<-paste0(serverDir,"/FischleLab_KarthikEswara/ribo0seq/")

removeOutliers<-T
removeRepeats<-T


## remove extreme outliers
genomeVer<-"WS298"
gtf<-import(paste0(serverDir,"/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.gtf"))

gtf <- gtf[gtf$type == "gene"]
mcols(gtf)<-mcols(gtf)[c("source","type","gene_id","gene_biotype","gene_name")]
gtf$source<-genomeVer
gtf<-sort(gtf)

# find major ribosomal RNA clusters
idx<-grep("rrn",gtf$gene_name)
rrr<-gtf[idx,]
rrr[seqnames(rrr)=="I"]
rrr[seqnames(rrr)=="V"]

riboRNA_I<-GRanges("I:15060299-15071033")
riboRNA_V<-GRanges("V:17115526-17132107")

rrI<-subsetByOverlaps(gtf,riboRNA_I,ignore.strand=T)
rrI
rrV<-subsetByOverlaps(gtf,riboRNA_V,ignore.strand=T)
rrV
sp<-gtf[grep("srpr",gtf$gene_name)]
sp

toRemove<-c(rrI$gene_id,rrV$gene_id,sp$gene_id)
length(toRemove) # 42 genes
write.csv(toRemove,paste0(workDir,"/geneIDs_toRemove.csv"))



counts<- read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv"))
dim(counts)
lengths<- read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_lengths.tsv"))
dim(lengths)

if(removeOutliers){
  counts<-counts[!(counts$gene_id %in% toRemove),]
  dim(counts)

  lengths<-lengths[!(lengths$gene_id %in% toRemove),]
  dim(lengths)
}

if(removeRepeats){
  counts<-counts[!grepl("_rpt",counts$gene_id),]
  dim(counts)
  lengths<-lengths[!grepl("_rpt",lengths$gene_id),]
  dim(lengths)
}


write.table(counts,paste0(workDir,"/star_salmon/salmon.merged.gene_counts",
                          ifelse(removeOutliers,"_noRR_noSP",""),
                          ifelse(removeRepeats,"_noRpts",""),
                          ".tsv"),sep="\t",row.names=F,quote=F)
write.table(lengths,paste0(workDir,"/star_salmon/salmon.merged.gene_lengths",
                           ifelse(removeOutliers,"_noRR_noSP",""),
                           ifelse(removeRepeats,"_noRpts",""),
                           ".tsv"),sep="\t",row.names=F,quote=F)
