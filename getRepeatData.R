library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)

genome<-Celegans
dfamVer="Dfam_3.5"

repeatDir<-"/Volumes/meister.data/publicData/genomes/dfam35"

dir.create(repeatDir)

print(paste0("workDir is: ", workDir))
setwd(workDir)


dfamURL=paste0("https://www.dfam.org/releases/",dfamVer,"/annotations/ce10/ce10.nrph.hits.gz")
download.file(url=dfamURL, destfile=paste0(repeatDir,"/ce10_",dfamVer,".nrph.hits.gz"),method="wget")
system(paste0("gunzip ",repeatDir,"/ce10_",dfamVer,".nrph.hits.gz"))

repeats_dfam <- read.delim(paste0(repeatDir,"/ce10_",dfamVer,".nrph.hits"))

repeats.gr <- GRanges(seqnames=repeats_dfam[,"X.seq_name"],
                      ranges=IRanges(start=apply(repeats_dfam[,c("env.st",
                                                                 "env.en")],1,FUN=min),
                                     end=apply(repeats_dfam[,c("env.st",
                                                               "env.en")],1,FUN=max)),
                      strand=repeats_dfam[,"strand"])

repeats.gr$source<-dfamVer
repeats.gr$type<-"repeat"

repeats.gr$Name<-repeats_dfam$family_name
repeats.gr$Alias<-repeats_dfam$family_acc


liftoverURL="http://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz"
download.file(url=liftoverURL,destfile=paste0(repeatDir,"/ce10ToCe11.over.chain.gz"),method="wget")
system(paste0("gunzip ",repeatDir,"/ce10ToCe11.over.chain.gz"))

chain <- import(paste0(repeatDir,"/ce10ToCe11.over.chain"))

repeats_ce11.gr <- unlist(range(liftOver(repeats.gr,chain)))
table(width(repeats_ce11.gr)-width(repeats.gr)) #checking the differences are minor
mcols(repeats_ce11.gr)<-mcols(repeats.gr)
uniqueID<-data.frame(repeats_ce11.gr) %>% group_by(Name) %>%
  mutate(uniqueID=paste0(Name,"_rpt",row_number())) %>% select(uniqueID)
repeats_ce11.gr$unique_id<-uniqueID$uniqueID
saveRDS(repeats_ce11.gr,paste0(repeatDir,"/repeats_ce11_",dfamVer,"_nr.rds"))



# bed
repeats_bed<-repeats_ce11.gr
mcols(repeats_bed)<-NULL
repeats_bed$score<-1
repeats_bed$name<-repeats_ce11.gr$Name
export.bed(repeats_bed,paste0(repeatDir,"/repeats_ce11_",dfamVer,"_nr.bed"))


# gtf
repeats_gtf<-repeats_ce11.gr
repeats_gtf$type<-"exon"
repeats_gtf$gene_id<-repeats_ce11.gr$Name
repeats_gtf$transcript_id<-repeats_ce11.gr$unique_id
repeats_gtf$gene_name<-repeats_ce11.gr$Name
repeats_gtf$gene_biotype<-"repeat"
repeats_gtf$Name<-NULL
repeats_gtf$Alias<-NULL
repeats_gtf$unique_id<-NULL

repeats_gtf

# gene version
repeats_as_genes<-repeats_gtf
repeats_as_genes$type<-"gene"
repeats_as_genes$transcript_id<-""

# transcript version
repeats_as_txpt<-repeats_gtf
repeats_as_txpt$type<-"transcript"


repeats_gtf<-c(repeats_as_genes,repeats_as_txpt,repeats_gtf) %>% sort()


export(repeats_gtf,paste0(repeatDir,"/repeats_ce11_",dfamVer,"_nr.gtf"))

seqinfo(repeats_gtf)
seqlevels(repeats_gtf)<-gsub("chr","",seqlevels(repeats_gtf))
seqinfo(repeats_gtf)
export(repeats_gtf,paste0(repeatDir,"/repeats_WB_",dfamVer,"_nr.gtf"))

# then manually combine this with the canonical genes file from worm base (removing headerlines)
# geneGTF="/mnt/meister.data/publicData/genomes/WS298/c_elegans.PRJNA13758.WS298.canonical_geneset.gtf"
# repeatGTF=/mnt/meister.data/publicData/genomes/dfam35/repeats_WB_Dfam_3.5_nr.gtf
# combinedGTF=/mnt/meister.data/publicData/genomes/dfam35/WS298_canonicalgenes_Dfam3.5nr_repeats.gtf
# grep '^#'  <(cat "$geneGTF" "$repeatGTF")  > combined_headers.txt
# grep -v '^#'  <(cat "$geneGTF" "$repeatGTF") > combined_lines.txt
# cat combined_headers.txt combined_lines.txt > $combinedGTF
# rm combined_headers.txt combined_lines.txt

# rename gene_id with unique ids so it doesn't give error for attribute
gtf<-import("/Volumes/meister.data/publicData/genomes/dfam35/WS298_canonicalgenes_Dfam3.5nr_repeats.gtf")

gtf_genes<-gtf[gtf$source!="Dfam_3.5"]

rpt_exon<-gtf[gtf$source=="Dfam_3.5" & gtf$type=="exon"]
rpt_txpt<-gtf[gtf$source=="Dfam_3.5" & gtf$type=="transcript"]
rpt_gene<-gtf[gtf$source=="Dfam_3.5" & gtf$type=="transcript"]
rpt_gene$type<-"gene"

rptgtf<-sort(c(rpt_exon,rpt_txpt,rpt_gene))
rptgtf$gene_id<-rptgtf$transcript_id

gtf_all<-c(gtf_genes,rptgtf)

export(gtf_all,paste0(repeatDir,"/WS298_canonicalgenes_Dfam3.5nr_repeats.uniq_gene_id.gtf"))

gtf_all<-import(paste0(repeatDir,"/WS298_canonicalgenes_Dfam3.5nr_repeats.uniq_gene_id.gtf"))

gtf_rpts<-gtf_all[gtf$source=="Dfam_3.5"]

export(gtf_rpts,paste0(repeatDir,"/Dfam3.5nr_repeats.uniq_gene_id.gtf"))
