library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(stringr)
library(dplyr)

workDir<-"/Volumes/meister.data/FischleLab_KarthikEswara/ribo0seq"
setwd(workDir)
dir.create(workDir,showWarnings=F, recursive=T)

# download metadata from ENA PRJEB102446
#df<-read.delim(paste0(workDir,"/filereport_read_run_PRJEB102446.tsv"))
ss<-read.csv(paste0(workDir,"/samplesheet/samplesheet.csv"))
#ids<-read.csv(paste0(workDir,"/samplesheet/id_mappings.csv"))
ss$fastq_1<-gsub("jsemple","FischleLab_KarthikEswara",ss$fastq_1)
ss$fastq_2<-gsub("jsemple","FischleLab_KarthikEswara",ss$fastq_2)


metadata<-ss[,c("sample_alias","sample_description")]
metadata$sample_name<- str_replace(metadata$sample_alias,"_Rep.?","") |>
  str_replace("_L00.?","")
metadata$replicate<- metadata$sample_alias |> str_replace("_L00.?","") |>
  str_extract("_Rep[:digit:]") |> str_replace("_Rep","") |> as.numeric()

strains<-data.frame(sample_name=unique(metadata$sample_name),
                    strain=NA)
strains$strain=c("EM88","EM90","EM92","EM91")

metadata<-left_join(metadata,strains,by="sample_name")

samplesheet<-left_join(ss[,c("sample","fastq_1","fastq_2","strandedness","sample_alias")],
          metadata[,c("sample_alias","sample_name","replicate","strain")])


samplesheet$sample_acc<-samplesheet$sample
samplesheet$sample<-paste0(samplesheet$strain,"_Rep",samplesheet$replicate)
samplesheet$sample_alias<-NULL

write.csv(samplesheet,paste0(workDir,"/samplesheet.csv"),row.names=F,
           quote=F)

# Data from Sumaiya Hasnain

ff<-read.delim(paste0(workDir,"/fileList_SH.txt"),header=F)
ff$sample<-sapply(strsplit(basename(ff$V1),"_"),"[[",1)
ff$replicate<-sapply(strsplit(basename(ff$V1),"_"),"[[",2) |> str_replace("Rep","")
ff$R1vR2<-gsub("\\.fq\\.gz","",sapply(strsplit(basename(ff$V1),"_"),"[[",3))

metadata<-data.frame(sample=c("EM88","EM90","EM92","EM91","Btm", "N2","EM38","hpl2tm","lin61tm"),
                     genotype=c("HPL2GFP_lin61","HPL2GFPneutHng_lin61","HPL2GFP3xhngI158A_lin61",
                                "HPL2GFPI158A_lin61","hpl2_lin61", "N2","HPL2GFP", "hpl2","lin61"),
                     phenotype=c("HPL2GFP_lin61","onlyDimer_lin61","onlyLLPS_lin61",
                                 "noDimerLLPS_lin61","hpl2_lin61", "N2","HPL2GFP", "hpl2","lin61"))
metadata
ff$genotype<-metadata$phenotype[match(ff$sample,metadata$sample)]

df<-data.frame(sample=paste0(ff$sample[ff$R1vR2=="1"],"_Rep",ff$replicate[ff$R1vR2=="1"]),
               fastq_1=ff$V1[ff$R1vR2=="1"],
               fastq_2=ff$V1[ff$R1vR2=="2"],
               strandedness="auto",
               sample_name=ff$genotype[ff$R1vR2=="1"],
               replicate=ff$replicate[ff$R1vR2=="1"],
               strain=ff$sample[ff$R1vR2=="1"],
               sample_acc=NA)

write.csv(df,paste0(workDir,"/samplesheet_SH.csv"),row.names=F,
          quote=F)


ss1<-read.csv(paste0(workDir,"/samplesheet.csv"))
ss2<-read.csv(paste0(workDir,"/samplesheet_SH.csv"))

ss<-rbind(ss1,ss2)

write.csv(ss,file=paste0(workDir,"/samplesheet_all.csv"),row.names=F,
            quote=F)

contrasts<-data.frame(id=c(paste0(c("Btm", "EM38", "hpl2tm","lin61tm"),"_vs_N2"),
                           paste0(c("EM90","EM91",
                                    "EM92"),"_vs_EM88")),
                      variable="strain",
                      reference=c(rep("N2",4),rep("EM88",3)),
                      target=c(c("Btm", "EM38", "hpl2tm","lin61tm"),
                               c("EM90","EM91","EM92")),
                      blocking="replicate")

write.table(contrasts,file="./contrasts.csv",sep=",",row.names=F, quote=F)


