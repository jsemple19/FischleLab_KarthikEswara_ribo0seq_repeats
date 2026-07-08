library(stringr)

setwd("/Volumes/meister.data/FieschleLab_KarthikEswara/ribo0seq")
# download metadata from ENA PRJEB102446
df<-read.delim("filereport_read_run_PRJEB102446.tsv")

head(df)
write.table(df$run_accession,"ids.txt",row.names=F, col.names=F,quote=F)

## sample sheet
oldfiles<-df$submitted_ftp |> str_split_i(";",1) |> str_split_i("\\/",6) |>
  str_replace("\\.fq\\.gz","")

workDir<-"/mnt/meister.data/jsemple/ribo0seq/"
samplesheet<-data.frame(sample= NA,
                        fastq_1=paste0(workDir,"/fastq/",df$experiment_accession,"_",df$run_accession,"_1.fastq.gz"),
           fastq_2=paste0(workDir,"/fastq/",df$experiment_accession,"_",df$run_accession,"_2.fastq.gz"),
      strandedness="auto",
      sampleName=oldfiles |> str_replace("_Rep.*",""),
      replicate=oldfiles |> str_extract("_Rep[:digit:]") |> str_replace("_",""))

