library(stringr)

setwd("/Volumes/meister.data/FieschleLab_KarthikEswara/ribo0seq")
# download metadata from ENA PRJEB102446
df<-read.delim(paste0(workDir,"/filereport_read_run_PRJEB102446.tsv"))

head(df)
write.table(df$run_accession,"ids.txt",row.names=F, col.names=F,quote=F)


