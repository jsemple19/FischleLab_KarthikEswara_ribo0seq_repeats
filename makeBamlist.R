library(dplyr)
library(rtracklayer)

workDir="/Volumes/meister.data/FischleLab_KarthikEswara/ribo0seq"
setwd(workDir)
bl<-read.csv(paste0(workDir,"/bamlist.txt"),header=F)
sampleName<-sapply(strsplit(basename(bl$V1),"\\."),"[[",1)
df<-data.frame(sampleName=sampleName,
               bamFile=bl$V1,
              replicate=sub(".*_(Rep[0-9]+)", "\\1", sampleName),
              group=sub("(.*)_Rep[0-9]+", "\\1", sampleName))
df

write.csv(df, file=paste0(workDir,"/bamlist.csv"),row.names=F,quote=F)

