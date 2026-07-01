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

# # rename gene_id with unique ids so it doesn't give error for attribute
# gtf<-import("/Volumes/meister.data/publicData/genomes/dfam35/WS295_canonicalgenes_Dfam3.5nr_repeats.gtf")
#
# gtf_genes<-gtf[gtf$source!="Dfam_3.5"]
#
# rpt_exon<-gtf[gtf$source=="Dfam_3.5" & gtf$type=="exon"]
# rpt_txpt<-gtf[gtf$source=="Dfam_3.5" & gtf$type=="transcript"]
# rpt_gene<-gtf[gtf$source=="Dfam_3.5" & gtf$type=="transcript"]
# rpt_gene$type<-"gene"
#
# rptgtf<-sort(c(rpt_exon,rpt_txpt,rpt_gene))
# rptgtf$gene_id<-rptgtf$transcript_id
#
# gtf_all<-c(gtf_genes,rptgtf)
#
# export(gtf_all,paste0(workDir,"/WS295_canonicalgenes_Dfam3.5nr_repeats.uniq_gene_id.gtf"))
#
#
#
#
# # contrasts<-read.csv(paste0(workDir,"/contrasts.csv"))
# # contrasts$treatmentBams<-NA
# # contrasts$controlBams<-NA
# # i=1
# # for(i in 1:nrow(contrasts)){
# #   contrastName<-contrasts[i,"id"]
# #   print(contrastName)
# #   contrasts$treatmentBams[i]<-paste(df$bamFile[df$sampleGroup==contrasts[i,"target"]],collapse=" ")
# #   contrasts$controlBams[i]<-paste(df$bamFile[df$sampleGroup==contrasts[i,"reference"]],collapse=" ")
# # }
# #
# # write.csv(contrasts[,c("id","treatmentBams","controlBams")],
# #             file=paste0(workDir,"/bamlist.csv"),row.names=F,quote=F)
