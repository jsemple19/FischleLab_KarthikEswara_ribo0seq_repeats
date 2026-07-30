
workDir<-"/Volumes/meister.data/FischleLab_KarthikEswara/ribo0seq"
setwd(workDir)
dir.create(paste0(workDir,"/forGEO"))

## canonocal genes
cnts<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv"))

cnts<-cnts[!grepl("_rpt[[:digit:]]*$", cnts$gene_id),]

idx1<-colnames(cnts) %in% c("Btm_Rep1", "Btm_Rep2", "Btm_Rep3", "EM38_Rep1", "EM38_Rep2", "EM38_Rep3", "N2_Rep1", "N2_Rep2", "N2_Rep3", "hpl2tm_Rep1", "hpl2tm_Rep2", "hpl2tm_Rep3", "lin61tm_Rep1", "lin61tm_Rep2", "lin61tm_Rep3")

idx2<-colnames(cnts) %in% c("EM88_Rep1", "EM88_Rep2", "EM88_Rep3", "EM90_Rep1", "EM90_Rep2", "EM90_Rep3", "EM92_Rep1", "EM92_Rep2", "EM92_Rep3", "EM91_Rep1", "EM91_Rep2", "EM91_Rep3")

salmon1<-cnts[,c("gene_id","gene_name",colnames(cnts)[idx1])]
write.table(salmon1,paste0(workDir,"/forGEO/canonical_genes_star_salmon_counts_set1.tsv"),
            sep="\t",row.names=F,quote=F)

salmon2<-cnts[,c("gene_id","gene_name",colnames(cnts)[idx2])]
write.table(salmon2,paste0(workDir,"/forGEO/canonical_genes_star_salmon_counts_set2.tsv"),
            sep="\t",row.names=F,quote=F)

### repeats
cnts<-read.delim(paste0(workDir,"/results_telescope/txt/allCounts.tsv"),sep=" ")

cnts<-cnts[grepl("_rpt[[:digit:]]*$", cnts$gene_id),]

idx1<-colnames(cnts) %in% c("Btm_Rep1", "Btm_Rep2", "Btm_Rep3", "EM38_Rep1", "EM38_Rep2", "EM38_Rep3", "N2_Rep1", "N2_Rep2", "N2_Rep3", "hpl2tm_Rep1", "hpl2tm_Rep2", "hpl2tm_Rep3", "lin61tm_Rep1", "lin61tm_Rep2", "lin61tm_Rep3")

idx2<-colnames(cnts) %in% c("EM88_Rep1", "EM88_Rep2", "EM88_Rep3", "EM90_Rep1", "EM90_Rep2", "EM90_Rep3", "EM92_Rep1", "EM92_Rep2", "EM92_Rep3", "EM91_Rep1", "EM91_Rep2", "EM91_Rep3")

telescope1<-cnts[,c("gene_id",colnames(cnts)[idx1])]
write.table(telescope1,paste0(workDir,"/forGEO/repeats_star_telescope_counts_set1.tsv"),
            sep="\t",row.names=F,quote=F)

telescope2<-cnts[,c("gene_id",colnames(cnts)[idx2])]
write.table(telescope2,paste0(workDir,"/forGEO/repeats_star_telescope_counts_set2.tsv"),
            sep="\t",row.names=F,quote=F)
