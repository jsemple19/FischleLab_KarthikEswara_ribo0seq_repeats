library(stringr)
library(dplyr)
library(tidyr)

options(width=10000)
serverPath="/Volumes/meister.data"

# summarise annotation table from dfam3.5
dfam<-readRDS(paste0(serverPath,"/publicData/genomes/dfam35/repeats_ce11_Dfam_3.5_nr.rds"))
dfamFam<-data.frame(dfam) %>% group_by(source,Name,Alias) %>% summarise(rptAnnCount=n(),avgLength=mean(width))

dfamFam

dfamClass<-read.delim(paste0(serverPath,"/publicData/genomes/dfam40/dfam4.0_family_summary_6239.tsv"),skip=30,sep=" ",header=F, strip.white=T)
colnames(dfamClass)<-c("dfamID","name","classification","consensusLength")

dfamClass$name<-dfamClass$name |> str_replace(":","") |> str_replace_all("'","")
dfamClass$consensusLength<-dfamClass$consensusLength |> str_replace("len=","") |> as.numeric()
dfamClass$type<-sapply(strsplit(dfamClass$classification, ";"), tail, 1)

dfamClass$type[dfamClass$type=="Unknown"]<-dfamClass$name[dfamClass$type=="Unknown"]

dfamClass$class<-dfamClass$classification |> str_extract("Class_I[^;]*") |> str_replace("ransposition","ransposon")

satellites<-dfamClass$classification |> str_extract("Satellite[^;]*")

dfamClass$class[ifelse(!is.na(satellites),T,F)]<-satellites[ifelse(!is.na(satellites),T,F)]

dfamClass$type[dfamClass$type=="Transposase"]<-dfamClass$name[dfamClass$type=="Transposase"]

#dfamClass$class[is.na(dfamClass$class)]<-dfamClass$name[is.na(dfamClass$class)]
head(dfamClass)
tail(dfamClass)

dfamFam<-left_join(dfamFam,dfamClass,by=join_by("Name"=="name"))

oldIDs<-dfamFam$dfamID |> str_replace("DF00","DF") |> str_replace("\\..?","")
dfamFam$Alias==oldIDs

colnames(dfamFam)[c(2,3,6)]<-c("name","Dfam3.5_ID","Dfam4.0_ID")
dfamFam

write.table(dfamFam,paste0(serverPath,"/publicData/genomes/dfam35/repeatFamilies.tsv"),row.names=F,quote=F,sep="\t")
