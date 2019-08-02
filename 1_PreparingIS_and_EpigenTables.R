library(rtracklayer)
library(data.table)
duljinekromosoma <- readRDS("duljinekromosoma.RDS")

# Preparing the table with values of all epigenetic marks:
# Marks are saved in files named Samplename_Epigenmark_sequence.bw
MakeAnEpigenScoreTable <- function(test_bw){
  gr <- import(test_bw)
  gr<-gr[seqnames(gr)%in%names(duljinekromosoma)]
  ddt <- as.data.table(gr)
  ddt[,region:=start%/%1000]
  ddt[,region_id:=paste(seqnames,region,sep="_")]
  ddt[,Score:=sum(score*width)/sum(width),by=region_id]
  ddt[,Epigen:=str_extract(test_bw, ".*(?=(_sequence))")]
  ddt[,.N, .(Epigen,region_id,Score)][,.(Epigen,region_id,Score)]
}

allEpigenomeTables<-lapply(list.files(pattern=".bw"), function(x)MakeAnEpigenScoreTable(x))
names(allEpigenomeTables) <- str_extract(list.files(pattern=".bw"), ".*(?=(_sequence))")
allEpigenomeTables <- do.call("rbind",allEpigenomeTables)
Epigens<- dcast(allEpigenomeTables, region_id ~ Epigen , value.var="Score", fill=0)
saveRDS(Epigens, "allEpigenomeTables.RDS")

#Preparing the Integration sites table:
bix<-fread("BIXhg38.bed")
bix[,bin:=V2%/%1000]
bix[,region:=paste(V1,bin, sep="_")]
colnames(bix)[7] <- "region_id"

setkey(Epigens,region_id)
setkey(bix,region_id)

bixepi<-merge(bix,epi, all.x=T)
saveRDS(bixepi, file="BIX_hg38_withEpigenValues.RDS")
