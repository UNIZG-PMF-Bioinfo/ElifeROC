library(rtracklayer)
library(data.table)
duljinekromosoma <- readRDS("duljinekromosoma.RDS")

#test_bw <- file.path("D75_bixH3K36me3_sequence.bw")

MakeAnEpigenScoreTable <- function(test_bw){
  gr <- import(test_bw)
  gr<-gr[seqnames(gr)%in%names(duljinekromosoma)]
  ddt <- as.data.table(gr)
  ddt[,region:=start%/%1000]
  ddt[,region_id:=paste(seqnames,region,sep="_")]
  ddt[,Score:=sum(score*width)/sum(width),by=region_id]
  ddt[,.N, .(region_id,Score)][,.(region_id,Score)]
}
#bixH3K36me3 <- MakeAnEpigenScoreTable(test_bw)

allEpigenomeTables<-lapply(list.files(pattern=".bw"), function(x)MakeAnEpigenScoreTable(x))
