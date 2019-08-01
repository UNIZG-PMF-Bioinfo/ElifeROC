library(rtracklayer)
library(data.table)
duljinekromosoma <- readRDS("duljinekromosoma.RDS")

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
                           
