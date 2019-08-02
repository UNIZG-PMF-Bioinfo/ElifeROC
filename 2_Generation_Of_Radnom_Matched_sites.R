#1. Create a table with 1kb regions on the genome. 
#2. Remove those which are blacklisted.
#3. Add DistanceToNearestGene

library(data.table)
library(GenomicRanges)
library(biomaRt)

genomeTiles <- readRDS("GenomeTiles.RDS")
genomeTiles <- as.data.table(genomeTiles)
genomeTiles[,strand:=NULL]
# distance to nearest protein coding gene is shown in kb:
genomeTiles[,distanceToNearestProteinCoding:=distanceToNearestProteinCoding%/%1000]

##3.3. adding epigenetic scores to tiles:
genomeTiles[,region_id:=paste(seqnames,start%/%1000,sep="_")]
Epigens <- readRDS("EpigenValues.RDS")
setkey(genomeTiles, region_id)
setkey(Epigens, region_id)
genomeTiles <- merge(genomeTiles, Epigens, all.x=T) # make sure that regions which do not appear in Epigens table have score 0!!!
genomeTiles[is.na(genomeTiles)]<-0

