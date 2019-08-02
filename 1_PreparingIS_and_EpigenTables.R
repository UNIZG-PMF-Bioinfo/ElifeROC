library(rtracklayer)
library(data.table)
library(parallel) #needs a lot of RAM, 200G is recommended, 80 not enough, run withsouth parallel option with less RAM
library(biomaRt)
library(stringr)
duljinekromosoma <- readRDS("duljinekromosoma.RDS")

# Preparing the table with values of all epigenetic marks:
# I will create a table with 1kb tiles with values of epigenetic marks on those tiles. 
# This table will exclude blacklisted regions.
# I will also add distance to nearest protein coding gene to each tile. 
# This table will be used later to create random matched controls and also to 
# add scores to integration sites based on tiles the integrations are in.

#1. Create a table with 1kb regions on the genome. 
cs <- readRDS("duljinekromosoma.RDS")
genomeTiles <- unlist(tileGenome(cs, tilewidth=1000))

#2. Remove those which are blacklisted.
blacklisted <- fread("hg38-blacklist.v2.bed")
blacklistedRegions <- GRanges(blacklisted$V1, IRanges(blacklisted$V2, blacklisted$V3), reason=blacklisted$V4)
rm(blacklisted)
genomeTiles <- genomeTiles[countOverlaps(genomeTiles, blacklistedRegions, ignore.strand=T)==0]

#3. Add DistanceToNearestGene
## 3.1. get genes from biomaRt

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
atts <- c("ensembl_gene_id","chromosome_name","start_position","end_position","strand","gene_biotype")
allGenes <- getBM(attributes=atts, mart=ensembl)
colnames(allGenes) <- c("ensembl_gene_id", "chr","start","end","strand","biotype")
allGenes$chr <- paste("chr", as.character(allGenes$chr), sep="")
allGenes$strand <- ifelse(allGenes$strand=="1","+",ifelse(allGenes$strand=="-1","-","*"))
allGenes <- GRanges(allGenes)
proteinCodingGenes <- allGenes[allGenes$biotype=="protein_coding"]
proteinCodingGenes <- proteinCodingGenes[seqnames(proteinCodingGenes)%in%names(cs)]

##3.2. Add DistanceToNearestGene to genomeTiles:
genomeTiles$distanceToNearestProteinCoding <- data.frame(distanceToNearest(genomeTiles, proteinCodingGenes, ignore.strand=T), 
                                                          stringsAsFactors=F)$distance
genomeTiles <- as.data.table(genomeTiles)
genomeTiles <- genomeTiles[,c(1:3,6)]
genomeTiles[,region_id:=paste(seqnames,start%/%1000,sep="_")]
genomeTiles<-genomeTiles[,.SD[1],region_id] # to make sure region_id is unique
setkey(genomeTiles,seqnames,start,end)
# Preparing the table with values of all epigenetic marks:
# Marks are saved in files named Samplename_Epigenmark_sequence.bw


MakeAnEpigenScoreTable <- function(test_bw){
  gr <- import(test_bw)
  gr<-gr[seqnames(gr)%in%names(duljinekromosoma)]
  ddt <- as.data.table(gr)
  # find overlaps of 1kb tiles with the bigWig epigenomic score track:
  # some regions in bw file are longer than 1kb. They will overlap multiple intervals. 
  # Some are smaller. Multiple of those will overlap one tile. So I will calculate the average score
  # on each tile.
  setkey(ddt,seqnames,start,end)
  genomeTiles <- foverlaps(genomeTiles,ddt)
  genomeTiles[,overlapStart:=ifelse(i.start>start,i.start,start),]
  genomeTiles[,overlapEnd:=ifelse(i.end>end,end,i.end),]
  genomeTiles[,overlapWidth:=overlapEnd-overlapStart+1,]
  genomeTiles[,Score:=sum(overlapWidth*score)/sum(overlapWidth),region_id]
  
  
  genomeTiles <- genomeTiles[,.N,.(region_id,distanceToNearestProteinCoding,Score)][,.(region_id,distanceToNearestProteinCoding,Score)]
  genomeTiles[,Epigen:=str_extract(test_bw, ".*(?=(_sequence))")]
  genomeTiles[,.(Epigen,region_id,Score,distanceToNearestProteinCoding)]
}

allEpigenomeTables<-mclapply(list.files(pattern=".bw"), function(x)MakeAnEpigenScoreTable(x), mc.cores=10)
# The option that runs slower but consumes less RAM:
#allEpigenomeTables<-lapply(list.files(pattern=".bw"), function(x)MakeAnEpigenScoreTable(x))
names(allEpigenomeTables) <- str_extract(list.files(pattern=".bw"), ".*(?=(_sequence))")
allEpigenomeTables <- do.call("rbind",allEpigenomeTables)
Epigens <- dcast(allEpigenomeTables, region_id + distanceToNearestProteinCoding ~ Epigen,value.var="Score")

#Preparing the Integration sites table:
bix<-fread("BIXhg38.bed")
bix[,bin:=V2%/%1000]
bix[,region:=paste(V1,bin, sep="_")]
colnames(bix)[7] <- "region_id"

setkey(Epigens,region_id)
setkey(bix,region_id)

bixepi<-merge(bix,Epigens, all.x=T)
bixepi<-bixepi[!is.na(bixepi$D75_Ledgf)] # excluding blacklisted integration sites!

saveRDS(bixepi, file="BIX_hg38_withEpigenValues.RDS")

                             
bixepi[,dist:=distanceToNearestProteinCoding%/%1000]
Epigens[,dist:=distanceToNearestProteinCoding%/%1000]
bixepisplit<-split(bixepi,bixepi$V5)
thismanytotakeBIX <- table(bixepisplit[["BIX"]]$dist)
thismanytotakeCONTROL <- table(bixepisplit[["CONTROL"]]$dist)
Epigens$takeThisManyBIX <- 0
Epigens$takeThisManyCONTROL <- 0
TAKETHATBIX<-thismanytotakeBIX[as.character(Epigens[dist%in%as.numeric(names(thismanytotakeBIX)),dist])]
Epigens[dist%in%as.numeric(names(thismanytotakeBIX)),"takeThisManyBIX"] <- TAKETHATBIX
TAKETHATCONTROL<-thismanytotakeCONTROL[as.character(Epigens[dist%in%as.numeric(names(thismanytotakeCONTROL)),dist])]
Epigens[dist%in%as.numeric(names(thismanytotakeCONTROL)),"takeThisManyCONTROL"] <- TAKETHATCONTROL

saveRDS(Epigens, "genomeTilesWithAddedEpigenScores.RDS")
set.seed(5555)
BIXmatched <- replicate(10,Epigens[, .SD[sample(.N, takeThisManyBIX, replace = T)], dist][order(-dist)])
CONTROLmatched <- replicate(10,Epigens[, .SD[sample(.N, takeThisManyCONTROL, replace = T)], dist][order(-dist)])


doItForCONTROLEpigen <- function(epiname){
  mean(do.call("cbind",CONTROLmatched[epiname,])>unlist(bixepisplit[["CONTROL"]][order(-dist)][,..epiname]))
}

doItForBIXEpigen <- function(epiname){
  mean(do.call("cbind",BIXmatched[epiname,])>unlist(bixepisplit[["BIX"]][order(-dist)][,..epiname]))
}
epinames <- c("D75_H3K36me3","D75_Input","D75_bixH3K36me3","D75_bixInput","D86_H3K36me3", "D86_Input","D86_bixH3K36me3","D86_bixInput","D75_H3K9me2","D75_Ledgf","D75_bixH3K9me2","D75_bixledgf","D86_H3K9me2","D86_Ledgf","D86_bixH3K9me2","D86_bixledgf" )
CTRL<-sapply(epinames,doItForCONTROLEpigen)
BIX<-sapply(epinames,doItForBIXEpigen)
x <- data.table(CTRL,BIX)
#####################################

# WHEN NOT MATCHED BY DIST TO NEAREST GENE:
BIXNOTmatched <- replicate(10,Epigens[, .SD[sample(.N, 635, replace = T)], dist][order(-dist)])
CONTROLNOTmatched <- replicate(10,Epigens[, .SD[sample(.N, 449, replace = T)], dist][order(-dist)])
UNMAdoItForCONTROLEpigen <- function(epiname){
  mean(do.call("cbind",CONTROLNOTmatched[epiname,])>unlist(bixepisplit[["CONTROL"]][order(-dist)][,..epiname]))
}

UNMAdoItForBIXEpigen <- function(epiname){
  mean(do.call("cbind",BIXNOTmatched[epiname,])>unlist(bixepisplit[["BIX"]][order(-dist)][,..epiname]))
}
epinames <- c("D75_H3K36me3","D75_Input","D75_bixH3K36me3","D75_bixInput","D86_H3K36me3", "D86_Input","D86_bixH3K36me3","D86_bixInput","D75_H3K9me2","D75_Ledgf","D75_bixH3K9me2","D75_bixledgf","D86_H3K9me2","D86_Ledgf","D86_bixH3K9me2","D86_bixledgf" )
UNMACTRL<-sapply(epinames,UNMAdoItForCONTROLEpigen)
UNMABIX<-sapply(epinames,UNMAdoItForBIXEpigen)
xUNMA <- data.table(UNMACTRL,UNMABIX)


