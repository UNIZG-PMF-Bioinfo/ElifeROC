#1. Create a table with 1kb regions on the genome. 
#2. Remove those which are blacklisted.
#3. Add DistanceToNearestGene

library(data.table)
library(GenomicRanges)
library(biomaRt)


#1. Create a table with 1kb regions on the genome. 
cs <- readRDS("duljinekromosoma.RDS")
genomeTiles <- unlist(tileGenome(cs, tilewidth=1000))

#2. Remove those which are blacklisted.
blacklisted <- fread("hg38-blacklist.v2.bed")
blacklistedRegions <- GRanges(blacklisted$V1, IRanges(blacklisted$V2, blacklisted$V3), reason=blacklisted$V4)
rm(blacklisted)
genomeTiles <- genomeTiles[countOverlaps(genomeTiles, blacklistedRegions, ignore.strand=T)>0]

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
saveRDS(genomeTiles, file="GenomeTiles.RDS")
