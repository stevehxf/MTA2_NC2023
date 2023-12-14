## This script take raw.count as input and output rpkm normalized counts
## Original count table was expected to have gene symbole as row names.

require(tidyverse)
require(data.table)
require(GenomicFeatures)
require(annotables)
require(edgeR)



ToRPKM <- function(data.raw, path.genelengths){
  
  
# read gtf file to calculate gene length
# first, import the GTF-file that you have also used as input for htseq-count
#txdb <- makeTxDbFromGFF("./data/GRCm38.M18/gencode.vM18.primary_assembly.annotation.gtf",format="gtf")
# then collect the exons per gene id
#exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
#exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
#unlist_geneLength <- unlist(exonic.gene.sizes)
#write.csv(unlist_geneLength, file = "results/geneLength.csv")
genelengths<-read.csv(file = path.genelengths)
colnames(genelengths) <- c("gene","length")

# strip the version value from the ID
genelengths$gene <- nth(tstrsplit(genelengths$gene, split ="\\."),n=1)

# Match ID with gene symbols and remove duplicates

## Return the IDs for the gene symbols 
ids <- grcm38[grcm38$symbol %in% rownames(data.raw), ]
## The gene names can map to more than one Ensembl ID (some genes change ID over time), so we need to remove duplicate IDs. Duplicated returns a logical vector indicating which rows of a data.table are duplicates of a row with smaller subscripts.
non_duplicates <- which(duplicated(ids$symbol) == FALSE)
ids <- ids[non_duplicates, ] 

## Create tibbles including row names
data.raw.tb <- data.raw %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble() # meta data in tibble format

## Merge the IDs with the results 
data.raw.tb <- inner_join(data.raw.tb, ids[,c(1,3)], by=c("gene"="symbol")) # merge
data.raw.tb <- dplyr::select(data.raw.tb, gene, ensgene, everything()) # order the colomn

## Merge the count table with the gene length based on ID 
data.raw.tb <- inner_join(data.raw.tb, genelengths, by=c("ensgene"="gene")) # merge
data.raw.tb <- dplyr::select(data.raw.tb, gene, ensgene, length, everything()) # order the colomn

## Compute rpkm matrix
data.raw.count.tb <- dplyr::select(data.raw.tb, -(gene:length))
gene.lengths <- data.raw.tb$length
rpkm.matrix<-rpkm(data.raw.count.tb,gene.length = gene.lengths)
rpkm.df <- as.data.frame(rpkm.matrix)
row.names(rpkm.df) <- data.raw.tb$gene
return(rpkm.df)

}