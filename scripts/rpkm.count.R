## This script take raw.count as input and output rpkm normalized counts
## Original count table was expected to have gene symbole as row names.


library(tidyverse)
library(dplyr)
library(data.table)
library(GenomicFeatures)
library(annotables)
library(edgeR)
# read the raw counts
data.raw <- read.table("data/data.txt", header=T, row.names=1, check.names = F)  # Put your unnormalized counts file in data directory
# read gtf file to calculate gene length
# first, import the GTF-file that you have also used as input for htseq-count
#txdb <- makeTxDbFromGFF("./data/GRCm38.M18/gencode.vM18.primary_assembly.annotation.gtf",format="gtf")
# then collect the exons per gene id
#exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
#exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
#unlist_geneLength <- unlist(exonic.gene.sizes)
#write.csv(unlist_geneLength, file = "results/geneLength.csv")
genelengths<-read.csv(file = "results/geneLength.csv")
colnames(genelengths) <- c("gene","length")

# strip the version value from the ID
genelengths$gene <- nth(tstrsplit(genelengths$gene, split ="\\."),n=1)

# Match ID with gene symbols and remove duplicates
## Explore the grcm38 table loaded by the annotables library
grcm38

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
rpkm.df
write.table(rpkm.df, "results/rpkm.txt", sep="\t", quote=F, col.names=NA)


plotHistFunc <- function(dataframe, pre, na.rm = TRUE, ...) {     # make a function
  nm <- colnames(dataframe)              # put column names in the variable nm
  for (i in seq_along(nm)) {            # for loop, extract the index number from nm
    plots <- ggplot(dataframe,aes_string(x = nm[i])) +    # ggplot(file,aes_from_string)
      geom_histogram(stat ="bin", bins = 200, fill = "#446179", alpha = 0.8) + # plot histogram
      xlim(-0.5, 2)  + # x axis limit
      ylim(0, 1000) +
      xlab("RPKM") + # x lable
      ylab("Number of genes") + # y lable
      theme_light() + # theme
      theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 1, colour = "black")) # border
    ggsave(plots,filename=paste("results/count_distribution/",pre,nm[i],".png",sep="")) # export
  }}

plotHistFunc(rpkm.df, "rpkm_")
