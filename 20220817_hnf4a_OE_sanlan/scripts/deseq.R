# This script is designed for differential expression (DE) analysis using DEseq2 as well as data visualization
# apeglm shrunkage used in this script

# Load libraries
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(DESeq2)
library(pheatmap)
library(ggrepel)
#library(DEGreport)
#library(BiocParallel)
#library(apeglm)
library(ggfortify)
library(data.table)
library(limma)
#library(GenomicFeatures)
library(annotables)
library(edgeR)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
#library(DOSE)
#library(pathview)
#library(purrr)
library(clusterProfiler)
library(enrichplot)
library(EnhancedVolcano)
#library(Rgraphviz)
library(cowplot)

# load feature counts
data <- read.table("data/hnf4a_cleancounts.txt", header=T, row.names=1, check.names = F)  # Put your unnormalized counts file in data directory



# extract protein_coding gene counts from the raw counts
coding_genes <- grch38 %>% dplyr::filter(biotype == "protein_coding")
coding_genes <- coding_genes$symbol

idx <- row.names(data) %in% coding_genes
data <- data[idx,]
# load and subset meta
meta <- read.csv("meta/meta.csv", row.names=1, header=T) # Put your meta file in meta directory

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

## pre-processing the data: filter and cleaning ----
# pre-filtering: only keep the genes that at least 3 (n.rep.min) samples have more than 5 counts
n.min <- 3
keep <- rowSums((data)>5) >= n.min
data <- data[keep,]

# create DESeq2Dataset object. Count normalization for QC, etc.
dds <- DESeqDataSetFromMatrix(countData = data,  
                              colData = meta,
                              design = ~donor+genotyping)  # Comparing factor 


# Transform counts for data visualization
vsd <- vst(dds, blind=F)
vsd.batchRM <- vsd
assay(vsd.batchRM) <- dds %>% vst(blind=T) %>% assay() %>% limma::removeBatchEffect(vsd$donor)

plotPCA(vsd, intgroup="genotyping")
plotPCA(vsd.batchRM, intgroup="genotyping")

vsd.df <- vsd %>% assay() %>% data.frame()
vsd.batchRM.df <- vsd.batchRM %>% assay() %>% data.frame()

## DE analysis ----
# change reference level
dds$genotyping <- dds$genotyping %>% relevel(ref = "ctrl")


## make directory for data storing
if(!dir.exists("results/DE_table")){
        dir.create("results/DE_table")}
if(!dir.exists("results/normalization")){
        dir.create("results/normalization")}
if(!dir.exists("dds")){
        dir.create("dds")}

## Convert the result table into a tibble
## write a function
res2tib <- function(result_table) {
        result_table %>% 
                data.frame() %>% 
                rownames_to_column(var="gene") %>% #convert rownames to a column
                as_tibble() # covert the dataframe to tibble
}  ## res2tib USAGE: res2tib(result_table)

# Extract normalized counts, check the size factors and total counts before and after normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds) # check the size factors
normalized_counts <- counts(dds, normalized=TRUE)
dir <-"results/normalization/"
# export normalized counts
write.table(normalized_counts, 
            file=paste0(dir, "normalized_counts.txt"),
            sep="\t", quote=F, col.names=NA)
## Take a look at the total number of reads for each sample using
colSums(counts(dds)) # Sum the unnormalized column values
colsum <- colSums(counts(dds)) 
write.csv(colsum, 
          file=paste0(dir, "colSum.csv") # Export
)
colSums(normalized_counts) # Sum the normalized column values
ncolsum <- colSums(normalized_counts) 
write.csv(ncolsum, 
          file=paste0(dir, "ncolSum.csv")
)

# Differential expression analysis
dds <- DESeq(dds)
sizeFactors(dds) %>% print()

# Examine dispersion estimates plot to ensure your data is a good fit for the DESeq2 model
plotDispEsts(dds) 

# Create contrasts to perform Wald testing on the shrunken log2 foldchanges between specific conditions:
# check the coef for shurnkage
resultsNames(dds) %>% print()
# Extract DE_table without or with shrunken
saveRDS(dds, "dds/dds.rds")
dds <- readRDS("dds/dds.rds")
## use coefficient
name <- resultsNames(dds)[3]
res <- results(dds, name = name[1], alpha = 0.1)  # defaut alpha = 0.1
resLFC.normal <- lfcShrink(dds, coef = name[1], res=res, type = 'normal')
resLFC.normal.tb <- resLFC.normal %>% res2tib()
resLFC.normal.tb %>% write_csv(paste0("results/DE_table/resLFC_normal_", name[1],".csv"))

resLFC.apeglm <- lfcShrink(dds, coef = name[1], res=res, type = 'apeglm')
resLFC.apeglm.tb <- resLFC.apeglm %>% res2tib()
resLFC.apeglm.tb %>% write_csv(paste0("results/DE_table/resLFC_apeglm_", name[1],".csv"))


# Subset the table to only keep the significant genes using our pre-defined thresholds
## write a function
sigDE <- function(res_table_tb){
        res_table_tb %>% 
                dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) # filter the genes by thresholds
}
# threshold
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
# subset DEG
DEG <- resLFC.normal.tb %>% sigDE()
View(DEG)

##################### Data Visualization ######################

# 1.1 filter with RPKM
#rpkm.df <- read.table("results/rpkm.txt", sep="\t", row.names = 1, header = T)
## set up filter threshold
#rpkm.filter <- 8 # rpkm threshold
#sample.filter <- 2 # sample number threshold
## filter
#keep <- rowSums(rpkm.df > rpkm.filter) >= sample.filter 
#rpkm.df.filtered <- rpkm.df[keep,]
#idx <- rownames(vsd.df) %in% rownames(rpkm.df.filtered)
#vsd.df.rpkm.filtered <- vsd.df[idx,]

# 2. Heatmap (Filter with RPKM counts)
# write a function to plot heatmap using the same arguments
plotHM <- function(data,title=""){
        pheatmap(data.frame(data), 
                 main = title,
                 #color = colorRampPalette(c("cornflowerblue","beige","brown1"))(50),
                 #color = magma(50),
                 color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(50),
                 #breaks = seq(-1.5,2.5,by=0.1),
                 scale="row", 
                 cluster_rows = T, cluster_cols = T,
                 show_colnames = T, show_rownames = T,
                 border_color = NA, 
                 fontsize = 10, 
                 fontsize_row = 10, 
                 #cellwidth = 10,
                 #cellheight= 10
                 )
}

# heatmap

p <- vsd.batchRM.df[DEG$gene,] %>% plotHM()
plot_grid(p$gtable)
ggsave("figures/heatmap_hnf4aOEVSctrl.pdf", width = 8, height = 20)


### DEG quantification and visualization
## plot with fixed gene list
#Absorption.df <- read.table("results/absorption.txt", header = T)
#feature.to.plot <- Absorption.df$gene
#sample.to.plot <- 
## which row exist in the given gene list?
#idx <- rownames(vsd.batchRM) %in% feature.to.plot
#vsd.batchRM.df.sub <- vsd.batchRM.df[idx, ]
## plot a heatmap
#p1 <- vsd.batchRM.df.sub %>% plotHM()
#plot_grid(p1$gtable)
#ggsave("figures/heatmap_absorption.pdf")
# fat genes
#fat.genes <- "Mttp/Npc1l1/Abcg8/Pla2g2a/Abcg5/Pla2g5/Pla2g2f/Plpp1/Mogat2/Pla2g12b/Pla2g3/Pla2g10/Acat2" %>% str_split("/") %>% unlist()
#idx <- rownames(vsd.batchRM) %in% fat.genes
#vsd.batchRM.df.sub <- vsd.batchRM.df[idx, ]
#p1 <- vsd.batchRM.df.sub %>% plotHM()
#plot_grid(p1$gtable)
#ggsave("figures/fat.genes.pdf")

## volcano plot
## adjust ceiling value
ceiling <- 100
table_tb <- resLFC.normal.tb
idx <- which(-log10(table_tb$padj) > ceiling)
table_tb[idx,]$padj <- 10^(-ceiling)
EnhancedVolcano(table_tb,
                lab = rownames(resLFC.normal),
                x = 'log2FoldChange',
                y = 'padj',    
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                #ylim = c(0, -log10(10e-50)),
                labSize = 6.0)
ggsave("figures/vol.pdf")

