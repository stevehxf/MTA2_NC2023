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
data <- read.table("data/3D_cleancounts.txt", header=T, row.names=1, check.names = F)  # Put your unnormalized counts file in data directory



# extract protein_coding gene counts from the raw counts
coding_genes <- grcm38 %>% dplyr::filter(biotype == "protein_coding")
coding_genes <- coding_genes$symbol

idx <- row.names(data) %in% coding_genes
data <- data[idx,]
# load meta
meta <- read.csv("meta/meta.csv", row.names=1, header=T) # Put your meta file in meta directory
# subset data
data <- data[,rownames(meta)]

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

## pre-processing the data: filter and cleaning ----
# pre-filtering: only keep the genes that at least 3 (n.rep.min) samples have more than 5 counts
n.min <- 2
keep <- rowSums((data)>5) >= n.min
data <- data[keep,]

# create DESeq2Dataset object. Count normalization for QC, etc.
dds <- DESeqDataSetFromMatrix(countData = data,  
                              colData = meta,
                              design = ~genotyping)  # Comparing factor 


# Transform counts for data visualization
vsd <- vst(dds, blind=F)
plotPCA(vsd, intgroup="genotyping")
vsd.df <- vsd %>% assay() %>% data.frame()

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
## use coefficient
# Extract DE_table without or with shrunken
name.list <- resultsNames(dds)[-1] %>% as.list()
names(name.list) <- resultsNames(dds)[-1] 
name.list
res.list <- name.list
resLFC.list <- name.list
resLFC.tb.list <- name.list

resLFC.apeglm.list <- name.list
resLFC.apeglm.tb.list <- name.list
resLFC.normal.list <- name.list
resLFC.normal.tb.list <- name.list

## use coefficient
for (i in seq(res.list)) {
        res.list[[i]] <- results(dds, name = name.list[[i]][1], alpha = 0.1)  # defaut alpha = 0.1
        
        resLFC.normal.list[[i]] <- lfcShrink(dds, coef = name.list[[i]][1], res=res.list[[i]], type = 'normal')
        resLFC.normal.tb.list[[i]] <- resLFC.normal.list[[i]] %>% res2tib()
        resLFC.normal.tb.list[[i]] %>% write_csv(paste0("results/DE_table/resLFC_normal_",name.list[[i]][1],".csv"))
        
        resLFC.apeglm.list[[i]] <- lfcShrink(dds, coef = name.list[[i]][1], res=res.list[[i]], type = 'apeglm')
        resLFC.apeglm.tb.list[[i]] <- resLFC.apeglm.list[[i]] %>% res2tib()
        resLFC.apeglm.tb.list[[i]] %>% write_csv(paste0("results/DE_table/resLFC_apeglm_",name.list[[i]][1],".csv"))
        
}

resLFC.apeglm.tb.list %>% saveRDS("rds/resLFC.apeglm.tb.list.rds")
resLFC.normal.tb.list %>% saveRDS("rds/resLFC.normal.tb.list.rds")
#resLFC.normal.tb.list <- readRDS("rds/resLFC.normal.tb.list.rds")

# Subset the table to only keep the significant genes using our pre-defined thresholds
## write a function
sigDE <- function(res_table_tb){
        res_table_tb %>% 
                dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) # filter the genes by thresholds
}
# threshold
padj.cutoff <- 0.05
lfc.cutoff <- 1
# subset DEG
DEG.list <- resLFC.normal.tb.list %>% lapply(sigDE)

##################### Data Visualization ######################
# convert data to rpkm 
source("scripts/ToRPKM.r")
rpkm <- data %>% ToRPKM("results/geneLength.csv")
## remove genes that are of low RPKM
filter.min <- function(data){
        n.min <- 2
        count.cutoff <- 1
        keep <- rowSums((data)>count.cutoff) >= n.min
        data <- data[keep,]
        return(data)
}
rpkm.filtered <- rpkm %>% filter.min()

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
table_tb <- resLFC.normal.tb.list$genotyping_hnf4a_oe_vs_ctrl
idx <- which(-log10(table_tb$padj) > ceiling)
table_tb[idx,]$padj <- 10^(-ceiling)
EnhancedVolcano(table_tb,
                lab = table_tb$gene,
                x = 'log2FoldChange',
                y = 'padj',    
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                #ylim = c(0, -log10(10e-50)),
                labSize = 6.0)
ggsave("figures/vol.hnf4avsctrl.pdf")

table_tb <- resLFC.normal.tb.list$genotyping_hnf4g_oe_vs_ctrl
idx <- which(-log10(table_tb$padj) > ceiling)
table_tb[idx,]$padj <- 10^(-ceiling)
EnhancedVolcano(table_tb,
                lab = table_tb$gene,
                x = 'log2FoldChange',
                y = 'padj',    
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                #ylim = c(0, -log10(10e-50)),
                labSize = 6.0)
ggsave("figures/vol.hnf4gvsctrl.pdf")

## mta2signature 
# heatmap
res.mta2 <- read.csv("../20220805_Mta2_2batches/results/DE_table/resLFC_normal_mta2_genotyping_MTA2_ko_vs_Ctrl.csv")
mta2ko.signature <- res.mta2 %>% filter(padj < 0.05 & log2FoldChange > 1)
idx <- rownames(vsd.df) %in% mta2ko.signature$gene
p <- vsd.df[idx,1:4] %>% plotHM()
plot_grid(p$gtable)
ggsave("figures/heatmap_hnf4aVSctrl_mta2kosignature.pdf", width = 8, height = 10)

# gsea
gsea.term <- data.frame(gs_name="mta2ko.signature", mta2ko.signature=mta2ko.signature$gene)
## fold change sorting ----
# hnf4a vs ctrl
fc <- resLFC.normal.tb.list$genotyping_hnf4a_oe_vs_ctrl$log2FoldChange
names(fc) <- resLFC.normal.tb.list$genotyping_hnf4a_oe_vs_ctrl$gene
fc <- sort(fc, decreasing = TRUE)

## run GSEA using germfree or germctrl signature on mta2 dataset
hnf4a.gsea <- GSEA(fc, TERM2GENE=gsea.term)
hnf4a.gsea@result %>% write.csv("results/functional_analysis/gsea/gsea_hnf4aVSctrl_mta2kosignature.csv")
gseaplot2(hnf4a.gsea, geneSetID = 1)
ggsave("figures/gsea_hnf4aVSctrl_mta2kosignature.pdf", width = 8, height = 10)

## mta2signature with rpkm filtered genes
# heatmap
res.mta2 <- read.csv("../20220805_Mta2_2batches/results/DE_table/resLFC_normal_mta2_genotyping_MTA2_ko_vs_Ctrl.csv")
mta2ko.signature <- res.mta2 %>% filter(padj < 0.05 & log2FoldChange > 1)

vsd.df.rpkm <- vsd.df[rownames(rpkm.filtered),]
idx <- rownames(vsd.df.rpkm) %in% mta2ko.signature$gene
vsd.df.rpkm[idx,1:4] %>% plotHM()
# gsea
gsea.term <- data.frame(gs_name="mta2ko.signature", mta2ko.signature=mta2ko.signature$gene)
## fold change sorting
# hnf4a vs ctrl

fc <- resLFC.normal.tb.list$genotyping_hnf4a_oe_vs_ctrl$log2FoldChange
names(fc) <- resLFC.normal.tb.list$genotyping_hnf4a_oe_vs_ctrl$gene
fc <- sort(fc, decreasing = TRUE)
idx <- names(fc) %in% rownames(rpkm.filtered)
fc.rpkmfiltered <- fc[idx]
## run GSEA using germfree or germctrl signature on mta2 dataset
hnf4a.gsea <- GSEA(fc.rpkmfiltered, TERM2GENE=gsea.term)
gseaplot2(hnf4a.gsea, geneSetID = 1)




