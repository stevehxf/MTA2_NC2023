# This script is designed for differential expression (DE) analysis using DEseq2 as well as data visualization
# apeglm shrunkage used in this script

######################## USAGE & SETUP ######################## 
## ***Please customize the content and run this part first ####

# 0. Create an R project and make a series of directories. 
# project_directory/data
# project_directory/meta
# project_directory/logs
# project_directory/results/count_distribution
# project_directory/results/meanVSvar
# project_directory/results/dispersion
# project_directory/results/DE_table
# project_directory/results/functional_analysis/ClusterProfiler
# project_directory/scripts
# project_directory/summary
# project_directory/results/DE_table/rpkm_filtered

# Load libraries
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(DEGreport)
library(BiocParallel)
library(apeglm)
library(ggfortify)
library(data.table)
library(limma)
#library(GenomicFeatures)
library(annotables)
#library(edgeR)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
#library(DOSE)
#library(pathview)
#library(purrr)
#library(clusterProfiler)
#library(annotables)
#library(enrichplot)
#library(Rgraphviz)
library(cowplot)
# load feature counts
data.raw <- read.table("data/satb2_cdx2_hnf4a_cleancounts_2.txt", header=T, row.names=1, check.names = F)  # Put your unnormalized counts file in data directory
# extract protein_coding gene counts from the raw counts
coding_genes <- grcm38 %>% dplyr::filter(biotype == "protein_coding")
coding_genes <- coding_genes$symbol
idx <- row.names(data.raw) %in% coding_genes
data <- data.raw[idx,]
data <- data[c(1:9)]
# load and subset meta
meta <- read.csv("meta/meta.csv", row.names=1, header=T) # Put your meta file in meta directory
idx <- row.names(meta) %in% colnames(data)
meta <- meta[idx,]
meta$batch <- meta$batch %>% as_factor()

# 3. Check the input files
##  (1) Check the data class to confirm they are "data.frame"
class(data)
class(meta)
##  (2) Check the colume names of data match the row names of meta: both should be TRUE
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))


# Genes that are barely detectable will contribute noise to the resulting clusters
# pre-filtering: only keep the genes that at least 3 (n.rep.min) samples have more than 5 counts
n.min <- 3
keep <- rowSums((data)>5) >= n.min
data <- data[keep,]
# create DESeq2Dataset object. Count normalization for QC, etc.
dds <- DESeqDataSetFromMatrix(countData = data,  
                              colData = meta,
                              design = ~ genotyping)  # Comparing factor 


# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE) 

# exploratory data analysis (PCA & heirarchical clustering) - identifying outliers and sources of variation in the data
## PCA ##
p1 <- plotPCA(rld, intgroup="genotyping", ntop=500, returnData = F) +
  ggtitle ("PCA") + # main title  
  theme_light() +  # Theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size = 1, colour = "black"),
        legend.title=element_blank()
  )+
  #geom_text_repel(aes(label = colnames(data))) +
  #coord_cartesian(ylim = c(-30,40), xlim = c(-60,40)) + # Coordinate system
  scale_color_manual(values = c("#1f78b4","#4daf4a","#e41a1c","gray","#e5d8bd","#a6cee3"))+
  coord_fixed()# Color 
p1

pdf("figures/pca.pdf")
p1
dev.off()


## ntop filtering function

ntop <- function(object, ntop = 500){
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  assay(object)[select, ]
}

rld_ntop <- ntop(rld) # Extract the rlog matrix from the object (top 500 var by default)
ntop_gene <- rownames(rld_ntop)


pca  <-  prcomp(t(rld_ntop), center=TRUE, scale=TRUE)



### extract PC1 and PC2 genes ###
PC1_top100 <- rownames(data.frame(sort(abs(pca$rotation[,"PC1"]), decreasing=TRUE)[1:100]))
PC2_top100 <- rownames(data.frame(sort(abs(pca$rotation[,"PC2"]), decreasing=TRUE)[1:100]))
write.table(PC1_top100, file="results/DE_table/PC1_top100_vivo_ileum_colon_satb2ko.txt", sep="\t", quote=F, col.names=NA)
write.table(PC2_top100, file="results/DE_table/PC2_top100.txt", sep="\t", quote=F, col.names=NA)

##################### Data Visualization ######################
# 1. extract vst assay 
# Transform counts for data visualization
vsd <- vst(dds, blind=F)

plotPCA(vsd, intgroup="genotyping")

vsd.df <- vsd %>% assay() %>% data.frame()

# 2. Heatmap
# write a function to plot heatmap using the same arguments
plotHM <- function(data,title=""){
  pheatmap(data.frame(data), 
           main = title,
           #color = colorRampPalette(c("blue","white","red"))(50),
           #color = magma(50),
           color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(50),
           #breaks = seq(-2,2,by=0.1),
           scale="row", 
           cluster_rows = T, cluster_cols = F,
           show_colnames = T, show_rownames = T,
           border_color = NA, 
           fontsize = 10, 
           fontsize_row = 10, 
           height=20)
}

### DEG quantification and visualization

## plot
p1 <- vsd.df[PC1_top100,] %>% plotHM("pc1_100")
p2 <- vsd.df[PC2_top100,] %>% plotHM("pc2_100")
plot_grid(p1$gtable,p2$gtable)
# Volcano plot

## Volcano plot function
plotVol <- function(table_tb,title="") {
  ggplot(table_tb) +
    
    geom_point(aes(x = log2FoldChange, 
                   y = -log10(padj), 
                   colour = threshold),
               size=1) +
    ggtitle(title) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    scale_y_continuous(limits=c(0,40)) +
    scale_x_continuous(limits=c(-5,5)) +
    theme_light() + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linetype = "solid", 
                                  fill = NA, 
                                  size = 1, 
                                  colour = "black"),
      legend.position = "none",
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.title = element_text(size = rel(1.25))
    )+
    scale_color_manual(values = c("black", "#F75554")) # Color 
}

## Obtain logical vector regarding whether padj values are less than 0.05 and fold change greater than 1.5 in either direction

resLFC.tb <- resLFC.tb %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 1)

# plot Volcano 
plotVol(resLFC.tb)

## average gene expression of each group
normalized_counts %>% head()
sample.list <- list(ileum=rownames(meta)[(1:3)],colon=rownames(meta)[(4:6)], satb2.ko=rownames(meta)[(7:9)])
output.list <- vector("list", 3)
names(output.list) <- c("ileum","colon","satb2.ko")
for (i in seq(sample.list)) {
  x <- normalized_counts[,sample.list[[i]]]
  output.list[[i]] <- rowMeans(x) %>% as.data.frame() %>% rownames_to_column("gene") %>% as_tibble()
  colnames(output.list[[i]])[2] <- names(output.list)[i]
}
output.list[[1]] %>% inner_join(output.list[[2]]) %>% inner_join(output.list[[3]]) %>% write_csv("results/mean_expression/mean.csv")

# Output the versions of all tools used in the DE analysis:
sessionInfo()
