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
#library(DEGreport)
#library(BiocParallel)
#library(apeglm)
library(ggfortify)
library(data.table)
library(limma)
#library(GenomicFeatures)
library(dplyr)
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

# load feature counts
data.raw <- read.table("data/mta2_cleancounts.txt", header=T, row.names=1, check.names = F)  # Put your unnormalized counts file in data directory
# extract protein_coding gene counts from the raw counts
coding_genes <- grcm38 %>% dplyr::filter(biotype == "protein_coding")
coding_genes <- coding_genes$symbol
idx <- row.names(data.raw) %in% c(coding_genes,"Inava")
data <- data.raw[idx,]
write.table(data, file="data/data.txt", sep="\t", quote=F, col.names=NA)
# load data
data <- read.table("data/data.txt", row.names=1, header=T)
# load and subset meta
meta <- read.csv("meta/meta.csv", row.names=1, header=T) # Put your meta file in meta directory
idx <- row.names(meta) %in% colnames(data)
#meta <- meta[idx,]
#meta$genotyping_region <- paste(meta$genotyping, meta$region, sep = "_")

#meta$batch <- meta$batch %>% as_factor()

# 3. Check the input files
##  (1) Check the data class to confirm they are "data.frame"
class(data)
class(meta)
##  (2) Check the colume names of data match the row names of meta: both should be TRUE
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

# Check if RNA distribution and statistical properties fit NB model (mean>variance)
# RNA distribution

#### ------ below is the function ------ ####
plotHistFunc <- function(dataframe, pre, na.rm = TRUE, ...) {     # make a function
  nm <- colnames(dataframe)              # put column names in the variable nm
  for (i in seq_along(nm)) {            # for loop, extract the index number from nm
    plots <- ggplot(dataframe,aes_string(x = nm[i])) +    # ggplot(file,aes_from_string)
      geom_histogram(stat ="bin", bins = 200, fill = "#446179", alpha = 0.8) + # plot histogram
      xlim(-5, 300)  + # x axis limit
      xlab("Raw expression counts") + # x lable
      ylab("Number of genes") + # y lable
      theme_light() + # theme
      theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 1, colour = "black")) # border
    ggsave(plots,filename=paste("results/count_distribution/",pre,nm[i],".png",sep="")) # export
  }}
#### ------ above is the function ------ ####
#plotHistFunc(data,"whole_") # plot counts distribution for RNA samples automatically. Call the function to check other dataset if you want.

##  (2) mean VS. variance # if var > mean, then model fits
mean_counts <- apply(data, MARGIN = 1, mean) # calculate mean for all genes
variance_counts <- apply(data, MARGIN = 1, var) # calculate variance for all genes
df <- data.frame(mean_counts, variance_counts) # make a data frame: x=mean, y=variance
meanVSvar <- ggplot(df) + # plot the dataframe
  geom_point(aes(x=mean_counts, y=variance_counts),alpha = 0.4) + #plot points with 60% transparency
  geom_line(aes(x=mean_counts, y=mean_counts), color = "red") +  # red line indicates mean=variance
  scale_y_log10() + scale_x_log10() +
  theme_light() +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 1, colour = "black")) # border
#ggsave(meanVSvar,filename=paste("results/meanVSvar/meanVSvar.png",sep="")) # export 

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
#PC1_top100 <- rownames(data.frame(sort(abs(pca$rotation[,"PC1"]), decreasing=TRUE)[1:100]))
#PC2_top100 <- rownames(data.frame(sort(abs(pca$rotation[,"PC2"]), decreasing=TRUE)[1:100]))
#write.table(PC1_top100, file="results/DE_table/PC1_top100.txt", sep="\t", quote=F, col.names=NA)
#write.table(PC2_top100, file="results/DE_table/PC2_top100.txt", sep="\t", quote=F, col.names=NA)

## plot heatmap ##
# Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
# Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames
# plot heatmap
p1 <- pheatmap(rld_cor) # Check the heatmap plot # Zoom and copy the vector figure
pdf("figures/heatmap_cor.pdf")
p1
dev.off()
######################### DE Analysis #########################
# **Optional step** - Re-create DESeq2 dataset if the design formula has changed after QC analysis in include other sources of variation
dds <- DESeqDataSetFromMatrix(countData = data,  
                              colData = meta,
                              design = ~ genotyping)  # Comparing factor 

# Extract normalized counts, check the size factors and total counts before and after normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds) # check the size factors
normalized_counts <- counts(dds, normalized=TRUE)

write.table(normalized_counts, file="results/normalization/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

## Take a look at the total number of reads for each sample using
colSums(counts(dds)) # Sum the unnormalized column values
colsum <- colSums(counts(dds)) 
write.csv(colsum, file="results/colSum.csv") # Export

colSums(normalized_counts) # Sum the normalized column values
ncolsum <- colSums(normalized_counts) 
write.csv(ncolsum, file = "results/nor_colSum.csv")

# change reference level
dds$genotyping <- dds$genotyping %>% relevel(ref = "Ctrl")
# Differential expression analysis
dds <- DESeq(dds)
sizeFactors(dds)

# Examine dispersion estimates plot to ensure your data is a good fit for the DESeq2 model
plotDispEsts(dds) 

# Create contrasts to perform Wald testing on the shrunken log2 foldchanges between specific conditions:
# check the coef for shurnkage
resultsNames(dds)
saveRDS(dds, "rds/dds.rds")
dds <- readRDS("rds/dds.rds")
## Convert the result table into a tibble
## write a function
res2tib <- function(result_table) {
  result_table %>% 
    data.frame() %>% 
    rownames_to_column(var="gene") %>% #convert rownames to a column
    as_tibble() # covert the dataframe to tibble
}  ## res2tib USAGE: res2tib(result_table)

# Extract DE_table without or with shrunken
name.list <- resultsNames(dds)[-1] %>% as.list()
names(name.list) <- resultsNames(dds)[-1] 
name.list
res.list <- name.list
resLFC.list <- name.list
resLFC.tb.list <- name.list
resLFC.normal.list <- name.list
resLFC.normal.tb.list <- name.list

## use coefficient
for (i in seq(res.list)) {
  res.list[[i]] <- results(dds, name = name.list[[i]][1], alpha = 0.1)  # defaut alpha = 0.1
  
#  resLFC.list[[i]] <- lfcShrink(dds, coef = name.list[[i]][1], res=res.list[[i]])
#  resLFC.tb.list[[i]] <- resLFC.list[[i]] %>% res2tib()
#  resLFC.tb.list[[i]] %>% write_csv(paste0("results/DE_table/resLFC_apeglm_",name.list[[i]][1],".csv"))
                                    
  resLFC.normal.list[[i]] <- lfcShrink(dds, coef = name.list[[i]][1], res=res.list[[i]], type = 'normal')
  resLFC.normal.tb.list[[i]] <- resLFC.normal.list[[i]] %>% res2tib()
  resLFC.normal.tb.list[[i]] %>% write_csv(paste0("results/DE_table/resLFC_normal_",name.list[[i]][1],".csv"))
}

#resLFC.tb.list %>% saveRDS("rds/resLFC.tb.list.rds")
#resLFC.tb.list <- readRDS("rds/resLFC.tb.list.rds")
resLFC.normal.tb.list %>% saveRDS("rds/resLFC.normal.tb.list.rds")

#-------------------------------------------------
#Run stop here and go to functional analysis for GSEA


## use contrast
#contrast.list <- list(cdx2ko=c("genotyping", "CHD4_KO_Crypt", "CHD4_CTRL_Crypt"),
#                      hnf4ako=c("genotyping", "hnf4ako", "wt_colon"),
#                      satb2_hnf4a=c("genotyping", "satb2_hnf4a_ko", "wt_colon"),
#                      satb2ko=c("genotyping", "satb2ko", "wt_colon"),
#                      wt_ileum=c("genotyping", "wt_ileum", "wt_colon"))

#res.contrast.list <- contrast.list
#resLFC.normal.list <- contrast.list
#resLFC.normal.tb.list <- contrast.list
#for (i in seq(contrast.list)) {
#  res.contrast.list[[i]] <- results(dds, contrast=contrast.list[[i]], alpha = 0.1)  # defaut alpha = 0.1
#  resLFC.normal.list[[i]] <- lfcShrink(dds, contrast=contrast.list[[i]], res=res.contrast.list[[i]], type = 'normal')
#  resLFC.normal.tb.list[[i]] <- resLFC.normal.list[[i]] %>% res2tib()
#  resLFC.normal.tb.list[[i]] %>% write.csv(paste0("results/DE_table/resLFC_",names(contrast.list)[i],"normal.csv"))
#}

## MA plot
#DESeq2::plotMA(res.list$genotyping_satb2ko_vs_wt_colon, alpha = 0.01, main = "unshrunken")
#DESeq2::plotMA(resLFC.list$genotyping_satb2ko_vs_wt_colon, alpha = 0.01, main = "apeglm")

## Summarize results
summary(resLFC.list$genotyping_satb2ko_vs_wt_colon)

# Subset the table to only keep the significant genes using our pre-defined thresholds
## write a function
sigDE <- function(res_table_tb){
  res_table_tb %>% 
    dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) # filter the genes by thresholds
}

# Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 1
#base.cutoff <- 10

# DEG research table
#DEG.apeglm.list <- resLFC.tb.list
DEG.normal.list <- resLFC.normal.tb.list

for (i in seq(DEG.normal.list)) {
#  DEG.apeglm.list[[i]] <- resLFC.tb.list[[i]] %>% sigDE()
  DEG.normal.list[[i]] <- resLFC.normal.tb.list[[i]] %>% sigDE()
#  DEG.apeglm.list[[i]] %>% write.csv(paste0("results/DE_table/DEG_apeglm_",name.list[[i]][1],".csv"))
  DEG.normal.list[[i]] %>% write.csv(paste0("results/DE_table/DEG_normal_",name.list[[i]][1],".csv"))
}


############################################################################################################################################
# DEG quantification

## write a function to extract important numbers with different cutoff parameters.
## No rpkm filtering ##
getEnrichedNumber <- function(table_1, padj.cutoff_1, lfc.cutoff_1, table_2, padj.cutoff_2, lfc.cutoff_2){
  
  DEG_1 <- table_1 %>% dplyr::filter(padj < padj.cutoff_1 & abs(log2FoldChange) > lfc.cutoff_1)
  DEG_2 <- table_2 %>% dplyr::filter(padj < padj.cutoff_2 & abs(log2FoldChange) > lfc.cutoff_2)
  
  DEG_1_up <- DEG_1 %>% dplyr::filter(log2FoldChange>0)
  DEG_1_down <- DEG_1 %>% dplyr::filter(log2FoldChange<0)
  DEG_2_up <- DEG_2 %>% dplyr::filter(log2FoldChange>0)
  DEG_2_down <- DEG_2 %>% dplyr::filter(log2FoldChange<0)
  
  n_up_overlap <- (DEG_1_up$gene %in% DEG_2_up$gene) %>% which() %>% length()
  n_down_overlap <- (DEG_1_down$gene %in% DEG_2_down$gene) %>% which() %>% length()
  
  cat(paste0("cutoff of table_1 (1) padj ", padj.cutoff_1, " (2) lfc ", lfc.cutoff_1, ".\n"))
  cat(paste0("cutoff of table_2 (1) padj ", padj.cutoff_2, " (2) lfc ", lfc.cutoff_2, ".\n"))
  
  cat(paste0("The number of genes up-regulated in table_1: ", DEG_1_up %>% nrow(), "\n"))
  cat(paste0("The number of genes up-regulated in table_2: ", DEG_2_up %>% nrow(), "\n"))
  cat(paste0("The number of genes up-regulated in both tables: ", n_up_overlap, "\n"))
  cat(paste0("overlap/(up-regulated gene in table_1): ", n_up_overlap/(DEG_1_up %>% nrow()), "\n"))
  cat(paste0("overlap/(up-regulated gene in table_2): ", n_up_overlap/(DEG_2_up %>% nrow()), "\n"))
  cat(paste0("The number of genes down-regulated in table_1: ", DEG_1_down %>% nrow(), "\n"))
  cat(paste0("The number of genes down-regulated in table_2: ", DEG_2_down %>% nrow(), "\n"))
  cat(paste0("The number of genes down-regulated in both tables: ", n_down_overlap, "\n"))
  cat(paste0("overlap/(down-regulated gene in table_1): ", n_down_overlap/(DEG_1_down %>% nrow()), "\n"))
  cat(paste0("overlap/(down-regulated gene in table_2): ", n_down_overlap/(DEG_2_down %>% nrow()), "\n"))
  cat("\n")
  
}


#### USAGE: getEnrichedNumber <- function(table_1, padj.cutoff_1, lfc.cutoff_1, table_2, padj.cutoff_2, lfc.cutoff_2)
#### Input: 
#### padj.cutoff_1, lfc.cutoff_1 for research table 1
#### padj.cutoff_2, lfc.cutoff_2 for research table 2
#### Output:
#### Basic statistics of gene number overlapped gene number.



############################################################################################################################################
  

## write a function to extract important numbers with different cutoff parameters.
## With rpkm filtering ##

getEnrichedNumber2 <- function(a_b_tb1, padj.cutoff_ab, lfc.cutoff_ab, c_d_tb2, padj.cutoff_cd, lfc.cutoff_cd, a, b, c, d, rpkm.df, rpkm_filter,sample.filter){
  
  sig_a_vs_b <- a_b_tb1 %>% dplyr::filter(padj < padj.cutoff_ab & abs(log2FoldChange) > lfc.cutoff_ab)
  sig_c_vs_d <- c_d_tb2 %>% dplyr::filter(padj < padj.cutoff_cd & abs(log2FoldChange) > lfc.cutoff_cd)
  
  ## filter
  keep <- rowSums(rpkm.df > rpkm_filter) >= sample.filter 
  rpkm_filtered <- rpkm.df[keep,]
  
  ## Create tibbles including row names
  rpkm_filtered.tb <- rpkm_filtered %>% 
    rownames_to_column(var="gene") %>% 
    as_tibble() # meta data in tibble format
  
  aEnriched_vs_b <- dplyr::filter(sig_a_vs_b,log2FoldChange>0)$gene %in% rpkm_filtered.tb$gene %>% which() %>% length()
  bEnriched_vs_a <- dplyr::filter(sig_a_vs_b,log2FoldChange<0)$gene %in% rpkm_filtered.tb$gene %>% which() %>% length()
  cEnriched_vs_d <- dplyr::filter(sig_c_vs_d,log2FoldChange>0)$gene %in% rpkm_filtered.tb$gene %>% which() %>% length()
  dEnriched_vs_c <- dplyr::filter(sig_c_vs_d,log2FoldChange<0)$gene %in% rpkm_filtered.tb$gene %>% which() %>% length()
  
  aEnriched_overlap_cEnriched <- dplyr::filter(sig_a_vs_b,log2FoldChange>0 & gene %in% rpkm_filtered.tb$gene)$gene %in% 
    dplyr::filter(sig_c_vs_d, log2FoldChange>0 & gene %in% rpkm_filtered.tb$gene)$gene %>% 
    which() %>% length()
  bEnriched_overlap_dEnriched <- dplyr::filter(sig_a_vs_b,log2FoldChange<0 & gene %in% rpkm_filtered.tb$gene)$gene %in% 
    dplyr::filter(sig_c_vs_d, log2FoldChange<0 & gene %in% rpkm_filtered.tb$gene)$gene %>% 
    which() %>% length()
  cat(paste0("cutoff of ", a, " and ", b,  ": padj ",padj.cutoff_ab," & lfc ", lfc.cutoff_ab,".\n"))
  cat(paste0("cutoff of ", c, " and ", d,  ": padj ",padj.cutoff_cd," & lfc ", lfc.cutoff_cd,".\n"))
  cat(paste0("RPKM cutoff: ", rpkm_filter,".\n"))
  cat(paste0("Sample number cutoff: ", sample.filter, ".\n"))
  cat(paste0(a,">",b," - number: ",aEnriched_vs_b,"\n"))
  cat(paste0(a,"<",b," - number: ",bEnriched_vs_a,"\n"))
  cat(paste0(a,">",b," overlapped with ",c,">",d," - number: ",aEnriched_overlap_cEnriched,"; ratio: ", round(aEnriched_overlap_cEnriched/aEnriched_vs_b*100, 2),"%\n"))
  cat(paste0(a,"<",b," overlapped with ",c,"<",d," - number: ",bEnriched_overlap_dEnriched,"; ratio: ", round(bEnriched_overlap_dEnriched/bEnriched_vs_a*100, 2),"%\n"))
  cat("\n")
  
  cat(paste0("Export DEG table list filtered by rpkm"))
  cat("\n")
  idx <- sig_a_vs_b$gene %in% rpkm_filtered.tb$gene
  a_vs_b <- sig_a_vs_b[idx,]
  idx <- sig_c_vs_d$gene %in% rpkm_filtered.tb$gene
  c_vs_d <- sig_c_vs_d[idx,]
  list <- list(a_vs_b=a_vs_b, c_vs_d=c_vs_d)
  return(list)
  
}

#### USAGE: getEnrichedNumber2 <- function(a_b_tb1, padj.cutoff_ab, lfc.cutoff_ab, c_d_tb2, padj.cutoff_cd, lfc.cutoff_cd, a, b, c, d, rpkm.df, rpkm_filter,sample.filter)

#### Input: 
#### padj.cufoff_ab for a vs b
#### lfc.cutoff_ab for a vs b
#### padj.cutoff_cd for c vs d
#### lfc.cutoff_cd for c vs d
#### a_b_tb1: result table tibble for a vs b
#### c_d_tb2: result table tibble for c vs d
#### a,b,c,d "title" for a,b,c and d
#### rpkm.df: rpkm data frame
#### rpkm_filter: filter value of rpkm
############################################################################################################################################

##################### Data Visualization ######################
# 1. extract vst assay 
# Transform counts for data visualization
vsd <- vst(dds, blind=F)
vsd.batchRM <- vsd
assay(vsd.batchRM) <- dds %>% vst(blind=T) %>% assay() %>% limma::removeBatchEffect(vsd$library)

plotPCA(vsd, intgroup="genotyping")
plotPCA(vsd.batchRM, intgroup="genotyping")

vsd.df <- vsd %>% assay() %>% data.frame()

# 1.1 filter with RPKM
rpkm.df <- read.table("results/rpkm.txt", sep="\t", row.names = 1, header = T)
## set up filter threshold
rpkm.filter <- 8 # rpkm threshold
sample.filter <- 2 # sample number threshold
## filter
keep <- rowSums(rpkm.df > rpkm.filter) >= sample.filter 
rpkm.df.filtered <- rpkm.df[keep,]
idx <- rownames(vsd.df) %in% rownames(rpkm.df.filtered)
vsd.df.rpkm.filtered <- vsd.df[idx,]
nrow(vsd.df.rpkm.filtered)
nrow(vsd.df)
# 2. Heatmap (Filter with RPKM counts)
# write a function to plot heatmap using the same arguments

  plotHM <- function(data,title=""){
  pheatmap(data.frame(data), 
           main = title,
           #color = colorRampPalette(c("cornflowerblue","beige","brown1"))(50),
           color = magma(50),
           #color = colorRampPalette(rev(brew er.pal(n = 11, name ="RdBu")))(50),
           #breaks = seq(-1.5,2.5,by=0.1),
           scale="row", 
           cluster_rows = T, cluster_cols = F,
           show_colnames = T, show_rownames = T,
           border_color = NA, 
           fontsize = 10, 
           fontsize_row = 10, 
           height=20)
}

  
### DEG quantification and visualization
## plot with fixed gene list
Absorption.df <- read.table("results/absorption.txt", header = T)
feature.to.plot <- Absorption.df$gene
sample.to.plot <- c("GS4","GS5","GS6",
                    "GS3","GS7")


#feature.to.plot <- filtered.DE.list$c_vs_d$gene
## subset a dataframe based on a given gene list
vsd.df.rpkm.filtered %>% head()
## which row exist in the given gene list?
idx <- rownames(vsd.df.rpkm.filtered) %in% feature.to.plot
vsd.df.rpkm.filtered.sub <- vsd.df.rpkm.filtered[idx, sample.to.plot]
## plot a heatmap
p1 <- vsd.df.rpkm.filtered.sub %>% plotHM()
pdf("figures/heatmap_absorption_no_cluster.pdf")
p1
dev.off()

microvillus_organization.df <- read.table("results/microvillus_organization.txt", header = T)
feature.to.plot <- microvillus_organization.df$gene
sample.to.plot <- c("GS4","GS5","GS6",
                    "GS3","GS7")

length(feature.to.plot)
vsd.df[feature.to.plot,sample.to.plot] %>% nrow()
#feature.to.plot <- filtered.DE.list$c_vs_d$gene
## subset a dataframe based on a given gene list
vsd.df %>% head()
## which row exist in the given gene list?
idx <- rownames(vsd.df) %in% feature.to.plot
vsd.df.sub <- vsd.df[idx, sample.to.plot]
## plot a heatmap
p1 <- vsd.df.sub %>% plotHM()
pdf("figures/heatmap_microvillus_organization.pdf")
p1
dev.off()

