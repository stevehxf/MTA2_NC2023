# This script is designed for differential expression (DE) analysis using DEseq2 as well as data visualization
# apeglm shrunkage can be used in this script

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
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(DEGreport)
library(apeglm)
library(ggfortify)
library(data.table)
library(GenomicFeatures)
library(annotables)
library(edgeR)
library(ggpubr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DOSE)
library(pathview)
library(purrr)
library(clusterProfiler)
library(enrichplot)
library(Rgraphviz)

# load feature counts
data.raw <- read.table("data/3D_cleancounts_4.txt", header=T, row.names=1, check.names = F)  # Put your unnormalized counts file in data directory

# extract protein_coding gene counts from the raw counts
coding_genes <- grcm38 %>% filter(biotype == "protein_coding")
coding_genes <- coding_genes$symbol

idx <- row.names(data.raw) %in% coding_genes
data <- data.raw[idx,]
write.table(data, file="data/data.txt", sep="\t", quote=F, col.names=NA)
# load data
data <- read.table("data/data.txt", row.names=1, header=T)
# load meta
meta <- read.table("meta/meta_all_4.txt", row.names=1, header=T) # Put your meta file in meta directory
data <- data[,c("A1","A2","A3",sprintf("A%d", 10:49))]
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


# Usage of the function: plotHistFunc(dataframe,"prefix_of_output_file_name")
plotHistFunc(data,"whole_") # plot counts distribution for RNA samples automatically. Call the function to check other dataset if you want.


##  (2) mean VS. variance # if var > mean, then model fits
mean_counts <- apply(data, MARGIN = 1, mean) # calculate mean for all genes
variance_counts <- apply(data, MARGIN = 1, var) # calculate variance for all genes
df <- data.frame(mean_counts, variance_counts) # make a data frame: x=mean, y=variance

meanVSvar <- ggplot(df) + # plot the dataframe
  geom_point(aes(x=mean_counts, y=variance_counts),alpha = 0.4) + #plot points with 60% transparency
  geom_line(aes(x=mean_counts, y=mean_counts), color = "red") +  # red line indicates mean=variance
  scale_y_log10() +
  scale_x_log10() +
  theme_light() + # theme
  theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 1, colour = "black")) # border

ggsave(meanVSvar,filename=paste("results/meanVSvar/meanVSvar.png",sep="")) # export 

# Genes that are barely detectable will contribute noise to the resulting clusters
# pre-filtering: only keep the genes that at least 3 (n.rep.min) samples have more than 5 counts
n.min <- 3
keep <- rowSums((data)>10) >= n.min
data <- data[keep,]

# create DESeq2Dataset object. Count normalization for QC, etc.
dds <- DESeqDataSetFromMatrix(countData = data,  
                              colData = meta,
                              design = ~ sgRNA)  # Comparing factor 


# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE) 

# exploratory data analysis (PCA & heirarchical clustering) - identifying outliers and sources of variation in the data
## PCA ##
plotPCA(rld, intgroup="sgRNA", ntop=500, returnData = F) +
  ggtitle ("PCA") + # main title  
  theme_light() +  # Theme
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size = 1, colour = "black"),
        legend.title=element_blank()
        )+
  geom_text_repel(aes(label = colnames(data)))
# coord_cartesian(ylim = c(-30,30), xlim = c(-30,70)) + # Coordinate system
# scale_color_manual(values = c("#244D6F","#793F36","#7CAABC","#FF968D","#929292","black")) +   # Color 
  
# PCA score for each gene
#pca  <-  prcomp(t(rld_ntop), center=TRUE, scale=TRUE)

### extract PC1 and PC2 genes ###
#PC1_top100 <- rownames(data.frame(sort(abs(pca$rotation[,"PC1"]), decreasing=TRUE)[1:100]))
#PC2_top100 <- rownames(data.frame(sort(abs(pca$rotation[,"PC2"]), decreasing=TRUE)[1:100]))
#write.table(PC1_top100, file="results/DE_table/PC1_top100.txt", sep="\t", quote=F, col.names=NA)
#write.table(PC2_top100, file="results/DE_table/PC2_top100.txt", sep="\t", quote=F, col.names=NA)

## ntop filtering function
#ntop <- function(object, ntop = 500){
  #rv <- rowVars(assay(object))
  #select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  #assay(object)[select, ]
#}
#rld_ntop <- ntop(rld) # Extract the rlog matrix from the object (top 500 var by default)

## plot heatmap ##
# Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
# Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
#rld_ntop_cor <- cor(rld_ntop) ## top genes correlation
head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames
# plot heatmap
pheatmap(rld_cor) # Check the heatmap plot # Zoom and copy the vector figure
#pheatmap(rld_ntop_cor) # top genes heatmap

# rpkm normalization ##############################################################################################################
data.raw <- read.table("data/data.txt", header=T, row.names=1, check.names = F)  # Put your unnormalized counts file in data directory
data.raw <- data.raw[,c("A1","A2","A3",sprintf("A%d", 10:49))]
###  Input GTF file to calculate gene length  ###
# First, import the GTF-file that you have also used as input for htseq-count
#txdb <- makeTxDbFromGFF("./data/gencode.v19.annotation.gtf",format="gtf")
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
write.table(rpkm.df, "results/normalization/rpkm.txt", sep="\t", quote=F, col.names=NA)
rpkm.df <- read.table("results/normalization/rpkm.txt", sep="\t", row.names = 1, header = T)
rpkm.tb <- rpkm.df %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble() # meta data in tibble format
#############################################################################################################################################

######################### DE Analysis #########################
# **Optional step** - Re-create DESeq2 dataset if the design formula has changed after QC analysis in include other sources of variation
dds <- DESeqDataSetFromMatrix(countData = data,  
                              colData = meta,
                              design = ~ sgRNA)  # Comparing factor 

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


# Differential expression analysis
dds <- DESeq(dds)
sizeFactors(dds)

# Examine dispersion estimates plot to ensure your data is a good fit for the DESeq2 model
plotDispEsts(dds) 
saveRDS(dds, "data/dds.rds")

#1 Satb2 gene DEG analysis
contrast_satb2 <- c("sgRNA", "satb2", "cas9")
resultsNames(dds)
res_satb2 <- results(dds, contrast=contrast_satb2, alpha = 0.1)  # defaut alpha = 0.1
resLFC_satb2 <- lfcShrink(dds, contrast=contrast_satb2, res=res_satb2, type="normal")
resLFC.ape_satb2 <- lfcShrink(dds, coef = "sgRNA_satb2_vs_cas9", type="apeglm") # shrink by apeglm
write.csv(resLFC_satb2, file = "results/DE_table/resLFC_satb2.csv")
write.csv(resLFC.ape_satb2, file = "results/DE_table/resLFC.ape_satb2.csv")
DESeq2::plotMA(res_satb2, ylim=c(-5,5), alpha = 0.05, main = "unshrunken")
DESeq2::plotMA(resLFC_satb2, ylim=c(-5,5),  alpha = 0.05, main = "normal")
DESeq2::plotMA(resLFC.ape_satb2, ylim=c(-5,5),  alpha = 0.05, main = "apeglm")

#1 ileum and colon gene DEG analysis
contrast_ileum<- c("sgRNA", "ileum", "cas9")
resultsNames(dds)
res_ileum <- results(dds, contrast=contrast_ileum, alpha = 0.1)  # defaut alpha = 0.1
resLFC_ileum <- lfcShrink(dds, contrast=contrast_ileum, res=res_ileum, type="normal")
resLFC.ape_ileum <- lfcShrink(dds, coef = "sgRNA_ileum_vs_cas9", type="apeglm") # shrink by apeglm
write.csv(resLFC_ileum, file = "results/DE_table/resLFC_ileum.csv")
write.csv(resLFC.ape_ileum, file = "results/DE_table/resLFC.ape_ileum.csv")
DESeq2::plotMA(res_ileum, ylim=c(-5,5), alpha = 0.05, main = "unshrunken")
DESeq2::plotMA(resLFC_ileum, ylim=c(-5,5),  alpha = 0.05, main = "normal")
DESeq2::plotMA(resLFC.ape_ileum, ylim=c(-5,5),  alpha = 0.05, main = "apeglm")

#3Smarcd2 gene DEG analysis
contrast_smarcd2 <- c("sgRNA", "smarcd2", "cas9")
resultsNames(dds)
res_smarcd2 <- results(dds, contrast=contrast_smarcd2, alpha = 0.1)  # defaut alpha = 0.1
resLFC_smarcd2 <- lfcShrink(dds, contrast=contrast_smarcd2, res=res_smarcd2, type="normal")
resLFC.ape_smarcd2 <- lfcShrink(dds, coef = "sgRNA_smarcd2_vs_cas9", type="apeglm") # shrink by apeglm
write.csv(resLFC_smarcd2, file = "results/DE_table/resLFC_smarcd2.csv")
write.csv(resLFC.ape_smarcd2, file = "results/DE_table/resLFC.ape_smarcd2.csv")
DESeq2::plotMA(res_smarcd2, ylim=c(-5,5), alpha = 0.05, main = "unshrunken")
DESeq2::plotMA(resLFC_smarcd2, ylim=c(-5,5),  alpha = 0.05, main = "normal")
DESeq2::plotMA(resLFC.ape_smarcd2, ylim=c(-5,5),  alpha = 0.05, main = "apeglm")
#3Smarcd2 gene DEG analysis
contrast_smarcd2 <- c("sgRNA", "smarcd2", "cas9")
resultsNames(dds)
res_smarcd2 <- results(dds, contrast=contrast_smarcd2, alpha = 0.1)  # defaut alpha = 0.1
resLFC_smarcd2 <- lfcShrink(dds, contrast=contrast_smarcd2, res=res_smarcd2, type="normal")
resLFC.ape_smarcd2 <- lfcShrink(dds, coef = "sgRNA_smarcd2_vs_cas9", type="apeglm") # shrink by apeglm
write.csv(resLFC_smarcd2, file = "results/DE_table/resLFC_smarcd2.csv")
write.csv(resLFC.ape_smarcd2, file = "results/DE_table/resLFC.ape_smarcd2.csv")
DESeq2::plotMA(res_smarcd2, ylim=c(-5,5), alpha = 0.05, main = "unshrunken")
DESeq2::plotMA(resLFC_smarcd2, ylim=c(-5,5),  alpha = 0.05, main = "normal")
DESeq2::plotMA(resLFC.ape_smarcd2, ylim=c(-5,5),  alpha = 0.05, main = "apeglm")

# Mta2 gene DEG analysis
# Create contrasts to perform Wald testing on the shrunken log2 foldchanges between specific conditions:
contrast_mta2 <- c("sgRNA", "mta2", "cas9")

# check the coef for shurnkage
resultsNames(dds)

# Extract DE_table without or with shrunken
res_mta2 <- results(dds, contrast=contrast_mta2, alpha = 0.1)  # defaut alpha = 0.1
resLFC_mta2 <- lfcShrink(dds, contrast=contrast_mta2, res=res_mta2, type="normal")

resLFC.ape_mta2 <- lfcShrink(dds, coef = "sgRNA_mta2_vs_cas9", type="apeglm") # shrink by apeglm

write.csv(resLFC_mta2, file = "results/DE_table/resLFC_mta2.csv")
write.csv(resLFC.ape_mta2, file = "results/DE_table/resLFC.ape_mta2.csv")

## MA plot
DESeq2::plotMA(res_mta2, ylim=c(-5,5), alpha = 0.05, main = "unshrunken")
DESeq2::plotMA(resLFC_mta2, ylim=c(-5,5),  alpha = 0.05, main = "normal")
DESeq2::plotMA(resLFC.ape_mta2, ylim=c(-5,5),  alpha = 0.05, main = "apeglm")

## Summarize results
#summary(resLFC_smarcd2)

## Convert the result table into a tibble
## write a function
df2tib <- function(result_table) {
  result_table %>% 
    data.frame() %>% 
    rownames_to_column(var="gene") %>% #convert rownames to a column
    as_tibble() # covert the dataframe to tibble
}  ## res2tib USAGE: res2tib(result_table)

resLFC.ape_ileum.tb <- resLFC.ape_ileum %>% df2tib()
##setting cutoff of resLFC.gene

resLFC.ape_ileum_base.tb <- resLFC.ape_ileum.tb %>% filter(baseMean >= base.cutoff)

#write.csv(resLFC.tb, file = "results/DE_table/resLFC.tb.csv")

# Subset the table to only keep the significant genes using our pre-defined thresholds
## write a function
sigDE <- function(res_table_tb){
  res_table_tb %>% 
    filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) # filter the genes by thresholds
}


# Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.5
base.cutoff <- 100
# DEG research table
ileum.ape.DEG <- sigDE(resLFC.ape_ileum.tb)
ileum.ape.DEG.baseFiltered <- ileum.ape.DEG %>% filter(baseMean >= base.cutoff)

write.csv(ileum.ape.DEG, file = "results/DE_table/smarcd2.ape.DEG.csv")
write.csv(ileum.ape.DEG.baseFiltered, file = "results/DE_table/smarcd2.ape.DEG.baseFiltered.csv")

##################### Data Visualization ######################

# 1. Create tibbles including row names
meta_tb <- meta %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble() # meta data in tibble format

normalized_counts_tb <- counts(dds, normalized=TRUE) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble() # normalized counts as tibble

# 2. Heatmap

## Write a function to extract normalized counts for significant DE genes and set the gene column to row names
SigNorCounts <- function(nor_counts=normalized_counts_tb, sig) {
  nor_counts %>% 
    filter(gene %in% sig$gene) %>% 
    data.frame() %>%
    column_to_rownames(var = "gene") 
}## SigNorCounts USAGE: SigNorCounts(normalized_counts, significant DE data)

## Extract using SigNorCounts function
norm.ileum.ape.DEG.baseFiltered <- SigNorCounts(normalized_counts_tb,ileum.ape.DEG.baseFiltered)

# display.brewer.all()
## Run pheatmap
# write a function to plot heatmap using the same arguments
plotHM <- function(data,title=""){
  pheatmap(data.frame(data), 
           main = title,
           #color = colorRampPalette(c("black","red"))(12),
           #color = inferno(50),
           color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(50),
           #breaks = seq(-1,1,by=0.1),
           #color = colorRampPalette(c("green", "black","red"))(40),
           
           cluster_rows = T, 
           cluster_cols = F,
           show_colnames = T,
           show_rownames = F,
           # annotation = annotation, 
           annotation_names_col = F,
           annotation_colors = ann_colors,
           border_color = NA, 
           fontsize = 10, 
           scale="row", 
           fontsize_row = 10, 
           height=20)
}


plotHM(norm.ileum.ape.DEG.baseFiltered)

# 1. Match ID with gene symbols and remove duplicates
## Explore the grch37 table loaded by the annotables library
grcm38
# extract protein coding
ref.coding <- grcm38 %>% filter(biotype == "protein_coding")
ref.coding <- ref.coding$ensgene

# Write a function to extract non_duplicated significant gene set (grch37)
#----------------------- Function begins -----------------------
SigGenSet <- function(res_table_tb, p, lfc, lfc2){
  ## Return the IDs for the gene symbols in the DE results
  ids <- grcm38[grcm38$symbol %in% res_table_tb$gene, ]
  ## The gene names can map to more than one Ensembl ID (some genes change ID over time), so we need to remove duplicate IDs. Duplicated returns a logical vector indicating which rows of a data.table are duplicates of a row with smaller subscripts.
  non_duplicates <- which(duplicated(ids$symbol) == FALSE)
  ids <- ids[non_duplicates, ] 
  ## Merge the IDs with the results 
  res_ids <- inner_join(res_table_tb, ids, by=c("gene"="symbol")) 
  ## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
  bg <- as.character(res_ids$ensgene)
  ## Extract significant results
  sigAll <- as.character(filter(res_ids, padj < p, abs(log2FoldChange) > lfc)$ensgene)
  sigUp <- as.character(filter(res_ids, padj < p & log2FoldChange >= lfc)$ensgene)
  sigDown <- as.character(filter(res_ids, padj < p & log2FoldChange <= lfc2)$ensgene)
  ## Output as a list
  list(sigAll, sigUp, sigDown,bg,res_ids)
}

#----------------------- Function ends -----------------------
## USAGE: output_list_variable <- SigGenSet(res_table_tb, padj_threshold, lfc of up-regualtion, lfc of down-regulation) ##

# Run the function to get significant gene set and background gene set
GS_List <- SigGenSet(resLFC.ape_smarcd2_base.tb, 0.01, 1.5, -1.5)

# Check the filtered gene number
length(GS_List[[1]]) # [[1]]: all the significant genes
length(GS_List[[2]]) # [[2]]: up-regulated
length(GS_List[[3]]) # [[3]]: down-regulated 



# What is the reason of inconsistent genes numberof up and down regulated genes between two methods?
## Some symbols don't have ID

## Function to run GO enrichment analysis by default setting MOUSE

get_ego_all <- function(Gene_Set_List){
  enrichGO(gene = Gene_Set_List[[1]], 
           universe = ref.coding,
           keyType = "ENSEMBL",
           OrgDb = org.Mm.eg.db, 
           ont = "BP", 
           pAdjustMethod = "BH", 
           qvalueCutoff = 0.05, 
           readable = TRUE)
}
get_ego_up <- function(Gene_Set_List){
  enrichGO(gene = Gene_Set_List[[2]], 
           universe = ref.coding,
           keyType = "ENSEMBL",
           OrgDb = org.Mm.eg.db, 
           ont = "BP", 
           pAdjustMethod = "BH", 
           qvalueCutoff = 0.05, 
           readable = TRUE)
}
get_ego_down <- function(Gene_Set_List){
  enrichGO(gene = Gene_Set_List[[3]], 
           universe = ref.coding,
           keyType = "ENSEMBL",
           OrgDb = org.Mm.eg.db, 
           ont = "BP", 
           pAdjustMethod = "BH", 
           qvalueCutoff = 0.05, 
           readable = TRUE)
}

ego <- get_ego_all(GS_List)

ego_up <- get_ego_up(GS_List)

ego_down <- get_ego_down(GS_List)



## Output results from GO analysis to a table
#write.csv(data.frame(ego_IvC), "results/functional_analysis/ClusterProfiler/GO_IvC.csv")
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.

## Dotplot 
dotplot(ego_up, showCategory=20, font.size=15)
dotplot(ego_down, showCategory=20, font.size=15)

## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
#enrichMap(ego, n=50, vertex.label.font=6)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
fc <- resLFC.ape_smarcd2_base.tb$log2FoldChange

#names(CeKW_fc) <- sigCeKW$gene

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=fc, 
         vertex.label.font=6)

## If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
#CeKW_fc <- ifelse(CeKW_fc > 2, 2, CeKW_fc)
#CeKW_fc <- ifelse(CeKW_fc < -2, -2, CeKW_fc)

#cnetplot(ego, 
#        categorySize="pvalue", 
#       showCategory = 5, 
#      foldChange=CeKW_fc, 
#     vertex.label.font=6)

# If you are interested in significant processes that are not among the top five, you can subset your ego dataset to only display these processes:

## Subsetting the ego results without overwriting original `ego` variable
#ego2 <- ego

#ego2@result <- ego@result[c(1,3,4,8,9),]

## Plotting terms of interest
#cnetplot(ego2, 
#         categorySize="pvalue", 
#         foldChange=fc, 
#         showCategory = 5, 
#         vertex.label.font=6)

# GSEA Gene set enrichment analysis using clusterProfiler and Pathview
## Remove any NA values
res_entrez <- filter(GS_List[[5]], entrez != "NA")

## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrez) == F), ]

## Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange

## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrez

## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)
length(foldchanges)
sum(foldchanges>0) 
sum(foldchanges<0) 

## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "mmu", # supported organisms listed https://www.genome.jp/kegg/catalog/org_list.html
                    nPerm = 1000, # default number permutations
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.03, # padj cutoff value
                    verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result
View(gseaKEGG_results)
## Write GSEA results to file
gseaKEGG_results <- gseaKEGG_results %>% arrange(desc(NES))
write.csv(gseaKEGG_results, "results/functional_analysis/gsea/gsea_kegg.csv", quote=F)
View(gseaKEGG_results)

# Covert ENTREZID to genesymbol if you want #
gseaKEGG_symbol <- setReadable(gseaKEGG, org.Mm.eg.db, keyType = "ENTREZID")
gseaKEGG_results_symbol <- gseaKEGG_symbol@result %>% as_tibble()
write.csv(gseaKEGG_results_symbol, "results/functional_analysis/gsea/gsea_kegg_symbol.csv", quote=F)
gseaKEGG_results_symbol %>% View()
# filtering the signaling pathway you want
#gseaKEGG_results_symbol_filter <- gseaKEGG_results_symbol[-grep("virus", gseaKEGG_results_symbol$Description),]

## Plot the GSEA plot for a single enriched pathway
gseaplot2(gseaKEGG, geneSetID = 'mmu04110', title = "Cell cycle")
gseaplot2(gseaKEGG, geneSetID = 'mmu03008', title = "Ribosome biogenesis in eukaryotes")

nrow(gseaKEGG_results)
p1 <- gseaplot2(gseaKEGG, geneSetID = gseaKEGG_results[1:5, "ID"], pvalue_table = FALSE)
p2 <- gseaplot2(gseaKEGG, geneSetID = gseaKEGG_results[6:10, "ID"], pvalue_table = FALSE)
p3 <- gseaplot2(gseaKEGG, geneSetID = gseaKEGG_results[11:15, "ID"], pvalue_table = FALSE)
p4 <- gseaplot2(gseaKEGG, geneSetID = gseaKEGG_results[16:19, "ID"], pvalue_table = FALSE)
p5 <- gseaplot2(gseaKEGG, geneSetID = gseaKEGG_results[20:26, "ID"], pvalue_table = FALSE)
ggarrange(p1, p2, p3, p4, p5,
          labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow = 2)

ggsave(p1,filename=paste("results/functional_analysis/gsea/p1",".png",sep=""))
ggsave(p2,filename=paste("results/functional_analysis/gsea/p2",".png",sep=""))
ggsave(p3,filename=paste("results/functional_analysis/gsea/p3",".png",sep=""))
ggsave(p4,filename=paste("results/functional_analysis/gsea/p4",".png",sep=""))
ggsave(p5,filename=paste("results/functional_analysis/gsea/p5",".png",sep=""))
ggsave(p6,filename=paste("results/functional_analysis/gsea/p6",".png",sep=""))
ggsave(p7,filename=paste("results/functional_analysis/gsea/p7",".png",sep=""))


write.csv(res_entrez, "results/functional_analysis/res_entrez.csv", quote=F)


#### http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html

## Output images for a single significant KEGG pathway
pathview(gene.data = foldchanges,
         kegg.dir = "results/functional_analysis/gsea/",
         kegg.native = TRUE,
         pathway.id = gseaKEGG_results[1:41, "ID"],
         species = "hsa",
         limit = list(gene = 1, # value gives the max/min limit for foldchanges
                      cpd = 1))

## Output images for all significant KEGG pathways
#get_kegg_plots <- function(x) {
# pathview(gene.data = foldchanges, pathway.id = gseaKEGG_results$ID[x], species = "mmu", 
#         limit = list(gene = 2, cpd = 1))
#}

#purrr::map(1:length(gsea_results$ID), get_kegg_plots)



# Output the versions of all tools used in the DE analysis:
sessionInfo()

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

resLFC.ape_smarcd2_base.tb <- resLFC.ape_smarcd2_base.tb %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 1)


# plot Volcano 
plotVol(resLFC.ape_smarcd2_base.tb)

# Output the versions of all tools used in the DE analysis:
sessionInfo()




