# This script is designed for differential expression (DE) analysis using DEseq2 as well as data visualization
# apeglm shrunkage used in this script

## USAGE & SETUP ----
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

library(EnhancedVolcano)
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
library(clusterProfiler)
library(annotables)
#library(enrichplot)
#library(Rgraphviz)

# load feature counts
data.raw <- read.table("data/3D_cleancounts_paper.txt", header=T, row.names=1, check.names = F)  # Put your unnormalized counts file in data directory
# extract protein_coding gene counts from the raw counts
coding_genes <- grcm38 %>% dplyr::filter(biotype == "protein_coding")
coding_genes <- coding_genes$symbol
idx <- row.names(data.raw) %in% c(coding_genes,"Inava")
data <- data.raw[idx,]
write.table(data, file="data/data.txt", sep="\t", quote=F, col.names=NA)
# load data
data <- read.table("data/data.txt", row.names=1, header=T)
# load and subset meta
meta <- read.table("meta/meta_paper.txt", row.names=1, header=T) # Put your meta file in meta directory
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
#### ------ above is the function --
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
                              design = ~ sgRNA)  # Comparing factor 


# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE) 

# exploratory data analysis (PCA & heirarchical clustering) - identifying outliers and sources of variation in the data
## PCA ##
p1 <- plotPCA(rld, intgroup="sgRNA", ntop=500, returnData = F) +
  ggtitle ("PCA") + # main title  
  theme_light() +  # Theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size = 1, colour = "black"),
        legend.title=element_blank()
  )+
  #geom_text_repel(aes(label = colnames(data))) +
  #coord_cartesian(ylim = c(-30,40), xlim = c(-60,40)) + # Coordinate system
  scale_color_manual(values = c("#1f78b4","#4daf4a","#e41a1c","gray","#e5d8bd","#a6cee3","#feb709","#fb027f","#8d5287"))+
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

rld_ntop <- ntop(rld, ntop = 2000) # Extract the rlog matrix from the object (top 500 var by default)
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
rld_top_cor <- rld_ntop %>% cor()
head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames
# plot heatmap
p1 <- pheatmap(rld_cor) # Check the heatmap plot # Zoom and copy the vector figure
p2 <- pheatmap(rld_top_cor)
pdf("figures/heatmap_cor.pdf")
p1
dev.off()
pdf("figures/heatmap_top_cor.pdf")
p2
dev.off()
## DEseq2 Analysis -----
# **Optional step** - Re-create DESeq2 dataset if the design formula has changed after QC analysis in include other sources of variation
dds <- DESeqDataSetFromMatrix(countData = data,  
                              colData = meta,
                              design = ~ sgRNA)  # Comparing factor 

# Extract normalized counts, check the size factors and total counts before and after normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds) # check the size factors
normalized_counts <- counts(dds, normalized=TRUE)

write.table(normalized_counts, file="results/normalization/normalized_counts_nurd.txt", sep="\t", quote=F, col.names=NA)

## Take a look at the total number of reads for each sample using
colSums(counts(dds)) # Sum the unnormalized column values
colsum <- colSums(counts(dds)) 
write.csv(colsum, file="results/colSum.csv") # Export

colSums(normalized_counts) # Sum the normalized column values
ncolsum <- colSums(normalized_counts) 
write.csv(ncolsum, file = "results/nor_colSum.csv")

# change reference level
dds$sgRNA <- dds$sgRNA %>% relevel(ref = "cas9")
# Differential expression analysis
dds <- DESeq(dds)
sizeFactors(dds)

# Examine dispersion estimates plot to ensure your data is a good fit for the DESeq2 model
plotDispEsts(dds) 

# Create contrasts to perform Wald testing on the shrunken log2 foldchanges between specific conditions:
# check the coef for shurnkage
resultsNames(dds)
saveRDS(dds, "rds/dds.rds")
#dds <- readRDS("rds/dds.rds")
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

## MA plot
#DESeq2::plotMA(res.list$genotyping_satb2ko_vs_wt_colon, alpha = 0.01, main = "unshrunken")
#DESeq2::plotMA(resLFC.list$genotyping_satb2ko_vs_wt_colon, alpha = 0.01, main = "apeglm")

## Summarize results
summary(resLFC.list$sgRNA_satb2_vs_cas9)

## DEG ----
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



## GSEA to show enrichment of satb2 signature ----
# load signature of ileum/colon
#ileum.signature <- read.table("gene_set/ileum_enriched.grp", header = TRUE)
#names(ileum.signature) <- "gene"
satb2.up.signature <- DEG.normal.list$sgRNA_satb2_vs_cas9 %>% filter(log2FoldChange > 2)
satb2.down.signature <- DEG.normal.list$sgRNA_satb2_vs_cas9 %>% filter(log2FoldChange < -2)
gsea.term <- rbind(#data.frame(gs_name="ileum.signature", symbol_gene=ileum.signature$gene),
                   data.frame(gs_name="satb2up.signature", symbol_gene=satb2.up.signature$gene),
                   data.frame(gs_name="satb2down.signature", symbol_gene=satb2.down.signature$gene))

## fold change sorting
fc.LFC.list <- 
  resLFC.normal.tb.list %>% lapply(function(x){
    res <- x
    fc <- res$log2FoldChange
    names(fc) <- res$gene
    fc <- sort(fc, decreasing = TRUE)
    
  })

str(fc.list)

## GSEA
gseaList <- fc.LFC.list %>% 
  lapply(GSEA, TERM2GENE= gsea.term)

gseaRes <- gseaList %>% 
  lapply(function(x){
    res <- x@result %>% 
    as_tibble()
    res <- res[,c(1,5,7)]
    return(res)
    })
idx <- gseaRes %>% lapply(nrow)
idx <- idx > 0
gseaList.fil <- gseaList[idx]

source("scripts/helper_functions.R")
gsea.plot.export(gseaList.fil, prefix = "", PATH = "figures/")

## heatmap visualization ----
# 2. Heatmap (Filter with RPKM counts)
# write a function to plot heatmap using the same arguments

plotHM <- function(data,title=""){
  pheatmap(data.frame(data), 
           main = title,
           #color = colorRampPalette(c("cornflowerblue","beige","brown1"))(50),
           #color = magma(50),
           color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(40),
           breaks = seq(-2.0,2.0,by=0.1),
           scale="row", 
           cluster_rows = T, cluster_cols = F,
           show_colnames = T, show_rownames = T,
           border_color = NA, 
           fontsize = 10, 
           fontsize_row = 10, 
           height=20)
}

feature.to.plot <- DEG.normal.list$sgRNA_satb2_vs_cas9$gene  
sample.to.plot <- c("Ctrl_1","Ctrl_2","Satb2_KO1","Satb2_KO2","Mta2_KO1","Mta2_KO2","Chd4_KO1","Chd4_KO2",
                    "Gatad2_KO1","Gatad2_KO2","Smarcd2_KO1","Smarcd2_KO2","Smarca4_KO1","Smarca4_KO2","Smarca5_KO1","Smarca5_KO2",
                    "Ctbp2_KO1","Ctbp2_KO2")

# extract vst assay 
# Transform counts for data visualization
vsd <- vst(dds, blind=F)
vsd.df <- vsd %>% assay() %>% data.frame()

## subset a dataframe based on a given gene list
vsd.df %>% head()
## which row exist in the given gene list?
idx <- rownames(vsd.df) %in% feature.to.plot
vsd.df.sub <- vsd.df[idx, sample.to.plot]
## plot a heatmap
p1 <- vsd.df.sub %>% plotHM()
p1
pdf("figures/heatmap_satb2_DEG.pdf")
p1
dev.off()

## volcano plot ----
## adjust ceiling value
ceiling <- 100
table_tb <- resLFC.normal.tb.list$sgRNA_mta2_vs_cas9
idx <- which(-log10(table_tb$padj) > ceiling)
table_tb[idx,]$padj <- 10^(-ceiling)
EnhancedVolcano(table_tb,
                lab = table_tb$gene,
                x = 'log2FoldChange',
                y = 'padj',    
                pCutoff = 10e-1,
                FCcutoff = 1,
                pointSize = 3.0,
                #ylim = c(0, -log10(10e-50)),
                labSize = 6.0)
ggsave("figures/vol.pdf")


## Obtain logical vector regarding whether padj values are less than 0.05 and fold change greater than 1.5 in either direction

resLFC_normal_sgRNA_mta2_vs_cas9.tb <- resLFC.normal.tb.list$sgRNA_mta2_vs_cas9 %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 1)



# plot Volcano (The picture showed in manuscript)
plotVol(resLFC_normal_sgRNA_mta2_vs_cas9.tb)

# We could label those dots with the gene name on the Volcano plot using geom_text_repel().
# Following is to generate gene labeled volcano plot.
## (1) Arrange the padj column and create a new column to indicate which genes to label
resLFC_normal_sgRNA_mta2_vs_cas9.tb <- resLFC_normal_sgRNA_mta2_vs_cas9.tb %>% arrange(padj) %>% mutate(genelabels = "")

resLFC_normal_sgRNA_mta2_vs_cas9.tb$genelabels[1:200] <- resLFC_normal_sgRNA_mta2_vs_cas9.tb$gene[1:200]

## (or) Arrange the lfc column and create a new column to indicate which genes to lablel

resLFC_normal_sgRNA_mta2_vs_cas9.tb <- resLFC_normal_sgRNA_mta2_vs_cas9.tb %>%
  arrange(log2FoldChange) %>%
  mutate(genelabels = "")
resLFC_normal_sgRNA_mta2_vs_cas9.tb$genelabels[1:10] <- resLFC_normal_sgRNA_mta2_vs_cas9.tb$gene[1:10]
resLFC_normal_sgRNA_mta2_vs_cas9.tb <- resLFC_normal_sgRNA_mta2_vs_cas9.tb %>%
  arrange(desc(log2FoldChange))
resLFC_normal_sgRNA_mta2_vs_cas9.tb$genelabels[1:10] <- resLFC_normal_sgRNA_mta2_vs_cas9.tb$gene[1:10]

## 2nd. Add the lables using geom_text_repel


plotVol(resLFC_normal_sgRNA_mta2_vs_cas9.tb, "cas9") +
  geom_text_repel(aes(x = log2FoldChange,
                      y = -log10(padj),
                      label = genelabels))

## 3rd. Move the points out of the range to the ceiling.
#res_tableCoKW_tb_padj_to_ceiling <- res_tableCoKW_tb
#idx <- which(-log10(res_tableCoKW_tb_padj_to_ceiling$padj)>50)
#res_tableCoKW_tb_padj_to_ceiling[idx,]$padj = 1e-51
#plotVol(res_tableCoKW_tb_padj_to_ceiling, "Colon")


##Pathway Analysis
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(annotables)
library(enrichplot)
library(Rgraphviz)
library(cowplot)
ref.coding <- grcm38 %>% filter(biotype == "protein_coding")
ref.coding <- ref.coding$ensgene %>% as.character()

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
  sigDown <- as.character(filter(res_ids, padj < p & log2FoldChange <= lfc)$ensgene)
  ## Output as a list
  list(sigAll, sigUp, sigDown,bg,res_ids)
}


## Function ends ##
## USAGE: output_list_variable <- SigGenSet(res_table_tb, padj_threshold, lfc of up-regualtion, lfc of down-regulation) ##

# Run the function to get significant gene set and background gene set
GS_List_mta2VSctrl <- SigGenSet(resLFC.normal.tb.list$sgRNA_mta2_vs_cas9, 0.05, 1, -1)
GS_List_satb2VSctrl <- SigGenSet(resLFC.normal.tb.list$sgRNA_satb2_vs_cas9, 0.05, 1, -1)
# Check the filtered gene number
length(GS_List_mta2VSctrl[[1]]) # [[1]]: all the significant genes

# What is the reason of inconsistent genes number of up and down regulated genes between two methods?
## Some symbols don't have ID

#### GO Term over representation test ####
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

#ego <- get_ego_all(GS_List_ileumVScolon)
ego_up_mta2VSctrl <- get_ego_up(GS_List_mta2VSctrl)
ego_down_mta2VSctrl <- get_ego_down(GS_List_mta2VSctrl)
ego_up_mta2VSctrl.sim <- ego_up_mta2VSctrl %>% simplify(cutoff=0.2)
ego_down_mta2VSctrl.sim <- ego_down_mta2VSctrl %>% simplify(cutoff=0.2)

p1 <- dotplot(ego_up_mta2VSctrl.sim, showCategory=20, font.size=10)
p2 <- dotplot(ego_down_mta2VSctrl.sim, showCategory=20, font.size=10)
cowplot::plot_grid(p1,p2, nrow = 2, labels = letters[1:2])



## Output results from GO analysis to a table
#write.csv(data.frame(ego_IvC), "results/functional_analysis/ClusterProfiler/GO_IvC.csv")
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.

## Dotplot 
up.gene.list <- list(mta2VSctrl = GS_List_mta2VSctrl[[2]],
                     satb2VSctrl = GS_List_satb2VSctrl[[2]])
compare.up <- compareCluster(up.gene.list, 
                             fun = "enrichGO", 
                             universe = ref.coding,
                             keyType = "ENSEMBL",
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = TRUE)

down.gene.list <- list(mta2VSctrl = GS_List_mta2VSctrl[[3]],
                       satb2VSctrl = GS_List_satb2VSctrl[[3]])
compare.down <- compareCluster(down.gene.list, 
                               fun = "enrichGO", 
                               universe = ref.coding,
                               keyType = "ENSEMBL",
                               OrgDb = org.Mm.eg.db, 
                               ont = "BP", 
                               pAdjustMethod = "BH", 
                               qvalueCutoff = 0.05, 
                               readable = TRUE)
compare.up.sim <- compare.up %>% simplify(cutoff=0.2)
compare.down.sim <- compare.down %>% simplify(cutoff=0.2)
p1 <- dotplot(compare.up.sim, showCategory=10, font.size=10)
p2 <- dotplot(compare.down.sim, showCategory=10, font.size=10)
plot_grid(p1,p2, nrow = 2, align = "v", labels = letters[1:2])


#### GSEA ####
# Gene set enrichment analysis using clusterProfiler and Pathview 
## Remove any NA values
res_entrez <- filter(GS_List_mta2VSctrl[[5]], entrez != "NA")

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
                    pvalueCutoff = 0.1, # padj cutoff value
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
gseaplot2(gseaKEGG, geneSetID = 'mmu00051', title = "Fructose and mannose metabolism")
gseaplot2(gseaKEGG, geneSetID = 'mmu00270', title = "Cysteine and methionine metabolism")
dotplot(gseaKEGG, x="NES", showCategory=34)

names(gseaKEGG_results)
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
