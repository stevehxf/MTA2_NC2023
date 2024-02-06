

# This script is designed for RNA-seq functional analysis
# Author: Xiaofeng Steve Huang
# reference: HBCtraining
# Date: May 11 2018
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

#### Match ID with gene symbols and remove duplicates ####

# extract protein coding
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
  sigDown <- as.character(filter(res_ids, padj < p & log2FoldChange <= lfc2)$ensgene)
  ## Output as a list
  list(sigAll, sigUp, sigDown,bg,res_ids)
}

#----------------------- Function ends -----------------------
## USAGE: output_list_variable <- SigGenSet(res_table_tb, padj_threshold, lfc of up-regualtion, lfc of down-regulation) ##
#Satb2_vs_cas9
# Run the function to get significant gene set and background gene set
GS_List_satb2_vs_cas9 <- SigGenSet(resLFC.normal.tb.list$sgRNA_satb2_vs_cas9, 0.05, 0.0, -0.0)

# Check the filtered gene number
length(GS_List_satb2_vs_cas9[[1]]) # [[1]]: all the significant genes
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
ego_up_satb2_vs_cas9 <- get_ego_up(GS_List_satb2_vs_cas9)
ego_down_satb2_vs_cas9 <- get_ego_down(GS_List_satb2_vs_cas9)


## Output results from GO analysis to a table
#write.csv(data.frame(ego_IvC), "results/functional_analysis/ClusterProfiler/GO_IvC.csv")
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.

## Dotplot 
up.gene.list <- list(satb2_vs_cas9 = GS_List_satb2_vs_cas9[[2]])
compare.up <- compareCluster(up.gene.list, 
                    fun = "enrichGO", 
                    universe = ref.coding,
                    keyType = "ENSEMBL",
                    OrgDb = org.Mm.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)

down.gene.list <- list(satb2_vs_cas9 = GS_List_satb2_vs_cas9[[3]])
compare.down <- compareCluster(down.gene.list, 
                             fun = "enrichGO", 
                             universe = ref.coding,
                             keyType = "ENSEMBL",
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = TRUE)

p1 <- dotplot(compare.up, showCategory=10, font.size=10)
p2 <- dotplot(compare.down, showCategory=10, font.size=10)
plot_grid(p1,p2, nrow = 2, align = "v", labels = letters[1:2])


## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
#enrichMap(ego, n=50, vertex.label.font=6)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
#fc <- resLFC.highRPKM.tb$log2FoldChange

#names(CeKW_fc) <- sigCeKW$gene

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
#cnetplot(ego, 
#   categorySize="pvalue", 
#   showCategory = 5, 
#   foldChange=fc, 
#   vertex.label.font=6)

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

#### GSEA ####
# Gene set enrichment analysis using clusterProfiler and Pathview 
## Remove any NA values
res_entrez <- filter(GS_List_satb2_vs_cas9[[5]], entrez != "NA")

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
write.csv(gseaKEGG_results_symbol, "results/functional_analysis/gsea/gsea_kegg_satb2_vs_cas9.csv", quote=F)
gseaKEGG_results_symbol %>% View()
# filtering the signaling pathway you want
#gseaKEGG_results_symbol_filter <- gseaKEGG_results_symbol[-grep("virus", gseaKEGG_results_symbol$Description),]

## Plot the GSEA plot for a single enriched pathway
gseaplot2(gseaKEGG, geneSetID = 'mmu04975', title = "Fat digestion and absorption")
gseaplot2(gseaKEGG, geneSetID = 'mmu00564', title = "Glycerophospholipid metabolism")

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

GS_List_mta2_vs_cas9 <- SigGenSet(resLFC.normal.tb.list$sgRNA_mta2_vs_cas9, 0.05, 0.0, -0.0)

# Check the filtered gene number
length(GS_List_mta2_vs_cas9[[1]]) # [[1]]: all the significant genes
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
ego_up_mta2_vs_cas9 <- get_ego_up(GS_List_mta2_vs_cas9)
ego_down_mta2_vs_cas9 <- get_ego_down(GS_List_mta2_vs_cas9)


## Output results from GO analysis to a table
#write.csv(data.frame(ego_IvC), "results/functional_analysis/ClusterProfiler/GO_IvC.csv")
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.

## Dotplot 
up.gene.list <- list(mta2_vs_cas9 = GS_List_mta2_vs_cas9[[2]])
compare.up <- compareCluster(up.gene.list, 
                             fun = "enrichGO", 
                             universe = ref.coding,
                             keyType = "ENSEMBL",
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = TRUE)

down.gene.list <- list(mta2_vs_cas9 = GS_List_mta2_vs_cas9[[3]])
compare.down <- compareCluster(down.gene.list, 
                               fun = "enrichGO", 
                               universe = ref.coding,
                               keyType = "ENSEMBL",
                               OrgDb = org.Mm.eg.db, 
                               ont = "BP", 
                               pAdjustMethod = "BH", 
                               qvalueCutoff = 0.05, 
                               readable = TRUE)

p1 <- dotplot(compare.up, showCategory=10, font.size=10)
p2 <- dotplot(compare.down, showCategory=10, font.size=10)
plot_grid(p1,p2, nrow = 2, align = "v", labels = letters[1:2])


## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
#enrichMap(ego, n=50, vertex.label.font=6)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
#fc <- resLFC.highRPKM.tb$log2FoldChange

#names(CeKW_fc) <- sigCeKW$gene

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
#cnetplot(ego, 
#   categorySize="pvalue", 
#   showCategory = 5, 
#   foldChange=fc, 
#   vertex.label.font=6)

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

#### GSEA ####
# Gene set enrichment analysis using clusterProfiler and Pathview 
## Remove any NA values
res_entrez <- filter(GS_List_mta2_vs_cas9[[5]], entrez != "NA")

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
write.csv(gseaKEGG_results_symbol, "results/functional_analysis/gsea/gsea_kegg_mta2_vs_cas9.csv", quote=F)
gseaKEGG_results_symbol %>% View()

# Output the versions of all tools used in the DE analysis:
sessionInfo()
