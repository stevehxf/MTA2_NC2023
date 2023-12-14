

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

# Run the function to get significant gene set and background gene set
GS.list <- resLFC.normal.tb.list %>% lapply(SigGenSet, 0.05, 1, -1)

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

#### GSEA ####
# Gene set enrichment analysis using clusterProfiler and Pathview 
foldchanges.list <- GS.list %>% lapply(function(x){
        ## Remove any NA values
        res_entrez <- filter(x[[5]], entrez != "NA")     
        ## Remove any Entrez duplicates
        res_entrez <- res_entrez[which(duplicated(res_entrez$entrez) == F), ]
        ## Extract the foldchanges
        foldchange <- res_entrez$log2FoldChange
        ## Name each fold change with the corresponding Entrez ID
        names(foldchange) <- res_entrez$entrez
        ## Sort fold changes in decreasing order
        foldchange <- sort(foldchange, decreasing = TRUE)
        x <- foldchange
})

## GSEA using gene sets from KEGG pathways
gseaKEGG.list <- foldchanges.list %>% 
        lapply(gseKEGG,
               organism = "mmu", # supported organisms listed https://www.genome.jp/kegg/catalog/org_list.html
               nPerm = 1000, # default number permutations
               minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
               pvalueCutoff = 0.1, # padj cutoff value
               verbose = FALSE)


## import the GSEA results
if(!dir.exists("results/functional_analysis/")){
        dir.create("results/functional_analysis/")}
if(!dir.exists("results/functional_analysis/gsea")){
        dir.create("results/functional_analysis/gsea")}

gseaKEGG.list.symbol <- gseaKEGG.list %>% lapply(setReadable, org.Mm.eg.db, keyType = "ENTREZID")

for(i in seq_along(gseaKEGG.list.symbol)){

        gseaKEGG_results_symbol <- gseaKEGG.list.symbol[[i]]@result 
        write.csv(gseaKEGG_results_symbol, paste0("results/functional_analysis/gsea/", names(gseaKEGG.list)[i],"_gsea_kegg_symbol.csv"))
        
}

# filtering the signaling pathway you want
#gseaKEGG_results_symbol_filter <- gseaKEGG_results_symbol[-grep("virus", gseaKEGG_results_symbol$Description),]

## Plot the GSEA plot for a single enriched pathway
gseaplot2(gseaKEGG.list$genotyping_hnf4a_oe_vs_ctrl, geneSetID = c('mmu00830','mmu04975','mmu04974','mmu04979','mmu00071'))
ggsave("figures/gsea_hnf4a_up.pdf", width = 8, height = 10)

#### http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html


# Output the versions of all tools used in the DE analysis:
sessionInfo()
