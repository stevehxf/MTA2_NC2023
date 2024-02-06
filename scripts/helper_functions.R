## to export gsea plot one by one. 
## user should provide the path of a csv file, where the first column contains the pathways of interest
require("clusterProfiler")
require("enrichplot")
gsea.plot.export <- function(gsea.filter.list, prefix, PATH="results/gsea/fig/"){
        
        
        for (i in seq_along(gsea.filter.list)){
                if(dim(gsea.filter.list[[i]]@result)[1] == 0){
                        gsea.filter.list[[i]] <- NA
                }
        }
        gsea.filter.list <-gsea.filter.list[!is.na(gsea.filter.list)]
        # visulization
        for(i in seq_along(gsea.filter.list)){
                n <- nrow(gsea.filter.list[[i]]@result)
                for (j in 1:n){
                        p <- gseaplot2(gsea.filter.list[[i]], geneSetID = j, title = paste0(prefix,names(gsea.filter.list)[[i]], "_", gsea.filter.list[[i]]$Description[j]))
                        filename <- paste0(PATH,prefix,names(gsea.filter.list)[[i]], "_",gsea.filter.list[[i]]$Description[j],".pdf")
                        ggsave(filename = filename, plot = p)
                }
        }
}
