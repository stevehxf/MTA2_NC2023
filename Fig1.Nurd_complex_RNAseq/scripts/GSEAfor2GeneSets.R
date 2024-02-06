## To plot GSEA using same expressoin table and two gene sets
## Script by Steve Huang
## path argument: path to output folder of analysis (e.g. PATH/my_analysis.GseaPreranked.1470948568349)
## gene.set argument: name of the gene set (e.g. V$AP1_Q2).
## It is used in a grep command, so multiple matching is possible.
## Also, R regular expressions can be handled, e.g. "IL2[0-9]$"
## Leading "V$" from gene set names are stripped to allow using the grep command.
## In case of multiple grep matches a warning is given and the first option is plotted.
## class.name: the name of the class / variable to which genes have been correlated (e.g. drug-treatment)
## metric.range: the range of the metric; defaults to [min(DEFINED RANGE), max(DEFINED RANGE)]
## enrichment.score.range: the range of the enrichment score; defaults to [min(ENRICHMENT SCORE), max(ENRICHMENT SCORE)]

replotGSEA <- function(path1, path2, gene.set1, gene.set2, class.name1, class.name2, metric.range1, metric.range2,
                       enrichment.score.range1, enrichment.score.range2) 
{
  
  ## GSEA data loading
  loadGSEA <- function(path, gene.set, class.name, metric.range,
                       enrichment.score.range) {
    
    if (missing(path)) {
      stop("Path argument is required")
    }
    if (!file.exists(path)) {
      stop("The path folder could not be found. Please change the path")
    }
    if (missing(gene.set)) {
      stop("Gene set argument is required")
    }
    
    ## Load .rnk data
    path.rnk <- list.files(path = file.path(path, "edb"),
                           pattern = ".rnk$", full.names = TRUE)
    gsea.rnk <- read.delim(file = path.rnk, header = FALSE)
    colnames(gsea.rnk) <- c("hgnc.symbol", "metric")
    if (missing(metric.range)) {
      metric.range <- c(min(gsea.rnk$metric), max(gsea.rnk$metric))
    }  
    
    ## Load .edb data
    path.edb <- list.files(path = file.path(path, "edb"),
                           pattern = ".edb$", full.names = TRUE)
    gsea.edb <- read.delim(file = path.edb,
                           header = FALSE, stringsAsFactors = FALSE)
    gsea.edb <- unlist(gsea.edb)
    gsea.metric <- gsea.edb[grep("METRIC=", gsea.edb)]
    gsea.metric <- unlist(strsplit(gsea.metric, " "))
    gsea.metric <- gsea.metric[grep("METRIC=", gsea.metric)]
    gsea.metric <- gsub("METRIC=", "", gsea.metric)
    gsea.edb <- gsea.edb[grep("  <DTG", gsea.edb)]
    
    # Select the right gene set
    if (length(gsea.edb) == 0) {
      stop(paste("The gene set name was not found, please provide",
                 "a correct name"))
    }
    if (length(grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), ""), gsea.edb)) > 1) {
      warning(paste("More than 1 gene set matched the gene.set",
                    "argument; the first match is plotted"))
    }
    gsea.edb <- gsea.edb[grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), ""), gsea.edb)[1]]
    
    # Get template name
    gsea.edb <- gsub(".*TEMPLATE=(.*)", "\\1", gsea.edb)
    gsea.edb <- unlist(strsplit(gsea.edb, " "))
    gsea.template <- gsea.edb[1]
    
    # Get gene set name
    gsea.gene.set <- gsea.edb[2]
    gsea.gene.set <- gsub("GENESET=gene_sets.gmt#", "", gsea.gene.set)
    
    
    # Get hit indices
    gsea.edb <- gsea.edb[grep("HIT_INDICES=", gsea.edb):length(gsea.edb)]
    gsea.hit.indices <- gsea.edb[seq_len(grep("ES_PROFILE=", gsea.edb) - 1)]
    gsea.hit.indices <- gsub("HIT_INDICES=", "", gsea.hit.indices)
    gsea.hit.indices <- as.integer(gsea.hit.indices)
    
    # Get ES profile
    gsea.edb <- gsea.edb[grep("ES_PROFILE=", gsea.edb):length(gsea.edb)]
    gsea.es.profile <- gsea.edb[seq_len(grep("RANK_AT_ES=", gsea.edb) - 1)]
    gsea.es.profile <- gsub("ES_PROFILE=", "", gsea.es.profile)
    gsea.es.profile <- as.numeric(gsea.es.profile)
    
    # Set enrichment score range
    if (missing(enrichment.score.range)) {
      enrichment.score.range <- c(min(gsea.es.profile), max(gsea.es.profile))
    }
    ## Output as a list
    list(gsea.hit.indices, gsea.es.profile, enrichment.score.range,gsea.rnk, metric.range, gsea.metric)
  }
  
  # assign 
  PATH1_GSEA <- loadGSEA(path=path1, gene.set=gene.set1, class.name=class.name1, metric.range=metric.range1, enrichment.score.range=enrichment.score.range1)
  PATH2_GSEA <- loadGSEA(path=path2, gene.set=gene.set2, class.name=class.name2, metric.range=metric.range2, enrichment.score.range=enrichment.score.range2)
  
  gsea.hit.indices1 <- PATH1_GSEA[[1]]
  gsea.es.profile1 <- PATH1_GSEA[[2]]
  enrichment.score.range1 <-  PATH1_GSEA[[3]]
  gsea.rnk1 <- PATH1_GSEA[[4]]
  metric.range1 <- PATH1_GSEA[[5]]
  gsea.metric1 <- PATH1_GSEA[[6]]
  
  gsea.hit.indices2 <- PATH2_GSEA[[1]]
  gsea.es.profile2 <- PATH2_GSEA[[2]]
  enrichment.score.range2 <- PATH2_GSEA[[3]]
  gsea.rnk2 <-  PATH2_GSEA[[4]]
  metric.range2 <- PATH2_GSEA[[5]]
  gsea.metric2 <- PATH2_GSEA[[6]]
  
  # Save default for resetting
  def.par <- par(no.readonly = TRUE)
  
  # Create a new device of appropriate size
  #dev.new(width = 3, height = 3)
  
  # Create a division of the device
  gsea.layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(1, 0.1, 0.1, 0.3))
  layout.show(gsea.layout)
  
  # Create plots
  par(mar = c(0, 5, 1, 2))
  plot(frame.plot=TRUE,
       c(1, gsea.hit.indices1, length(gsea.rnk1$metric)),
       c(0, gsea.es.profile1, 0), type = "l", col = "firebrick2", lwd = 3, xaxt = "n",
       xaxs = "i", xlab = "", ylab = "",
       ylim = c(min(enrichment.score.range1,enrichment.score.range2)-0.1, max(enrichment.score.range1,enrichment.score.range2)+0.1)
       
  )
  
  abline(h = 0, lty = 2)
  
  lines(c(1, gsea.hit.indices2, length(gsea.rnk2$metric)),
        c(0, gsea.es.profile2, 0), type = "l", col = "dodgerblue2", lwd = 3, xaxt = "n",
        xaxs = "i", xlab = "", ylab = ""
  )
  
  
  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(1, length(gsea.rnk1$metric)))
  abline(v = gsea.hit.indices1, lwd = 1, col = "firebrick2")
  
  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(1, length(gsea.rnk2$metric)))
  abline(v = gsea.hit.indices2, lwd = 1, col = "dodgerblue2")
  
  par(mar = c(1, 5, 0, 2))
  rank.metric <- rle(round(gsea.rnk1$metric, digits = 2))
  plot(gsea.rnk1$metric, type = "n", xaxs = "i",
       xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk1$metric)),
       ylim = metric.range1, yaxs = "i",
       ylab = if(gsea.metric1 == "None") {"Ranking metric"} else {gsea.metric1},
       panel.first = abline(h = seq(metric.range1[1] / 2,
                                    metric.range1[2] - metric.range1[1] / 4,
                                    metric.range1[2] / 2), col = "gray95", lty = 2))
  
  barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
          xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk1$metric)),
          ylim = c(-1, 1), yaxs = "i", width = rank.metric$lengths, border = NA,
          ylab = ifelse(gsea.metric1 == "None", "Ranking metric", gsea.metric1), space = 0, add = TRUE)
  box()
  
  # Reset to default
  par(def.par)
}

replotGSEA(path1="GSEA/results/colon_enriched_lfc2.5", path2="GSEA/results/ileum_enriched_lfc2.5", gene.set1="colon_enriched", gene.set2="ileum_enriched", class.name1=" ", class.name2=" ")


replotGSEA(path1="GSEA/results/ileum_enriched_normal_ranked", path2="GSEA/results/colon_enriched_normal_ranked", gene.set1="ileum_enriched", gene.set2="colon_enriched", class.name1=" ", class.name2=" ")
