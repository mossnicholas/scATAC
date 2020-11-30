library(ggplot2)
library(scales)
library(RColorBrewer)

# generates overlaid nucleosomal cluster plots, with associated cell % representation
ClusterNucs = function(coverage_matrix, # matrix - columns are cells, rows are single-base coverages
                       clusters = 6, # number of cluster, Kmeans
                       center = FALSE # center coverages?
                       ){
  if (center == TRUE){
    coverage_matrix = scale(coverage_matrix, center=TRUE, scale=TRUE)
  } else {coverage_matrix = scale(coverage_matrix, center=FALSE, scale=TRUE)}
  
  coverage_matrix = as.data.frame(coverage_matrix[,colSums(is.na(coverage_matrix)) == 0])
  t = t(coverage_matrix)
  result = kmeans(t, clusters)
  
  cl_order = as.numeric(names(sort(table(result$cluster), decreasing = TRUE)))
  cluster_coverages = as.data.frame(matrix(ncol = clusters, nrow = 600))
  ratios = c()
  for (i in 1:clusters){
    cluster = result$cluster[result$cluster == cl_order[i]]
    cluster_coverages[,i] = smooth.spline(rowMeans(coverage_matrix[,names(cluster)]), spar = 0.7)$y
    ratio = round(length(cluster)/length(result$cluster)*100, 1)
    ratios = c(ratios, ratio)
  }
  colnames(cluster_coverages) = ratios
  melt = suppressWarnings(reshape2::melt(cluster_coverages))
  melt$bp = rep(-299:300, clusters)
  colnames(melt) = c("percent", "coverage", "bp")
  
  g = ggplot(melt, aes(x = bp, y = coverage, group = percent)) + 
        geom_line(aes(color = percent)) + 
        scale_color_manual(values=rev(brewer.pal(clusters, "RdPu")), name = "Cell %") +
        theme_classic() + 
        scale_x_continuous(breaks = seq(-300, 300, by = 300)) + 
        ylab("Normalized coverage") +
        xlab("Distance from motif (bp)")
  return(g)
}

# generates single-cell nucleosome profile heatmaps, ordered by NFR metric for assessing bimodality 
HeatmapNucs = function(coverage_matrix,
                       subset = FALSE){
  coverage_matrix = as.data.frame(coverage_matrix)
  
  if (subset == TRUE){
    if (ncol(coverage_matrix)>500){
      coverage_matrix = coverage_matrix[, sample(colnames(coverage_matrix), 500)]
    }
  }
  
  coverage_matrix = coverage_matrix[, colSums(coverage_matrix != 0) > 0]
  central_region = colSums(coverage_matrix[275:325,])
  flanking_regions = colSums(coverage_matrix[c(175:274, 326:425),])
  
  nfr = flanking_regions/200 - central_region/51
  nfr_rep = sapply(nfr, function (x) rep(x, 600))
  cellMaxs = apply(coverage_matrix, 2, function(x) max(x, na.rm = TRUE))
  normalized = sweep(coverage_matrix, 2, cellMaxs, FUN = '/')
  
  melt_oc = reshape2::melt(normalized)
  melt_oc$bp = rep(-300:299, ncol(coverage_matrix))
  melt_oc$nfr = as.vector(nfr_rep)
  
  colnames(melt_oc) = c("cell", "occupancy", "bp", "nfr")
  g = ggplot(melt_oc, aes(bp, reorder(cell, nfr), fill = occupancy)) + 
    geom_tile() +
    scale_fill_gradient2(low="white",high="blue", limits = c(0,1), breaks = c(0,1), oob = scales::squish) +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    theme(text = element_text(size = 20)) +
    scale_x_continuous(breaks = seq(-300, 300, by = 300)) +
    theme(legend.position = "none") + 
    theme(plot.title = element_text(size = 1))
  return(g)
}


