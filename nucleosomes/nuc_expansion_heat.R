library(stringr)

fpath = "~/nmoss/nuc_position/sc_atac/"
sc_files = list.files(path = fpath, 
                      pattern = ".*mono.*promoter.*", full.names = TRUE)
names = list.files(path = fpath, 
                   pattern = ".*mono.*promoter.*")

motifs = read.table(paste(fpath, "chip_nuc"))[,1]
metadata = read.table(paste(fpath, "revision/chrom_remod_41_plates/barcode_guide_data_sorted_w_frags.txt"))
ko_genes = unique(metadata$gene)
barcodes = substr(metadata$DNA_BARCODE, 1,24)
print(table(pbmc_two$gene))

motifs_n = length(motifs)
ko_n = length(ko_genes)

nfr_matrix = matrix(NA, nrow = ko_n, ncol = motifs_n)
pval_matrix = matrix(NA, nrow = ko_n, ncol = motifs_n)
expansion_matrix = matrix(NA, nrow = ko_n, ncol = motifs_n)
distance_matrix = matrix(NA , nrow = ko_n, ncol = motifs_n)

# generate local scATAC nucleosome heatmaps around motifs
for (a in 1:length(motifs)){
  retool = FALSE
  for (b in 1:length(ko_genes)){
    ko_barcodes = barcodes[pbmc$gene==ko_genes[b]]
    occupancies = matrix(NA, nrow = 600, ncol = 2000)
    count = 1
    for (i in 1:length(sc_files)){
      barcode = substring(names[i], 1, 24)
      
      if (is.na(str_match(names[i], motifs[a])[1])==FALSE && barcode %in% ko_barcodes){
        table = read.table(sc_files[i])
        occupancies[,count] = table$V2[1:600]
        count = count + 1
        dim(occupancies)
      }
    } 
    occupancies = t(na.omit(t(occupancies)))
    
    if (retool = TRUE){
      occ_count
      final_occ = matrix(NA, nrow = 600, ncol = ncol(occupancies))
      for (i in seq(1, ncol(occupancies), 2)){
        final_occ[,occ_count] = rowMeans(cbind(occupancies[,i], occupancies[,i+1]))
        occ_count = occ_count + 1
      }
    }
    
    final_df = as.data.frame(occupancies)
    final_df = final_df[, colSums(final_df != 0) > 0]
    
    central_region = colSums(final_df[275:325,])
    flanking_regions = colSums(final_df[c(175:274, 326:425),])
    
    nfr = flanking_regions/200 - central_region/51
    nfr_rep = sapply(nfr, function (x) rep(x, 600))
    total = colMeans(final_df)
    total_rep = sapply(total, function (x) rep(x, 600))
    
    melt_oc = melt(final_df)
    melt_oc$bp = rep(-300:299, ncol(final_df))
    melt_oc$nfr = as.vector(nfr_rep)
    melt_oc$total = as.vector(total_rep)
    colnames(melt_oc) = c("cell", "occupancy", "bp", "nfr", "total")
    cutoff = quantile(melt_oc$occupancy, 0.99) # test rough outlier cutoff for upper read bound
    melt_oc$normalized = melt_oc$occupancy/cutoff
    ggplot(melt_oc, aes(bp, reorder(cell, nfr), fill = normalized)) + 
      geom_tile() +
      scale_fill_gradient2(low="white",high="blue", limits = c(0,1), breaks = c(0,1), oob = scales::squish) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      theme(axis.title.x=element_blank()) +
      scale_x_continuous(breaks = seq(-300, 300, by = 300)) +
      theme(legend.title = element_blank())
      
    
    coverage = rowMeans(occupancies)
    coverage = as.data.frame(coverage)
    coverage$bp = c(-300:299)
    coverage = coverage[,c(2,1)]
    
    coverage$smoothed = (smooth.spline(coverage$coverage, spar = 0.7))$y
    upstream = coverage[1:300,]
    downstream = coverage[301:600,]
    downstream_max = downstream$bp[downstream$smoothed==max(downstream$smoothed)]
    upstream_max = upstream$bp[upstream$smoothed==max(upstream$smoothed)]
    zeros_matrix[b,a] = 0
    
    if (downstream_max == 0){
      for (position in 1:(length(downstream$smoothed)-1)){
        if (downstream$smoothed[position+1]>downstream$smoothed[position]){break}
      } 
      max = max(downstream$smoothed[position:length(downstream$smoothed)])
      if (position == length(downstream$smoothed)-1){
        max = max(downstream$smoothed)
        zeros_matrix[b,a] = 1
      }
      downstream_max = downstream$bp[downstream$smoothed==max]
    }
    
    if (upstream_max == -1){
      rev = rev(upstream$smoothed)
      for (position in 1:(length(rev)-1)){
        if (rev[position+1]>rev[position]){break}
      } 
      max = max(rev[position:length(rev)]) 
      if (position == length(upstream$smoothed)-1){
        max = max(upstream$smoothed)
        zeros_matrix[b,a] = 1
      }
      upstream_max = upstream$bp[upstream$smoothed==max]
    }

    pdf(paste(motifs[a], ko_genes[b], "enhancer.pdf", sep = "_"), width = 6, height = 5)
    ylim = c(min(coverage$smoothed), max(rbind(coverage$smoothed)))
    plot(coverage$smoothed, type = "l", col="darkblue", xaxt = "n", ylim=ylim, ylab="Normalized coverage",
         xlab=sprintf("Distance from %s binding site", motifs[a]), 
         main=sprintf("KO=%s", ko_genes[b]), cex=15, cex.axis=1.5, cex.lab=1.5)
    lines(coverage$smoothed, col="darkblue", lwd =2)
    abline(v=upstream_max+300,col="darkblue",lty=3)
    abline(v=downstream_max+300,col="darkblue",lty=3)
    axis(1, at=c(0, 300, 600), labels=c(-300, 0, 300), cex.axis=1.5)
    dev.off()

    central_region = coverage$smoothed[275:325]
    flanking_regions = coverage$smoothed[c(175:274, 326:425)]
    nfr = sum(flanking_regions)/length(flanking_regions) - sum(central_region)/length(central_region)
    nfr_matrix[b,a] = nfr

  }
}
