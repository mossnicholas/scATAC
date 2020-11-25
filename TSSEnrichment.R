library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(data.table)
library(rtracklayer)

# TSS + promoters file to intersect bed
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
transcripts <- unique(transcripts(txdb)) 

promoters <- promoters(transcripts, downstream = 2000, upstream = 2000)
promoters <- keepStandardChromosomes(promoters, pruning.mode = "coarse")
promoters <- dropSeqlevels(promoters, c("chrX", "chrY", "chrM"), pruning.mode = "coarse")
export(promoters, con = "promoters2000.bed", format = "BED")

# ATAC-seq signal-to-noise QC: coveraged at TSS normalized by up/downstream flank baseline
# Great TSS Enrichment scores are considered above 7
TSSEnrichment <- function(TSS_coverage, # single-base read coverage bed file for regions around TSS
                          range = 4000, 
                          flank = 100){
  
  reads <- fread(TSS_coverage)
  num_promoters <- nrow(reads)/range
  tss_cov <- matrix(0, num_promoters, range)
  
  for(i in 0:num_promoters-1){
    tss_cov[i+1,] <- reads$V8[(i*range+1) : (i*range+range)]
  }
  tss_cov <- as.data.frame(tss_cov)
  
  flank_up <- as.numeric(lapply(tss_cov[,1:flank], mean, na.rm = TRUE))
  flank_down <- as.numeric(lapply(tss_cov[,(range-flank):range], mean, na.rm = TRUE))

  avg_flank <- mean(c(flank_up, flank_down))
  tss_cov_norm <- tss_cov/avg_flank
  tss_values <- colMeans(tss_cov_norm)
  tss_score <- mean(tss_values[range/2], tss_values[range/flank+1])
  
  return(tss_score)
}


