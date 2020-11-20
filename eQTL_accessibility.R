#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

# scATAC eQTL accessibility 
# This script uses scATAC readouts to assess significant differential accessibility at disease eQTL sites, between 
# CRISPR sgRNA-perturbed cell populations and non-targeting (control) cells
ko_gene <- args[1] # pooled perturbed population to analyze
sites_n <- args[2] # consider only top n sites by coverage across KO and NT cell populations
outdir <- args[3] 

fpath = "/gpfs/commons/home/nmoss/sc_ataq/revision/human_genetics/"
fragments <- fread(paste(fpath, "snATAC_HARM_sorted_frag.bed"))
colnames(fragments) = c("seqnames", "start", "end", "barcode")
metadata <- read.csv(paste(fpath, "FRAGS_TFBS.csv"))
metadata$sgRNA <- lapply(metadata$sgRNA, gsub, pattern = "non_targeting", replacement = "non", fixed = TRUE)
metadata$gene <- gsub("\\_.*","", metadata$sgRNA)

# format fragment dataframe as GRanges object
GRangesFormat <- function(fragments_df){
  fragments_sub <- fragments_df[,c("seqnames", "start", "end")]
  colnames(fragments_sub) <- c("chr", "start", "end")
  granges <- makeGRangesFromDataFrame(fragments_sub)
  return(granges)
}

# create binary accessibility matrix - rows are eQTL sites and columns are individual cells
BinaryAccessibilty = function(ko_gene, eqtl_sites){
  if (ko_gene == "NT"){
    cells <- NT_cells
  } else {
    cells <- filter(metadata, gene == ko_gene)$DNA_BARCODE
  }
  binary_matrix <- matrix(NA, ncol = length(cells), nrow = length(eqtl_sites))
  for (i in 1:length(cells)){
    print(i)
    cell_frags <- filter(fragments, barcode == cells[i])[,1:3]
    cell_frags_granges <- GRangesFormat(cell_frags)
    overlap_cell <- suppressWarnings(countOverlaps(eqtl_sites, cell_frags_granges))
    binary_matrix[,i] <- overlap_cell
  }
  binary_matrix[binary_matrix > 1] <- 1
  sites_df <- as.data.frame(eqtl_sites)
  colnames(binary_matrix) <- cells
  rownames(binary_matrix) <- paste(sites_df$seqnames, sites_df$end, sep="_")
  return(as.data.frame(binary_matrix))
}

# read in eQTL sites and format 
eqtls <- fread("/gpfs/commons/home/nmoss/eqtl/final_liftover")
granges_eqtl <- GRangesFormat(unique(eqtls))

# select ATAC fragments 
NT_cells <- filter(metadata, gene == "non")$DNA_BARCODE
KO_cells <- filter(metadata, gene == ko_gene)$DNA_BARCODE
NT_frags <- filter(fragments, barcode %in% NT_cells)
KO_frags <- filter(fragments, barcode %in% KO_cells)

# downsample NT or KO fragment population
downsample_number <- min(nrow(KO_frags), nrow(NT_frags))
NT_sample <- sample_n(NT_frags, downsample_number)
KO_sample <- sample_n(KO_frags, downsample_number)
combined_frags <- rbind(NT_sample, KO_sample)

combined_frags_grange <- GRangesFormat(combined_frags)
overlap_eqtls <- countOverlaps(granges_eqtl, combined_frags_grange)
eqtl_ranges_df <- cbind(as.data.frame(granges_eqtl)[,c(1:3)], overlap_eqtls)

# many eQTL sites empty, select top n sites by coverage across downsampled pooled KO and NT cell populations
top_sites <- eqtl_ranges_df[order(-overlap_eqtls),][1:sites_n,]
top_sites_granges <- GRangesFormat(top_sites)
sites_df <- as.data.frame(top_sites_granges)

# calculate single-cell site-specific eQTL accessibility across KO and NT populations
binary_nt <- BinaryAccessibilty("NT", top_sites_granges)
binary_ko <- BinaryAccessibilty(ko_gene, top_sites_granges)

# calculate p-values (via proportion test) and logFC at each site
n <- c(ncol(binary_ko), ncol(binary_nt))
site_results <- as.data.frame(matrix(NA, nrow = nrow(binary_ko), ncol = 3))
for (i_region in 1:nrow(binary_ko)){
  x <- c(sum(binary_ko[i_region,]), sum(binary_nt[i_region,]))
  
  prop_res <- suppressWarnings(prop.test(x,n, alternative = "two.sided")) 
  ko_norm <- (x[1]+1)/n[1]
  nt_norm <- (x[2]+1)/n[2]
  
  site_results[i_region, 1] <- rownames(binary_ko[i_region,])
  site_results[i_region, 2] <- as.numeric(prop_res$p.value)
  site_results[i_region, 3] <- log10(abs(ko_norm/nt_norm-1))
}

colnames(site_results) <- c("eqtl_site", "pval", "logFC")
threshold_p <- site_results$pval < 0.05
site_results$significant <- threshold_p

# volcano plot: fold change vs p-value at each site, red-labeled indicate significant sites
ggplot(site_results) +
  geom_point(aes(x=logFC, y=-log10(pval), color=significant)) +
  scale_color_manual(values=c("gray", "red")) + 
  xlab(sprintf("Relative %s eQTL accessibility", ko_gene)) + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  +
  theme_classic() +
  theme(legend.position = "none")

ggsave(sprintf("%s_eQTL_top%s_volcano.pdf", ko_gene, n_sites), device = "pdf", path = outdir, width = 4, height = 3, dpi=200, useDingbats=FALSE)
write.table(site_results, paste(outdir, sprintf("%s_eQTL_top%s.txt", ko_gene, sites_n), sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")

