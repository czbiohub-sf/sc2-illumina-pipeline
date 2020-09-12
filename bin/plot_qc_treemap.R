#!/usr/bin/Rscript

library(magrittr)
library(ggplot2)
library(ape)
library(ggtree)
library(aplot)

args <- commandArgs(trailingOnly=TRUE)

snps_csv <- args[1]
stats_tsv <- args[2]
out_svg_prefix <- args[3]

snps_df <- read.csv(snps_csv, stringsAsFactors=F)
stats_df <- read.csv(stats_tsv, sep="\t", stringsAsFactors=F)

snps_df %>%
  dplyr::select(sample, pos, gt) %>%
  tidyr::spread(pos, gt) %>%
  tibble::column_to_rownames("sample") %>%
  as.matrix() ->
  gt_mat

gt_mat_imputed <- gt_mat
for (i in 1:ncol(gt_mat_imputed)) {
  gt_mat_imputed[is.na(gt_mat_imputed[,i]), i] <- mean(
    gt_mat_imputed[,i], na.rm=TRUE)
}

gt_mat_imputed %>%
  dist() %>%
  nj() ->
  nj_tree

plot_qc_treemap <- function(keep_samples, out_svg) {
  if (length(keep_samples) <= 1) {
    # empty file
    write.table(data.frame(), out_svg, col.names=F)
    return()
  }

  gt_mat[keep_samples,] %>%
    .[,colSums(., na.rm=TRUE) > 1] %>%
    colnames() ->
    keep_pos

  if (length(keep_pos) == 0) {
    # empty file
    write.table(data.frame(), out_svg, col.names=F)
    return()
  }

  p_tree <- ggtree(keep.tip(nj_tree, keep_samples)) +
    geom_tiplab(aes(label=""), align=TRUE)

  snps_df %>%
    dplyr::filter(pos %in% keep_pos) %>%
    dplyr::filter(sample %in% keep_samples) %>%
    dplyr::mutate(alt_freq=ad1/(ad0+ad1), depth=ad0+ad1) %>%
    dplyr::mutate(pos=as.factor(pos)) %>%
    ggplot(aes(y=sample, x=pos, fill=alt_freq, alpha=depth)) +
    theme(axis.text.x=element_text(angle=90, hjust=1), axis.title.y=element_blank()) +
    scale_fill_viridis_c() +
    scale_alpha_continuous(trans="log10") +
    geom_tile() ->
    p_heatmap

  svg(out_svg,
      width=min(19, max(12, .14 * length(keep_pos))),
      height=max(8, .14 * length(keep_samples)))
  p_heatmap %>% insert_left(p_tree, width=.1) %>% print()
  dev.off()
}

plot_qc_treemap(stats_df$sample_name, paste(out_svg_prefix, "_filtered.svg", sep=""))

gt_mat %>%
  .[,colSums(., na.rm=TRUE) > 1] ->
  gt_mat_subset

if (ncol(gt_mat_subset) == 0) {
  keep_combined <- c()
} else {
  keep_combined <- rownames(gt_mat_subset)[rowMeans(is.na(gt_mat_subset)) < 1]
}

plot_qc_treemap(keep_combined,
                paste(out_svg_prefix, "_combined.svg", sep=""))
