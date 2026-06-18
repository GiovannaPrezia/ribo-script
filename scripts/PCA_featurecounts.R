#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

count_dir <- args[1]
outdir <- args[2]
project_name <- args[3]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Load featureCounts files
# =========================

files <- list.files(
  count_dir,
  pattern = "CDS_counts.txt$",
  full.names = TRUE
)

if (length(files) < 2) {
  stop("Need at least 2 count files for PCA.")
}

count_list <- lapply(files, function(f) {
  df <- read_tsv(f, comment = "#", show_col_types = FALSE)

  sample <- basename(f) |>
    str_remove(".CDS_counts.txt$")

  counts <- df[, c("Geneid", ncol(df))]
  colnames(counts) <- c("gene_id", sample)

  counts
})

count_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), count_list)

count_matrix[is.na(count_matrix)] <- 0

gene_ids <- count_matrix$gene_id
count_matrix <- as.data.frame(count_matrix[, -1])
rownames(count_matrix) <- gene_ids

count_matrix <- round(as.matrix(count_matrix))

# =========================
# Metadata
# =========================

samples <- colnames(count_matrix)

meta <- data.frame(
  sample = samples,
  row.names = samples
)

meta$day <- str_extract(meta$sample, "DAY_?[0-9]+|D[0-9]+")
meta$day <- str_replace(meta$day, "DAY_", "D")

meta$replicate <- case_when(
  str_detect(meta$sample, "rep1|R1") ~ "R1",
  str_detect(meta$sample, "rep2|R2") ~ "R2",
  TRUE ~ "Rep"
)

# fallback para Germany/D21 se só tiver DAY21
if (all(is.na(meta$day))) {
  meta$day <- "D21"
}

# =========================
# DESeq2 VST + PCA
# =========================

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = meta,
  design = ~ 1
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

vsd <- vst(dds, blind = TRUE)

pca <- prcomp(t(assay(vsd)))

percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  day = meta$day,
  replicate = meta$replicate
)

write_tsv(
  pca_df,
  file.path(outdir, paste0(project_name, "_PCA_coordinates.tsv"))
)

# =========================
# Plot
# =========================

p <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = day, shape = replicate)
) +
  geom_point(size = 5, alpha = 0.9) +
  theme_minimal(base_size = 18) +
  labs(
    title = "PCA - Ribo-seq",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "day",
    shape = "replicate"
  ) +
  theme(
    plot.title = element_text(size = 26),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

ggsave(
  file.path(outdir, paste0(project_name, "_PCA_Riboseq.png")),
  p,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(outdir, paste0(project_name, "_PCA_Riboseq.pdf")),
  p,
  width = 8,
  height = 7
)
