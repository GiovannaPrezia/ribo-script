#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

# =========================================================
# ARGS
# =========================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript 02_psite_region.R <bam_dir> <gtf> <outdir> [project_name] [offset]\n",
    "Example: Rscript 02_psite_region.R 05_alignment/ribo_seq/all_lengths 08_annotation/gencode.v45.annotation.gtf 12_QC_Figures/all_lengths PRJEB29208 12"
  )
}

bam_dir <- args[1]
gtf <- args[2]
outdir <- args[3]
project_name <- ifelse(length(args) >= 4, args[4], "RiboSeq")
offset <- ifelse(length(args) >= 5, as.integer(args[5]), 12)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# BAM DISCOVERY
# =========================================================

bam_pattern <- "_Aligned.sortedByCoord.out.bam$"

bam_files <- list.files(
  bam_dir,
  pattern = bam_pattern,
  full.names = TRUE
)

if (length(bam_files) == 0) {
  stop("No BAM files were found in: ", bam_dir)
}

sample_names <- basename(bam_files) |>
  str_remove(bam_pattern)

meta <- data.frame(
  sample = sample_names,
  bam = bam_files,
  stringsAsFactors = FALSE
)

meta$day <- str_extract(meta$sample, "DAY_?[0-9]+|D[0-9]+")
meta$day <- str_replace(meta$day, "DAY_", "D")
meta$day[is.na(meta$day)] <- "Unknown"

meta$replicate <- case_when(
  str_detect(meta$sample, "rep1|Rep1|R1") ~ "Rep1",
  str_detect(meta$sample, "rep2|Rep2|R2") ~ "Rep2",
  TRUE ~ "Rep"
)

# =========================================================
# ANNOTATION
# =========================================================

txdb <- makeTxDbFromGFF(gtf)

cds  <- unlist(cdsBy(txdb, by = "tx"))
utr5 <- unlist(fiveUTRsByTranscript(txdb))
utr3 <- unlist(threeUTRsByTranscript(txdb))

# =========================================================
# RNA CONTROL
# =========================================================

cds_len  <- sum(width(cds))
utr5_len <- sum(width(utr5))
utr3_len <- sum(width(utr3))

total <- cds_len + utr5_len + utr3_len

rna_control <- data.frame(
  region = c("5UTR", "CDS", "3UTR"),
  percent = c(
    utr5_len / total * 100,
    cds_len  / total * 100,
    utr3_len / total * 100
  )
) %>%
  mutate(day = "RNA") %>%
  tidyr::crossing(replicate = unique(meta$replicate))

# =========================================================
# P-SITE REGION
# =========================================================

get_psite_region <- function(bam, sample, day, replicate, offset = 12) {

  cat("Processing:", sample, "\n")

  aln <- readGAlignments(bam)
  aln <- aln[as.character(strand(aln)) %in% c("+", "-")]

  psite_pos <- ifelse(
    as.character(strand(aln)) == "+",
    start(aln) + offset,
    end(aln) - offset
  )

  gr_psite <- GRanges(
    seqnames = seqnames(aln),
    ranges = IRanges(psite_pos, width = 1),
    strand = strand(aln)
  )

  region <- rep("Other", length(gr_psite))

  region[overlapsAny(gr_psite, utr5, ignore.strand = FALSE)] <- "5UTR"
  region[overlapsAny(gr_psite, utr3, ignore.strand = FALSE)] <- "3UTR"
  region[overlapsAny(gr_psite, cds,  ignore.strand = FALSE)] <- "CDS"

  data.frame(
    region = region,
    sample = sample,
    day = day,
    replicate = replicate
  )
}

df_region_list <- mapply(
  get_psite_region,
  meta$bam,
  meta$sample,
  meta$day,
  meta$replicate,
  MoreArgs = list(offset = offset),
  SIMPLIFY = FALSE
)

df_region <- bind_rows(df_region_list)

write_tsv(
  df_region,
  file.path(outdir, paste0(project_name, "_psite_region_raw.tsv"))
)

# =========================================================
# SUMMARY
# =========================================================

df_region_sum <- df_region %>%
  group_by(replicate, day, region) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(replicate, day) %>%
  mutate(percent = n / sum(n) * 100)

write_tsv(
  df_region_sum,
  file.path(outdir, paste0(project_name, "_psite_region_summary.tsv"))
)

# =========================================================
# PLOT
# =========================================================

df_plot <- bind_rows(df_region_sum, rna_control)

df_plot$day <- factor(
  df_plot$day,
  levels = unique(c(sort(unique(meta$day)), "RNA"))
)

p_region <- ggplot(
  df_plot,
  aes(x = day, y = percent, fill = region)
) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  facet_wrap(~replicate) +
  scale_fill_manual(values = c(
    "5UTR" = "#10B981",
    "CDS" = "#7C3AED",
    "3UTR" = "#F59E0B",
    "Other" = "gray80"
  )) +
  theme_classic(base_size = 14) +
  labs(
    title = "P-site regional distribution",
    x = "",
    y = "% P-sites",
    fill = "Region"
  )

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_psite_region.png")),
  p_region,
  width = 10,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_psite_region.pdf")),
  p_region,
  width = 10,
  height = 5
)

cat("P-site region figure saved in:", outdir, "\n")
