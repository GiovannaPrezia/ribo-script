#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(readr)
})

# =========================================================
# ARGS
# =========================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript 05_frame_preference.R <bam_dir> <gtf> <outdir> [project_name] [offset]\n",
    "Example: Rscript 05_frame_preference.R 05_alignment/ribo_seq/28_36 08_annotation/gencode.v45.annotation.gtf 12_QC_Figures/28_36 PRJEB29208 12"
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
cds <- unlist(cdsBy(txdb, by = "tx"))

# =========================================================
# FRAME PREFERENCE
# =========================================================

get_frame_preference <- function(bam, sample, day, replicate, offset = 12) {

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

  hits <- findOverlaps(gr_psite, cds, ignore.strand = FALSE)

  if (length(hits) == 0) {
    return(NULL)
  }

  psite_hits <- gr_psite[queryHits(hits)]
  cds_hits <- cds[subjectHits(hits)]

  cds_strand <- as.character(strand(cds_hits))
  psite_start <- start(psite_hits)

  position_from_start <- ifelse(
    cds_strand == "+",
    psite_start - start(cds_hits),
    end(cds_hits) - psite_start
  )

  frame <- position_from_start %% 3

  data.frame(
    sample = sample,
    day = day,
    replicate = replicate,
    frame = paste0("Frame ", frame)
  )
}

df_list <- mapply(
  get_frame_preference,
  meta$bam,
  meta$sample,
  meta$day,
  meta$replicate,
  MoreArgs = list(offset = offset),
  SIMPLIFY = FALSE
)

df <- bind_rows(df_list)

if (nrow(df) == 0) {
  stop("No P-sites overlapped CDS regions.")
}

df$frame <- factor(df$frame, levels = c("Frame 0", "Frame 1", "Frame 2"))

write_tsv(
  df,
  file.path(outdir, paste0(project_name, "_frame_preference_raw.tsv"))
)

df_sum <- df %>%
  group_by(sample, day, replicate, frame) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(percent = n / sum(n) * 100)

write_tsv(
  df_sum,
  file.path(outdir, paste0(project_name, "_frame_preference_summary.tsv"))
)

# =========================================================
# PLOT
# =========================================================

frame_colors <- c(
  "Frame 0" = "#7C3AED",
  "Frame 1" = "#10B981",
  "Frame 2" = "#F59E0B"
)

p <- ggplot(
  df_sum,
  aes(x = sample, y = percent, fill = frame)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = frame_colors, drop = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "P-site frame preference",
    x = "",
    y = "% P-sites in CDS",
    fill = "Frame"
  )

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_frame_preference.png")),
  p,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_frame_preference.pdf")),
  p,
  width = 10,
  height = 6
)

cat("Frame preference figure saved in:", outdir, "\n")
