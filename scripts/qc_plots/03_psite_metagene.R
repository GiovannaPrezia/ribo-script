#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# =========================================================
# ARGS
# =========================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript 03_psite_metagene.R <bam_dir> <gtf> <outdir> [project_name] [offset]\n",
    "Example: Rscript 03_psite_metagene.R 05_alignment/ribo_seq/all_lengths gencode.v45.annotation.gtf 12_QC_Figures/all_lengths PRJEB29208 12"
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

if (any(is.na(meta$day))) {
  warning("Could not extract day from some sample names. Setting day = 'Unknown'.")
  meta$day[is.na(meta$day)] <- "Unknown"
}

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
# METAGENE
# =========================================================

get_metagene <- function(bam, sample, day, replicate, offset = 12) {

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

  dist_start <- ifelse(
    cds_strand == "+",
    psite_start - start(cds_hits),
    end(cds_hits) - psite_start
  )

  dist_stop <- ifelse(
    cds_strand == "+",
    psite_start - end(cds_hits),
    start(cds_hits) - psite_start
  )

  df_start <- data.frame(
    position = dist_start,
    type = "start",
    sample = sample,
    day = day,
    replicate = replicate
  )

  df_stop <- data.frame(
    position = dist_stop,
    type = "stop",
    sample = sample,
    day = day,
    replicate = replicate
  )

  rbind(df_start, df_stop)
}

df_meta_list <- mapply(
  get_metagene,
  meta$bam,
  meta$sample,
  meta$day,
  meta$replicate,
  MoreArgs = list(offset = offset),
  SIMPLIFY = FALSE
)

df_meta <- bind_rows(df_meta_list)

if (nrow(df_meta) == 0) {
  stop("No P-sites overlapped CDS regions.")
}

write.table(
  df_meta,
  file.path(outdir, paste0(project_name, "_psite_metagene_raw.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

df_meta_filtered <- df_meta %>%
  filter(position >= -50, position <= 50)

df_meta_sum <- df_meta_filtered %>%
  group_by(replicate, type, day, position) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(replicate, type, day) %>%
  mutate(freq = n / sum(n))

write.table(
  df_meta_sum,
  file.path(outdir, paste0(project_name, "_psite_metagene_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# =========================================================
# PLOT
# =========================================================

p_meta <- ggplot(
  df_meta_sum,
  aes(x = position, y = freq, color = day)
) +
  geom_line(linewidth = 0.5) +
  facet_grid(replicate ~ type) +
  scale_color_brewer(palette = "Set2") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    title = "P-site metagene profile",
    x = "Distance from start/stop codon (nt)",
    y = "Normalized P-site frequency",
    color = "day"
  )

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_psite_metagene.png")),
  p_meta,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_psite_metagene.pdf")),
  p_meta,
  width = 10,
  height = 6
)

cat("P-site metagene figure saved in:", outdir, "\n")
