#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript 06_periodicity.R <bam_dir> <gtf> <outdir> [project_name] [offset]")
}

bam_dir <- args[1]
gtf <- args[2]
outdir <- args[3]
project_name <- ifelse(length(args) >= 4, args[4], "RiboSeq")
offset <- ifelse(length(args) >= 5, as.integer(args[5]), 12)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bam_pattern <- "_Aligned.sortedByCoord.out.bam$"

bam_files <- list.files(bam_dir, pattern = bam_pattern, full.names = TRUE)

if (length(bam_files) == 0) {
  stop("No BAM files found in: ", bam_dir)
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

txdb <- makeTxDbFromGFF(gtf)
cds <- unlist(cdsBy(txdb, by = "tx"))

get_periodicity <- function(bam, sample, day, replicate, offset = 12) {
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

  if (length(hits) == 0) return(NULL)

  psite_hits <- gr_psite[queryHits(hits)]
  cds_hits <- cds[subjectHits(hits)]

  cds_strand <- as.character(strand(cds_hits))
  psite_start <- start(psite_hits)

  position_from_start <- ifelse(
    cds_strand == "+",
    psite_start - start(cds_hits),
    end(cds_hits) - psite_start
  )

  data.frame(
    sample = sample,
    day = day,
    replicate = replicate,
    position = position_from_start
  )
}

df_list <- mapply(
  get_periodicity,
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

df <- df %>%
  filter(position >= 0, position <= 150)

df_sum <- df %>%
  group_by(sample, day, replicate, position) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = n / sum(n))

write.table(
  df_sum,
  file.path(outdir, paste0(project_name, "_periodicity_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p <- ggplot(
  df_sum,
  aes(x = position, y = freq, color = day)
) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~sample, scales = "free_y") +
  theme_classic(base_size = 14) +
  labs(
    title = "3-nt periodicity near CDS start",
    x = "Position from CDS start (nt)",
    y = "Normalized P-site frequency",
    color = "day"
  )

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_periodicity.png")),
  p,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_periodicity.pdf")),
  p,
  width = 10,
  height = 6
)

cat("Periodicity figure saved in:", outdir, "\n")
