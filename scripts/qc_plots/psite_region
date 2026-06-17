#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

meta_file <- args[1]
gtf <- args[2]
outdir <- args[3]
offset <- as.integer(args[4])

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read_tsv(meta_file, show_col_types = FALSE)

txdb <- makeTxDbFromGFF(gtf)

cds <- unlist(cdsBy(txdb, by = "tx"))
utr5 <- unlist(fiveUTRsByTranscript(txdb))
utr3 <- unlist(threeUTRsByTranscript(txdb))

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

get_psite_region <- function(bam, sample, day, replicate, offset = 12) {

  cat("Processing:", sample, "\n")

  aln <- readGAlignments(bam)

  aln <- aln[!is.na(strand(aln)) & strand(aln) %in% c("+", "-")]

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

write_tsv(df_region, file.path(outdir, "psite_region_raw.tsv"))

df_region_sum <- df_region %>%
  group_by(replicate, day, region) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(replicate, day) %>%
  mutate(percent = n / sum(n) * 100)

write_tsv(df_region_sum, file.path(outdir, "psite_region_summary.tsv"))

df_plot <- bind_rows(df_region_sum, rna_control)

df_plot$day <- factor(
  df_plot$day,
  levels = unique(c(meta$day, "RNA"))
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
    y = "% P-sites"
  )

ggsave(
  file.path(outdir, "Fig_psite_region.png"),
  p_region,
  width = 10,
  height = 5,
  dpi = 300
)
