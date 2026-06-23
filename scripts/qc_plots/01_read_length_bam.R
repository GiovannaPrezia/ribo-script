#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript 01_read_length_bam.R <bam_dir> <outdir> [project_name]\n",
    "Example: Rscript 01_read_length_bam.R 05_alignment/ribo_seq/all_lengths 12_QC_Figures/all_lengths PRJNA544411"
  )
}

bam_dir <- args[1]
outdir <- args[2]
project_name <- ifelse(length(args) >= 3, args[3], "RiboSeq")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

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
  str_detect(meta$sample, "rep3|Rep3|R3") ~ "Rep3",
  TRUE ~ "Rep"
)

get_lengths <- function(bam, sample, day, replicate) {
  cat("Processing:", sample, "\n")

  aln <- readGAlignments(bam)
  read_lengths <- qwidth(aln)

  data.frame(
    length = read_lengths,
    sample = sample,
    day = day,
    replicate = replicate
  )
}

df_list <- mapply(
  get_lengths,
  meta$bam,
  meta$sample,
  meta$day,
  meta$replicate,
  SIMPLIFY = FALSE
)

df <- bind_rows(df_list)

df <- df[df$length >= 20 & df$length <= 40, ]

write.table(
  df,
  file.path(outdir, paste0(project_name, "_read_length_raw.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

df_sum <- df %>%
  group_by(sample, day, replicate, length) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample, day, replicate) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

unique_days <- sort(unique(df_sum$day))

if (length(unique_days) == 1) {
  df_sum <- df_sum %>%
    mutate(plot_label = replicate)

  plot_title <- paste0(
    "Read length distribution — ",
    project_name,
    " (",
    unique_days,
    ")"
  )

  legend_title <- "Replicate"

} else {
  df_sum <- df_sum %>%
    mutate(plot_label = paste(day, replicate))

  plot_title <- paste0(
    "Read length distribution — ",
    project_name
  )

  legend_title <- "Sample"
}

write.table(
  df_sum,
  file.path(outdir, paste0(project_name, "_read_length_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_len <- ggplot(
  df_sum,
  aes(x = length, y = percent, color = plot_label)
) +
  geom_line(linewidth = 1) +
  coord_cartesian(xlim = c(20, 40)) +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 14) +
  labs(
    title = plot_title,
    x = "Read length (nt)",
    y = "% reads",
    color = legend_title
  )

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_read_length.png")),
  p_len,
  width = 10,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_read_length.pdf")),
  p_len,
  width = 10,
  height = 5
)

cat("Read length figure saved in:", outdir, "\n")
