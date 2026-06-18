#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript 08_contaminant_summary.R <log_dir> <outdir> <project_name>\n",
    "Example: Rscript 08_contaminant_summary.R logs 12_QC_Figures PRJEB29208"
  )
}

log_dir <- args[1]
outdir <- args[2]
project_name <- args[3]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

parse_bowtie_log <- function(file) {
  lines <- readLines(file)

  sample_mode <- basename(file) |>
    str_remove(".bowtie.log$")

  mode <- case_when(
    str_detect(sample_mode, "all_lengths") ~ "all_lengths",
    str_detect(sample_mode, "28_36") ~ "28_36",
    TRUE ~ "unknown"
  )

  sample <- sample_mode |>
    str_remove("\\.all_lengths$") |>
    str_remove("\\.28_36$")

  day <- str_extract(sample, "DAY_?[0-9]+|D[0-9]+")
  day <- str_replace(day, "DAY_", "D")
  if (is.na(day)) day <- "Unknown"

  replicate <- case_when(
    str_detect(sample, "rep1|Rep1|R1") ~ "Rep1",
    str_detect(sample, "rep2|Rep2|R2") ~ "Rep2",
    TRUE ~ "Rep"
  )

  processed <- lines[str_detect(lines, "# reads processed:")] |>
    str_extract("[0-9]+") |>
    as.numeric()

  contaminant <- lines[str_detect(lines, "# reads with at least one reported alignment:")] |>
    str_extract("[0-9]+") |>
    as.numeric()

  clean <- lines[str_detect(lines, "# reads that failed to align:")] |>
    str_extract("[0-9]+") |>
    as.numeric()

  tibble(
    sample = sample,
    mode = mode,
    day = day,
    replicate = replicate,
    processed_reads = processed,
    contaminant_reads = contaminant,
    clean_reads = clean,
    contaminant_pct = contaminant / processed * 100,
    retention_pct = clean / processed * 100
  )
}

bowtie_logs <- list.files(
  log_dir,
  pattern = ".bowtie.log$",
  full.names = TRUE
)

if (length(bowtie_logs) == 0) {
  stop("No Bowtie logs found in: ", log_dir)
}

df <- bind_rows(lapply(bowtie_logs, parse_bowtie_log))

write_tsv(
  df,
  file.path(outdir, paste0(project_name, "_contaminant_summary.tsv"))
)

df_long <- df %>%
  select(sample, mode, day, contaminant_reads, clean_reads) %>%
  pivot_longer(
    cols = c(contaminant_reads, clean_reads),
    names_to = "class",
    values_to = "reads"
  ) %>%
  mutate(
    class = recode(
      class,
      contaminant_reads = "Removed contaminants",
      clean_reads = "Retained clean reads"
    )
  )

p_counts <- ggplot(
  df_long,
  aes(x = sample, y = reads / 1e6, fill = class)
) +
  geom_col(width = 0.75) +
  facet_wrap(~mode, scales = "free_x") +
  scale_fill_manual(values = c(
    "Removed contaminants" = "#DC2626",
    "Retained clean reads" = "#10B981"
  )) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Bowtie contaminant filtering",
    x = "",
    y = "Reads (millions)",
    fill = ""
  )

p_retention <- ggplot(
  df,
  aes(x = sample, y = retention_pct, fill = day)
) +
  geom_col(width = 0.75) +
  facet_wrap(~mode, scales = "free_x") +
  scale_fill_manual(values = c(
    "D8" = "#2563EB",
    "D10" = "#7C3AED",
    "D18" = "#DC2626",
    "D21" = "#059669",
    "Unknown" = "gray40"
  ), drop = FALSE) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Retention after contaminant filtering",
    x = "",
    y = "% retained reads",
    fill = "day"
  )

dashboard <- p_counts / p_retention +
  plot_annotation(
    title = paste0("Contaminant QC — ", project_name),
    theme = theme(plot.title = element_text(size = 22, face = "bold"))
  )

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_contaminant_summary.png")),
  dashboard,
  width = 14,
  height = 10,
  dpi = 300
)

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_contaminant_summary.pdf")),
  dashboard,
  width = 14,
  height = 10
)

cat("Contaminant summary saved in:", outdir, "\n")
