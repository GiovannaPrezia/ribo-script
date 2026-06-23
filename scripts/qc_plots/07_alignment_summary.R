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
    "Usage: Rscript 07_alignment_summary.R <star_qc_dir> <outdir> <project_name>\n",
    "Example: Rscript 07_alignment_summary.R 06_star_qc/ribo_seq 12_QC_Figures PRJEB29208"
  )
}

star_qc_dir <- args[1]
outdir <- args[2]
project_name <- args[3]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# PALETTE
# =========================================================

metric_colors <- c(
  "Unique" = "#7C3AED",
  "Multimapping" = "#F59E0B",
  "Too_many_loci" = "#DC2626",
  "Too_short" = "#10B981",
  "Other_unmapped" = "gray70"
)

# =========================================================
# FUNCTIONS
# =========================================================

extract_star_value <- function(lines, pattern) {
  line <- lines[str_detect(lines, fixed(pattern))]
  if (length(line) == 0) return(NA_character_)
  value <- str_split(line[1], "\\|", simplify = TRUE)[, 2]
  str_trim(value)
}

parse_star_log <- function(file) {
  lines <- readLines(file)

  sample_mode <- basename(file) |>
    str_remove("_Log.final.out$")

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

  tibble(
    sample = sample,
    mode = mode,
    day = day,
    replicate = replicate,
    input_reads = as.numeric(extract_star_value(lines, "Number of input reads")),
    unique_reads = as.numeric(extract_star_value(lines, "Uniquely mapped reads number")),
    unique_pct = as.numeric(str_remove(extract_star_value(lines, "Uniquely mapped reads %"), "%")),
    multimapping_reads = as.numeric(extract_star_value(lines, "Number of reads mapped to multiple loci")),
    multimapping_pct = as.numeric(str_remove(extract_star_value(lines, "% of reads mapped to multiple loci"), "%")),
    too_many_loci_reads = as.numeric(extract_star_value(lines, "Number of reads mapped to too many loci")),
    too_many_loci_pct = as.numeric(str_remove(extract_star_value(lines, "% of reads mapped to too many loci"), "%")),
    too_short_reads = as.numeric(extract_star_value(lines, "Number of reads unmapped: too short")),
    too_short_pct = as.numeric(str_remove(extract_star_value(lines, "% of reads unmapped: too short"), "%")),
    other_unmapped_reads = as.numeric(extract_star_value(lines, "Number of reads unmapped: other")),
    other_unmapped_pct = as.numeric(str_remove(extract_star_value(lines, "% of reads unmapped: other"), "%")),
    avg_input_length = as.numeric(extract_star_value(lines, "Average input read length")),
    avg_mapped_length = as.numeric(extract_star_value(lines, "Average mapped length"))
  )
}

# =========================================================
# LOAD STAR LOGS
# =========================================================

star_logs <- list.files(
  star_qc_dir,
  pattern = "_Log.final.out$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(star_logs) == 0) {
  stop("No STAR Log.final.out files found in: ", star_qc_dir)
}

df <- bind_rows(lapply(star_logs, parse_star_log))

write_tsv(
  df,
  file.path(outdir, paste0(project_name, "_alignment_summary.tsv"))
)

# =========================================================
# LONG TABLE FOR PLOTS
# =========================================================

df_long <- df %>%
  select(
    sample, mode, day, replicate,
    Unique = unique_pct,
    Multimapping = multimapping_pct,
    Too_many_loci = too_many_loci_pct,
    Too_short = too_short_pct,
    Other_unmapped = other_unmapped_pct
  ) %>%
  pivot_longer(
    cols = c(Unique, Multimapping, Too_many_loci, Too_short, Other_unmapped),
    names_to = "metric",
    values_to = "percent"
  )

# =========================================================
# PLOT 1 — UNIQUE MAPPING
# =========================================================

p_unique <- ggplot(
  df,
  aes(x = sample, y = unique_pct, fill = day)
) +
  geom_col(width = 0.75) +
  facet_wrap(~mode, scales = "free_x") +
  scale_fill_viridis_d() +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "STAR unique mapping",
    x = "",
    y = "% uniquely mapped reads",
    fill = "day"
  )

# =========================================================
# PLOT 2 — ALIGNMENT COMPOSITION
# =========================================================

p_composition <- ggplot(
  df_long,
  aes(x = sample, y = percent, fill = metric)
) +
  geom_col(width = 0.75) +
  facet_wrap(~mode, scales = "free_x") +
  scale_fill_manual(values = metric_colors) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "STAR alignment composition",
    x = "",
    y = "% input reads",
    fill = "Mapping class"
  )

# =========================================================
# PLOT 3 — READ DEPTH
# =========================================================

df_depth <- df %>%
  select(sample, mode, day, input_reads, unique_reads) %>%
  pivot_longer(
    cols = c(input_reads, unique_reads),
    names_to = "class",
    values_to = "reads"
  ) %>%
  mutate(class = recode(
    class,
    input_reads = "Input STAR",
    unique_reads = "Unique mapped"
  ))

p_depth <- ggplot(
  df_depth,
  aes(x = sample, y = reads / 1e6, fill = class)
) +
  geom_col(position = "dodge", width = 0.75) +
  facet_wrap(~mode, scales = "free_x") +
  scale_fill_manual(values = c(
    "Input STAR" = "gray70",
    "Unique mapped" = "#7C3AED"
  )) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "STAR read depth",
    x = "",
    y = "Reads (millions)",
    fill = ""
  )

# =========================================================
# DASHBOARD
# =========================================================

dashboard <- (p_unique / p_composition / p_depth) +
  plot_annotation(
    title = paste0("Alignment QC — ", project_name),
    theme = theme(
      plot.title = element_text(size = 22, face = "bold")
    )
  )

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_alignment_summary.png")),
  dashboard,
  width = 14,
  height = 14,
  dpi = 300
)

ggsave(
  file.path(outdir, paste0(project_name, "_Fig_alignment_summary.pdf")),
  dashboard,
  width = 14,
  height = 14
)

cat("Alignment summary saved in:", outdir, "\n")
