#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(writexl)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 09_qc_master_table.R <project_root> <project_name>")
}

project_root <- args[1]
project_name <- args[2]

raw_dir <- file.path(project_root, "02_fastq/ribo_seq")
trim_dir <- file.path(project_root, "03_trimmed/ribo_seq")
clean_dir <- file.path(project_root, "04_cleaned/ribo_seq")
star_qc_dir <- file.path(project_root, "06_star_qc/ribo_seq")
ribotricer_dir <- file.path(project_root, "10_Ribotricer")
log_dir <- file.path(project_root, "logs")
outdir <- file.path(project_root, "13_Report/tables")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(outdir, paste0(project_name, "_QC_master_table.tsv"))
xlsx_outfile <- file.path(outdir, paste0(project_name, "_QC_master_table.xlsx"))
csv_outfile <- file.path(outdir, paste0(project_name, "_QC_master_table.csv"))

get_fastq_stats <- function(fastq) {
  if (!file.exists(fastq)) return(tibble(Reads = NA, Length = NA, GC = NA))

  x <- read_tsv(pipe(paste("seqkit stats -T", shQuote(fastq))), show_col_types = FALSE)

  tibble(
    Reads = x$num_seqs[1],
    Length = x$avg_len[1],
    GC = if ("%GC" %in% colnames(x)) x$`%GC`[1] else NA
  )
}

count_fastq_reads <- function(fastq) {
  if (!file.exists(fastq)) return(NA_real_)
  as.numeric(system2("bash", c("-c", shQuote(paste0("zcat ", fastq, " | wc -l"))), stdout = TRUE)) / 4
}

extract_star_value <- function(file, pattern) {
  if (!file.exists(file)) return(NA_character_)
  lines <- readLines(file, warn = FALSE)
  line <- lines[str_detect(lines, fixed(pattern))]
  if (length(line) == 0) return(NA_character_)
  str_trim(str_split(line[1], "\\|", simplify = TRUE)[, 2])
}

parse_bowtie <- function(file) {
  if (!file.exists(file)) {
    return(tibble(Post_Contaminant_Reads = NA, Retention_Percent = NA))
  }

  lines <- readLines(file, warn = FALSE)

  processed <- lines[str_detect(lines, "# reads processed:")] |>
    str_extract("[0-9]+") |>
    as.numeric()

  clean <- lines[str_detect(lines, "# reads that failed to align:")] |>
    str_extract("[0-9]+") |>
    as.numeric()

  tibble(
    Post_Contaminant_Reads = clean,
    Retention_Percent = round(clean / processed * 100, 2)
  )
}

parse_cutadapt_removed <- function(file) {
  if (!file.exists(file)) return(NA_real_)

  lines <- readLines(file, warn = FALSE)

  total_line <- lines[str_detect(lines, "Total reads processed:")]
  written_line <- lines[str_detect(lines, "Reads written \\(passing filters\\):")]

  if (length(total_line) == 0 || length(written_line) == 0) return(NA_real_)

  total_reads <- total_line[1] |>
    str_extract("[0-9,]+$") |>
    str_replace_all(",", "") |>
    as.numeric()

  written_reads <- written_line[1] |>
    str_extract("[0-9,]+") |>
    str_replace_all(",", "") |>
    as.numeric()

  round((1 - written_reads / total_reads) * 100, 2)
}

count_rows_file <- function(file) {
  if (!file.exists(file)) return(NA_integer_)
  nrow(read_tsv(file, show_col_types = FALSE))
}

get_ribotricer_counts <- function(ribotricer_dir, sample, mode) {
  mode_dir <- file.path(ribotricer_dir, mode)

  if (!dir.exists(mode_dir)) {
    return(tibble(
      ORFs_detected = NA,
      smORFs_detected = NA,
      lncRNA_smORFs = NA,
      High_confidence_smORFs = NA
    ))
  }

  orf_file <- file.path(mode_dir, paste0(sample, ".", mode, "_translating_ORFs.tsv"))
  smorf_file <- file.path(mode_dir, paste0(sample, ".", mode, "_smorfs_20_150aa.tsv"))
  lncrna_file <- file.path(mode_dir, paste0(sample, ".", mode, "_lncrna_smorfs.tsv"))
  hc_file <- file.path(mode_dir, paste0(mode, "_high_confidence_candidates.tsv"))

  hc_count <- NA_integer_

  if (file.exists(hc_file)) {
    hc <- read_tsv(hc_file, show_col_types = FALSE)

    if ("samples" %in% colnames(hc)) {
      hc_count <- sum(str_detect(hc$samples, fixed(paste0(sample, ".", mode))))
    } else {
      hc_count <- nrow(hc)
    }
  }

  tibble(
    ORFs_detected = count_rows_file(orf_file),
    smORFs_detected = count_rows_file(smorf_file),
    lncRNA_smORFs = count_rows_file(lncrna_file),
    High_confidence_smORFs = hc_count
  )
}

raw_fastqs <- list.files(raw_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)

if (length(raw_fastqs) == 0) {
  stop("No raw FASTQ files found in: ", raw_dir)
}

rows <- list()

for (raw_fastq in raw_fastqs) {
  sample <- basename(raw_fastq) |> str_remove("\\.fastq\\.gz$")

  raw_stats <- get_fastq_stats(raw_fastq)

  trimmed_fastq <- file.path(trim_dir, paste0(sample, ".trim.fastq.gz"))
  trimmed_stats <- get_fastq_stats(trimmed_fastq)

  cutadapt_log <- file.path(log_dir, paste0(sample, ".cutadapt.log"))
  percent_removed <- parse_cutadapt_removed(cutadapt_log)

  post_n_polyg_fastq <- file.path(
    clean_dir,
    paste0(sample, ".trim.noN.noPolyG.all_lengths.fastq.gz")
  )

  post_n_polyg_reads <- count_fastq_reads(post_n_polyg_fastq)

  for (mode in c("all_lengths", "28_36")) {
    clean_fastq <- file.path(
      clean_dir,
      mode,
      paste0(sample, ".", mode, ".clean.fastq.gz")
    )

    bowtie_log <- file.path(log_dir, paste0(sample, ".", mode, ".bowtie.log"))
    star_log <- file.path(star_qc_dir, paste0(sample, ".", mode, "_Log.final.out"))

    clean_stats <- get_fastq_stats(clean_fastq)
    bowtie_stats <- parse_bowtie(bowtie_log)
    ribo_stats <- get_ribotricer_counts(ribotricer_dir, sample, mode)

    rows[[length(rows) + 1]] <- tibble(
      Sample = sample,
      Run = NA_character_,
      Mode = mode,

      Raw_Reads = raw_stats$Reads,
      Raw_Length = raw_stats$Length,
      Raw_GC = raw_stats$GC,

      Trimmed_Reads = trimmed_stats$Reads,
      Percent_Removed = percent_removed,
      Trimmed_Length = trimmed_stats$Length,
      Trimmed_GC = trimmed_stats$GC,
      Trim_Info = "Cutadapt_-u3_polyA_A10_min17",

      Post_N_PolyG_Reads = post_n_polyg_reads,
      Post_Contaminant_Reads = bowtie_stats$Post_Contaminant_Reads,
      Retention_Percent = bowtie_stats$Retention_Percent,

      Clean_GC = clean_stats$GC,
      Clean_Length = clean_stats$Length,
      Clean_Info = "noN_noPolyG_Bowtie_contaminant_filter",

      Input_STAR = extract_star_value(star_log, "Number of input reads"),
      Unique_Reads_STAR = extract_star_value(star_log, "Uniquely mapped reads number"),
      Unique_Mapping_Percent = extract_star_value(star_log, "Uniquely mapped reads %"),

      ORFs_detected = ribo_stats$ORFs_detected,
      smORFs_detected = ribo_stats$smORFs_detected,
      lncRNA_smORFs = ribo_stats$lncRNA_smORFs,
      High_confidence_smORFs = ribo_stats$High_confidence_smORFs
    )
  }
}

qc_table <- bind_rows(rows)

write_tsv(qc_table, outfile)
write_csv(qc_table, csv_outfile)
write_xlsx(qc_table, xlsx_outfile)

cat("QC master table saved:\n")
cat(outfile, "\n")
cat(csv_outfile, "\n")
cat(xlsx_outfile, "\n")
