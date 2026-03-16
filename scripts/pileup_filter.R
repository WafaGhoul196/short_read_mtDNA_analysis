# =============================================================================
#  pileup_filter.R
#  Filter and transform pileup heteroplasmy CSV files.
#  Accepts a results directory as a command-line argument.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ── Parse argument ────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "."
cat(sprintf("Working directory: %s\n", results_dir))

# ── Variant transformation function ──────────────────────────────────────────
transform_to_csv_format <- function(df) {
  # Handle the '*' unmapped column (R reads it as 'X.' or literal '*')
  if ("X." %in% colnames(df)) {
    df <- df %>% rename(unmapped = X.)
  } else if ("*" %in% colnames(df)) {
    df <- df %>% rename(unmapped = `*`)
  } else {
    df$unmapped <- 0
  }

  df <- df %>%
    mutate(total_nucleotides = A + T + G + C + ins + del + unmapped)

  # ── Insertions ──────────────────────────────────────────────────────────────
  ins_data <- df %>%
    select(Pos, RefCall,
           ins_type1, ins_count1, ins_type2, ins_count2, ins_type3, ins_count3,
           Depth, Heteroplasmy, total_nucleotides, unmapped) %>%
    pivot_longer(cols = starts_with("ins_type"),
                 names_to = "ins_type", values_to = "ins_val") %>%
    mutate(ins_count = case_when(
      ins_type == "ins_type1" ~ ins_count1,
      ins_type == "ins_type2" ~ ins_count2,
      ins_type == "ins_type3" ~ ins_count3
    )) %>%
    filter(ins_val != "") %>%
    mutate(
      ins_val  = toupper(ins_val),
      Type     = "insertion",
      AF       = ins_count / total_nucleotides,
      Allele_Count = ins_count
    ) %>%
    group_by(Pos, RefCall, ins_val, Depth, Heteroplasmy, Type, unmapped) %>%
    summarise(AF = sum(AF), Allele_Count = sum(Allele_Count), .groups = "drop") %>%
    rename(ALT = ins_val)

  # ── Deletions ───────────────────────────────────────────────────────────────
  del_data <- df %>%
    select(Pos, RefCall,
           del_type1, del_count1, del_type2, del_count2, del_type3, del_count3,
           Depth, Heteroplasmy, total_nucleotides, unmapped) %>%
    pivot_longer(cols = starts_with("del_type"),
                 names_to = "del_type", values_to = "del_val") %>%
    mutate(del_count = case_when(
      del_type == "del_type1" ~ del_count1,
      del_type == "del_type2" ~ del_count2,
      del_type == "del_type3" ~ del_count3
    )) %>%
    filter(del_val != "") %>%
    mutate(
      del_val  = toupper(del_val),
      Type     = "deletion",
      AF       = del_count / total_nucleotides,
      Allele_Count = del_count
    ) %>%
    group_by(Pos, RefCall, del_val, Depth, Heteroplasmy, Type, unmapped) %>%
    summarise(AF = sum(AF), Allele_Count = sum(Allele_Count), .groups = "drop") %>%
    rename(ALT = del_val)

  # ── SNVs ────────────────────────────────────────────────────────────────────
  snv_data <- df %>%
    select(Pos, RefCall, A, T, G, C, Depth, Heteroplasmy, total_nucleotides, unmapped) %>%
    pivot_longer(cols = c(A, T, G, C), names_to = "Variant", values_to = "Count") %>%
    filter(!is.na(Count) & Count > 0 & Variant != RefCall) %>%
    mutate(
      Type     = "SNV",
      ALT      = Variant,
      AF       = Count / total_nucleotides,
      Allele_Count = Count
    ) %>%
    select(Pos, RefCall, ALT, AF, Depth, Heteroplasmy, Type, Allele_Count, unmapped)

  bind_rows(ins_data, del_data, snv_data) %>%
    arrange(Pos) %>%
    select(Pos, RefCall, ALT, AF, Depth, Heteroplasmy, Type, Allele_Count, unmapped)
}

# ── Process files ─────────────────────────────────────────────────────────────
pattern <- file.path(results_dir, "output_pileup_analysis_.*_heteroplasmy\\.csv$")
pileup_files <- list.files(path = results_dir,
                           pattern = "^output_pileup_analysis_.*_heteroplasmy\\.csv$",
                           full.names = TRUE)

if (length(pileup_files) == 0) {
  cat(sprintf("No heteroplasmy CSV files found in: %s\n", results_dir))
  quit(status = 1)
}

cat(sprintf("Found %d file(s) to process.\n", length(pileup_files)))

all_filtered_data <- data.frame()

for (file in pileup_files) {
  patient_name <- basename(file) %>%
    gsub("^output_pileup_analysis_", "", .) %>%
    gsub("_heteroplasmy\\.csv$", "", .)

  cat(sprintf("\nProcessing: %s\n", patient_name))
  data <- tryCatch(read.csv(file), error = function(e) {
    cat(sprintf("  ERROR reading file: %s\n", e$message)); NULL
  })
  if (is.null(data)) next

  num_nonzero <- sum(data$Heteroplasmy != 0, na.rm = TRUE)
  cat(sprintf("  Positions with heteroplasmy ≠ 0 : %d\n", num_nonzero))

  data_final    <- transform_to_csv_format(data)

  data_filtered <- subset(data_final,
    ((Heteroplasmy > 0.05 & AF > 0.05) | AF > 0.08) &
    (unmapped / Depth < 0.7))

  cat(sprintf("  Positions after filtering       : %d\n", nrow(data_filtered)))

  out_file <- file.path(results_dir,
                        paste0(patient_name, "_data_filtred_heteroplasmy.csv"))
  write.csv(data_filtered, out_file, row.names = FALSE)
  cat(sprintf("  Saved → %s\n", basename(out_file)))

  data_filtered$Patient <- patient_name
  all_filtered_data <- bind_rows(all_filtered_data, data_filtered)
}

# ── Save merged output ────────────────────────────────────────────────────────
merged_out <- file.path(results_dir, "all_patients_filtered_data.csv")
write.csv(all_filtered_data, merged_out, row.names = FALSE)
cat(sprintf("\nAll patients merged → %s\n", merged_out))

# ── Visualisations ────────────────────────────────────────────────────────────
plots_dir <- file.path(results_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Bar plot: raw vs filtered variant counts
filter_stats <- data.frame()
for (file in pileup_files) {
  patient_name <- basename(file) %>%
    gsub("^output_pileup_analysis_", "", .) %>%
    gsub("_heteroplasmy\\.csv$", "", .)
  raw_data         <- read.csv(file)
  num_raw          <- sum(raw_data$Heteroplasmy != 0, na.rm = TRUE)
  filtered_n       <- nrow(subset(all_filtered_data, Patient == patient_name))
  filter_stats     <- rbind(filter_stats, data.frame(
    Patient    = patient_name,
    Raw        = num_raw,
    Filtered   = filtered_n,
    FilterRate = ifelse(num_raw > 0, filtered_n / num_raw * 100, 0)
  ))
}

p1 <- ggplot(filter_stats, aes(x = Patient)) +
  geom_bar(aes(y = Raw,      fill = "Raw"),      stat = "identity", alpha = 0.7) +
  geom_bar(aes(y = Filtered, fill = "Filtered"), stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("Raw" = "grey60", "Filtered" = "#2166AC")) +
  labs(title = "Variant counts before and after filtering",
       y = "Number of variants", x = "Patient", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plots_dir, "filter_stats_barplot.png"), p1,
       width = 10, height = 6, dpi = 150)

# 2. Heatmap: heteroplasmy by genomic position
if (nrow(all_filtered_data) > 0) {
  all_filtered_data <- all_filtered_data %>%
    mutate(bin = floor(Pos / 100) * 100)

  p2 <- ggplot(all_filtered_data, aes(x = bin, y = Patient, fill = Heteroplasmy)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "plasma") +
    labs(title = "Heteroplasmy heatmap by genomic position",
         x = "Position (100 bp bins)", y = "Patient") +
    theme_minimal()

  ggsave(file.path(plots_dir, "heteroplasmy_heatmap.png"), p2,
         width = 14, height = max(4, nrow(filter_stats) * 0.5 + 2), dpi = 150)
}

cat(sprintf("\nPlots saved to: %s\n", plots_dir))
cat("Filtering complete.\n")
