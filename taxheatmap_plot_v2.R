#!/usr/bin/env Rscript

cat("#########################################################################################\n",
    "Welcome to TSANG script\n",
    "sangtq@ntu.edu.vn\n",
    "Nha Trang University\n",
    "Running taxonomy and Heatmap autoplot\n",
    "Version 1.1 (Italic + Clean Tax Label)\n",
    "#########################################################################################\n\n")

library(tidyverse)
library(scales)
library(readr)
library(ggplot2)
library(pheatmap)
library(vegan)
library(ggtext)
library(htmltools)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript Taxonomy_and_Heatmap_auto.R <taxonomy_file.csv> <top_n> [tax_levels] [reorder_samples] [low_abundance_threshold] [metadata_file.tsv]\n")
  cat("Example: Rscript Taxonomy_and_Heatmap_auto.R 'species_trans.csv' 30 'species' 'sort' 0.005 'metadata_16S.tsv'\n")
  quit(save = "no", status = 1)
}

# Input variables
tax_file <- args[1]
top_n <- as.integer(args[2])
tax_levels <- ifelse(length(args) >= 3, strsplit(args[3], ",")[[1]], c("phylum", "class", "order", "family", "genus", "species"))
reorder_mode <- ifelse(length(args) >= 4, args[4], "none")
low_thresh <- ifelse(length(args) >= 5, as.numeric(args[5]), 0.001)
metadata_file <- ifelse(length(args) >= 6, args[6], "metadata_16S.tsv")

# Output directory
out_dir <- "TaxPlot_Output"
dir.create(out_dir, showWarnings = FALSE)

# Color palette
palette_colors <- c("gray", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
                    "#7f7f7f", "#bcbd22", "#17becf", "#393b79", "#637939", "#8c6d31", "#843c39",
                    "#ad494a", "#e7ba52", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef", "#e6550d",
                    "#fd8d3c", "#fdae6b", "#fdd0a2", "#31a354", "#74c476", "#a1d99b", "#c7e9c0",
                    "#756bb1", "#9e9ac8", "blue")

# Read taxonomy data
tax_data <- read_csv(tax_file)
colnames(tax_data)[1] <- "index"

# Get sample names
sample_names <- if (reorder_mode == "none") {
  names(tax_data)[!grepl("index|\\.\\.\\.", names(tax_data))]
} else if (reorder_mode == "sort") {
  sort(names(tax_data)[!grepl("index|\\.\\.\\.", names(tax_data))])
} else {
  unlist(strsplit(reorder_mode, ","))
}
sample_names <- intersect(sample_names, names(tax_data))

# Function: Plot taxonomy barplot
plot_taxonomy <- function(tax_level) {
  tax_prefix <- switch(tax_level,
                       phylum = ";p__", class = ";c__", order = ";o__",
                       family = ";f__", genus = ";g__", species = ";s__",
                       stop("Unsupported taxonomic level"))

  out_name <- paste0(out_dir, "/", tax_level)
  plot_file <- paste0(out_name, ".png")
  table_file <- paste0(out_name, "_topBar.tsv")

  filtered_data <- tax_data %>%
    filter(str_detect(index, tax_prefix)) %>%
    filter(!str_detect(index, paste0(tax_prefix, "$"))) %>%
    filter(!str_detect(index, "(?i)uncultured|unassigned|__$|_bacterium|primary|Homo|Chordata")) %>%
    rowwise() %>%
    mutate(sum = sum(c_across(all_of(sample_names)), na.rm = TRUE)) %>%
    ungroup() %>%
    filter(sum / sum(sum) > low_thresh) %>%
    arrange(desc(sum))

  top_data <- slice_head(filtered_data, n = top_n)

  if (nrow(filtered_data) > top_n) {
    other_data <- slice_tail(filtered_data, n = nrow(filtered_data) - top_n) %>%
      select(all_of(sample_names)) %>%
      summarise(across(everything(), sum, na.rm = TRUE)) %>%
      mutate(index = "Others")
  } else {
    other_data <- tibble(index = "Others") %>%
       bind_cols(as_tibble(setNames(rep(0, length(sample_names)), sample_names)))
  }

  merged_data <- bind_rows(top_data, other_data) %>%
    separate(index, into = c("junk", "taxname"), sep = paste0(tax_prefix), fill = "right") %>%
    mutate(taxname = ifelse(is.na(taxname), "Others", taxname),
           taxname = str_replace_all(taxname, "^_+|_+$", ""),
           taxname = str_remove(taxname, "^s__|^g__|^f__|^o__|^c__|^p__"),
           taxname = str_replace_all(taxname, "_", " "),
           taxname = factor(taxname),
           taxname = fct_relevel(taxname, "Others", after = 0)) %>%
    pivot_longer(cols = all_of(sample_names), names_to = "Sample", values_to = "Abundance") %>%
    mutate(Sample = factor(Sample, levels = sample_names))

  label_parser <- function(labels) {
    if (tax_level %in% c("genus", "species")) {
      return(sapply(labels, function(x) {
        if (x == "Others") return("Others")
        paste0("<i>", htmlEscape(x), "</i>")
      }))
    } else {
      return(labels)
    }
  }

  p <- ggplot(merged_data, aes(x = Sample, y = Abundance, fill = taxname)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = label_percent()) +
    labs(x = "Samples", y = "Relative Abundance", fill = NULL, title = paste("Top", top_n, tax_level)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      legend.text = ggtext::element_markdown(size = 12)
    ) +
     scale_fill_manual(values = palette_colors) +
    guides(fill = guide_legend(
  label.theme = element_text(
    face = ifelse(tax_level %in% c("genus", "species"), "italic", "plain"),
    size = 12
  )
)) +
theme_classic() +
theme(
  axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold"),
  axis.text.y = element_text(size = 12, face = "bold")
)

  ggsave(plot_file, plot = p, width = 20, height = 8)
  write_tsv(merged_data, table_file)
  message(glue::glue("📊 Generating Barplot for taxonomic level: {tax_level}"))
  message(paste("✅ Saved:", plot_file, "&", table_file))
}

# Function: Generate heatmap for taxonomy level
generate_heatmap_for_level <- function(tax_level) {
  message(glue::glue("📊 Generating heatmap for taxonomic level: {tax_level}"))

  tax_prefix <- switch(tax_level,
                       phylum = ";p__", class = ";c__", order = ";o__",
                       family = ";f__", genus = ";g__", species = ";s__",
                       stop("❌ Unsupported taxonomic level"))

  tax_data <- read_csv(tax_file, show_col_types = FALSE) %>% as.data.frame()

  # Ensure first column is unique before setting row names
  if (anyDuplicated(tax_data[[1]])) {
    stop("❌ Duplicate row names detected in taxonomy data.")
  }
  
  rownames(tax_data) <- tax_data[[1]]
  tax_data <- tax_data[, -1]  # Remove the first column after setting row names

  # Ensure `sample_names` only contains valid column names
  parsed_samples <- if (reorder_mode == "none") {
    names(tax_data)[!grepl("index|\\.\\.\\.", names(tax_data))]
  } else if (reorder_mode == "sort") {
    sort(names(tax_data)[!grepl("index|\\.\\.\\.", names(tax_data))])
  } else {
    unlist(strsplit(reorder_mode, ","))
  }

  sample_names <- intersect(parsed_samples, colnames(tax_data))  # Ensure alignment
  tax_data <- tax_data[, sample_names, drop = FALSE]

  # Load metadata and ensure correct row names
  if (!file.exists(metadata_file)) stop("❌ Metadata file not found.")
  metadata <- read_tsv(metadata_file, show_col_types = FALSE)

  if (ncol(metadata) < 3) {
    stop("❌ Metadata file must contain at least three columns.")
  }

  metadata_fix <- metadata[, 2:3] %>% as.data.frame()
  rownames(metadata_fix) <- metadata[[1]]  # Set first column as row names

  # Ensure metadata rows match sample names
  metadata_fix <- metadata_fix[sample_names, , drop = FALSE]

  metadata_fix$Days <- as.character(metadata_fix$Days)
  metadata_fix$Probiotics <- as.character(metadata_fix$Probiotics)

  filtered_data <- tax_data %>%
    rownames_to_column("taxonomy") %>%
    filter(str_detect(taxonomy, tax_prefix)) %>%
    filter(!str_detect(taxonomy, "(?i)uncultured|unassigned|__$|_bacterium|primary|Homo|Chordata")) %>%
    rowwise() %>%
    mutate(total = sum(c_across(all_of(sample_names)), na.rm = TRUE)) %>%
    ungroup() %>%
    filter(total / sum(total) > low_thresh) %>%
    arrange(desc(total)) %>%
    slice_head(n = top_n) %>%
    select(-total)

  # Clean label (remove s__, g__, etc.)
  short_names <- str_remove(filtered_data$taxonomy, ".*__")
  filtered_data$short_name <- short_names
  filtered_data <- filtered_data %>% filter(!is.na(short_name))

  row_labels <- filtered_data$short_name
  if (tax_level %in% c("genus", "species")) {
    row_labels <- sapply(row_labels, function(x) parse(text = paste0("italic('", x, "')")))
  }

  heat_data <- filtered_data %>%
    select(-taxonomy) %>%
    group_by(short_name) %>%
    summarise(across(all_of(sample_names), ~ sum(.x, na.rm = TRUE))) %>%
    column_to_rownames("short_name") %>%
    as.data.frame()

  # Remove zero-sum rows and columns
  heat_data <- heat_data[rowSums(heat_data) > 0, , drop = FALSE]
  heat_data <- heat_data[, colSums(heat_data) > 0, drop = FALSE]

  write_tsv(heat_data %>% rownames_to_column("index"),
            file.path(out_dir, paste0("heatmap_", tax_level, ".tsv")))

  heat_norm <- apply(heat_data, 2, function(x) (x / sum(x)) * 100)
  bc_dissimilarity_col <- vegdist(t(heat_data), method = "bray")
  bc_dissimilarity_row <- vegdist(heat_data, method = "bray")

  pheatmap::pheatmap(as.matrix(heat_norm),
                     annotation_col = metadata_fix,
                     clustering_distance_cols = bc_dissimilarity_col,
                     clustering_distance_rows = bc_dissimilarity_row,
                     clustering_method = "average",
                     cellwidth = 15, cellheight = 12, fontsize = 8,
                    # labels_row = if (tax_level %in% c("genus", "species")) unlist(row_labels) else NULL,
                     filename = file.path(out_dir, paste0("heatmap_", tax_level, ".pdf")))

  message(glue::glue("✅ Heatmap saved to: heatmap_{tax_level}.pdf"))
}

# Run
walk(tax_levels, plot_taxonomy)
walk(tax_levels, generate_heatmap_for_level)
message("🎉 All taxonomy barplots and heatmaps generated successfully!")
