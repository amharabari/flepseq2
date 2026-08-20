# ==============================================================================
# POLYA frequency in degradation fragments
# ==============================================================================

library(data.table)
library(stringr)
library(arrow)   # For Out-of-Memory processing
library(dplyr)   # For Arrow "remote control"
library(ggplot2)
library(viridis)
library(viridisLite)

# Define Professional Color Palette


Genotype_colors <- c("Col0"      = "#888888",
                 "xrn4"      = "#E6A817",
                 "xrn4_dne1" = "#5B9BD5")

# ==============================================================================
# 2. DATA LOADING & INITIAL FILTERING (ARROW)
# ==============================================================================
ds <- open_dataset("parquet_dataset")
SampleTable <- fread("1_Mapping_Flepseq2/barcode_correspondance.tsv")


# Identify Target Samples
target_barcodes <- SampleTable[
  Genotype %in% c("Col0", "xrn4", "xrn4_dne1") & Condition == "Seedling", Barcode]


colnames(ds)
# Memory Cleanup before loading large table

# Extract and Filter Data using Arrow Pipeline
df_targets <- ds %>%
  filter(
    Barcode %in% target_barcodes, 
    !is.na(PolyA_length_signal),
    PolyA_length_basecall > 1,
    PolyA_length_signal > 1
  ) %>%
  dplyr::select(Barcode, Feature.name, PolyA, Add_tail, PolyA_length_signal, Degradation_Tag, RA5_tag) %>%
  collect() %>% 
  setDT()

# ==============================================================================
# 3. METADATA MERGING & RENAMING
# ==============================================================================
# Clean Gene Names (Remove transcript versioning)

# Merge with Experimental Metadata
df_targets <- merge(SampleTable, df_targets, by = "Barcode", all.y = TRUE)

# Standardize Column Names
setnames(df_targets,
         old = c("Feature.name", "Rep", "PolyA_length_signal"), 
         new = c("AGI", "rep", "polya_size" ))


# ==============================================================================
# 4. QUALITY CONTROL (C50 FILTER)
# ==============================================================================
# Count reads per Gene/Rep/Genotype
grouped_counts <- df_targets[, .N, by = .(AGI, rep, Genotype)]
# Keep only genes with at least 50 reads
filtered_counts <- grouped_counts[N >= 50]

# Identify genes present in EVERY biological replicate
unique_combinations <- uniqueN(grouped_counts, by = c("Genotype", "rep"))
list_C50_df_targets <- unique(filtered_counts[, if (.N == unique_combinations) .SD, by = .(AGI)]$AGI)

print(paste("Number of High-Confidence Genes (C50):", length(list_C50_df_targets)))

# Cleanup temporary objects
rm(filtered_counts, unique_combinations)

# Create the Final High-Quality Dataset (C50)
C50 <- df_targets[AGI %in% list_C50_df_targets]
nrow(C50)
# Calculate Weights for Normalization
C50[, weight := 1/.N, by = .(AGI, rep, Genotype)]

# ==============================================================================
# 5. MATHEMATICAL PROCESSING: WEIGHTED DISTRIBUTION
# ==============================================================================
# Create full grid to include zeros for missing tail sizes
all_combinations <- CJ(
  Genotype = unique(C50$Genotype),
  rep = unique(C50$rep),
  polya_size = unique(C50$polya_size))
  
# Aggregate weighted counts
weighted_data <- C50[, .(weighted_n = sum(weight)), by = .(Genotype, polya_size, rep)]
weighted_data <- merge(all_combinations, weighted_data, by = c("Genotype", "rep", "polya_size"), all.x = TRUE)
weighted_data[is.na(weighted_n), weighted_n := 0]

# Calculate Percentages and Mean across replicates
weighted_data[, `:=` (total_weighted = sum(weighted_n), 
                      perc = 100 * (weighted_n) / sum(weighted_n)), by = .(Genotype, rep)]

weighted_data[, mean := mean(perc, na.rm = T), by = .(Genotype, polya_size)]
weighted_data <- weighted_data[!is.na(mean)]

# ==============================================================================
# 6. SMOOTHING FUNCTION (LOESS)
# ==============================================================================
precalc_smooth <- function(data, x_col, y_col, group_cols, span = 0.025) {
  smoothed_data <- data[, {
    fit <- loess(get(y_col) ~ get(x_col), 
                 data = .SD, 
                 span = span, 
                 control = loess.control(surface = "direct")) 
    smooth <- predict(fit)
    .(x = get(x_col), y = smooth)
  }, by = group_cols]
  return(smoothed_data)
}

grouping_vars <- c("Genotype")  
weighted_data_smooth <- precalc_smooth(weighted_data, "polya_size", "mean", grouping_vars)

# ==============================================================================
# 7. VISUALIZATION PREPARATION
# ==============================================================================
custom_breaks <- seq(0, 200, by = 10)

# Helper function to clean up axis labels
every_nth <- function(x, n, inverse = FALSE) {
  if (!inverse) {
    ifelse(seq_along(x) %% n == 0, x, "")
  } else {
    ifelse(seq_along(x) %% n == 0, "", x)
  }
}

major_breaks <- seq(0, 200, by = 50)
minor_grid   <- seq(0, 200, by = 10)

# ==============================================================================
# 8. GGPLOT: POLY(A) DISTRIBUTION
# ==============================================================================
weighted <- ggplot() +
  geom_point(data = weighted_data, aes(x = polya_size, y = perc, color = Genotype), size = 0.05) +
  geom_area(
    data = weighted_data_smooth,
    aes(x = x, y = y, fill = Genotype, group = Genotype),
    alpha = 0.3,
    position = "identity"
  ) +
  scale_fill_manual(values = Genotype_colors) +
  scale_color_manual(values = Genotype_colors) +
  scale_x_continuous(
    limits = c(9, 200),
    breaks = custom_breaks,
    labels = every_nth(custom_breaks, 3, inverse = TRUE)
  )+
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, by = 0.5)) +
  guides(fill = guide_legend(override.aes = list(size = 3, alpha = 0.4))) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Poly(A) tails distribution",
    x = "polyA tail size",
    y = "Frequency"
  ) +
  facet_wrap(~ Genotype, ncol = 1)
weighted
# Save Final Output
ggsave(filename = "r_plots/polyA_distribution_weighted.svg", weighted, width = 6, height = 6, dpi = 1200)




#ticks and breaks for polyA profiles
custom_breaks <- c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200)
every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}