# ==============================================================================
# POLYA frequency in degradation fragments
# ==============================================================================
setwd("~/AMHLABO/")

library(data.table)
library(stringr)
library(arrow)   # For Out-of-Memory processing
library(dplyr)   # For Arrow "remote control"
library(ggplot2)
library(viridis)
library(viridisLite)

# Optimization for data.table
setDTthreads(0)

# Define Professional Color Palette
Genotype_colors <- c(
  "Col0"      = "#696969", 
  "xrn4"      = "#FF1493", 
  "xrn4_dne1" = "#00BFFF"
)

# ==============================================================================
# 2. DATA LOADING & INITIAL FILTERING (ARROW)
# ==============================================================================
parquet_dir <- "flepseq2/parquet_dataset"
ds <- open_dataset(parquet_dir)
SampleTable <- fread("flepseq2/bamparsed/barcode_correspondance.tsv")

# Identify Target Samples
target_barcodes <- SampleTable[
  Genotype %in% c("Col0", "xrn4", "xrn4_dne1") & Condition == "Seedling", 
  Barcode]

RA5_tags   <- c("NO_RA5", "RA5")
degradation_tags   <- c("Intact", "5'_truncated", "3'_truncated", "Both_truncated")


# Memory Cleanup before loading large table
gc()

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
df_targets[, Feature.name := sub("\\..*", "", Feature.name)]

# Merge with Experimental Metadata
df_targets <- merge(SampleTable, df_targets, by = "Barcode", all.y = TRUE)

df_targets <- df_targets[
  Degradation_Tag %in% c("Intact", "5'_truncated", "3'_truncated", "Both_truncated") | RA5_tag == "RA5"]

# Standardize Column Names
setnames(df_targets,
         old = c("Feature.name", "Rep", "PolyA_length_signal"), 
         new = c("AGI", "rep", "polya_size" ))

# Build a unified category variable for faceting
df_targets[, read_class := fcase(
  RA5_tag == "RA5",              "RA5",
  Degradation_Tag == "Intact",   "Intact",
  Degradation_Tag == "5'_truncated", "5'_truncated",
  Degradation_Tag == "3'_truncated", "3'_truncated",
  Degradation_Tag == "Both_truncated", "Both_truncated",
  default = NA_character_
)]
df_targets <- df_targets[!is.na(read_class)]
# ==============================================================================
# 4. QUALITY CONTROL (C50 FILTER)
# ==============================================================================
# Count reads per Gene/Rep/Genotype
grouped_counts <- df_targets[, .N, by = .(AGI, rep, Genotype, read_class)]

# Keep only genes with at least 50 reads
filtered_counts <- grouped_counts[N >= 50]

# Identify genes present in EVERY biological replicate
unique_combinations <- uniqueN(grouped_counts, by = c("Genotype", "rep"))
list_C50_df_targets <- unique(filtered_counts[, if (.N == unique_combinations) .SD, by = .(AGI, read_class)]$AGI)

print(paste("Number of High-Confidence Genes (C50):", length(list_C50_df_targets)))

# Cleanup temporary objects
rm(filtered_counts, unique_combinations)

# Create the Final High-Quality Dataset (C50)
C50 <- df_targets[AGI %in% list_C50_df_targets]
nrow(C50)
# Calculate Weights for Normalization
C50[, weight := 1/.N, by = .(AGI, rep, Genotype, read_class)]

# ==============================================================================
# 5. MATHEMATICAL PROCESSING: WEIGHTED DISTRIBUTION
# ==============================================================================
# Create full grid to include zeros for missing tail sizes
all_combinations <- CJ(
  Genotype = unique(C50$Genotype),
  rep = unique(C50$rep),
  polya_size = unique(C50$polya_size),
  read_class = unique(C50$read_class))

# Aggregate weighted counts
weighted_data <- C50[, .(weighted_n = sum(weight)), by = .(Genotype, polya_size, rep, read_class)]
weighted_data <- merge(all_combinations, weighted_data, by = c("Genotype", "rep", "polya_size", "read_class"), all.x = TRUE)
weighted_data[is.na(weighted_n), weighted_n := 0]

# Calculate Percentages and Mean across replicates
weighted_data[, `:=` (total_weighted = sum(weighted_n), 
                      perc = 100 * (weighted_n) / sum(weighted_n)), by = .(Genotype, rep, read_class)]

weighted_data[, mean := mean(perc, na.rm = T), by = .(Genotype, polya_size, read_class)]
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

grouping_vars <- c("Genotype", "read_class")  
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
  # Raw data points
  geom_point(data = weighted_data, aes(x = polya_size, y = perc, color = Genotype), size = 0.05) +
  # Smoothed distribution area
  geom_area(data = weighted_data_smooth, aes(x = x, y = y, fill = Genotype, 
                                             group = interaction(Genotype, read_class)), 
            alpha = 0.3, position = "identity") +  
  
  # Color Mapping
  scale_fill_manual(values = Genotype_colors) +
  scale_color_manual(values = Genotype_colors) +
  
  # Axis Scaling
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 200),
                     breaks = custom_breaks,        
                     minor_breaks = minor_grid, 
                     labels = every_nth(custom_breaks, 100, inverse = TRUE)) +
  
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1)) +
  
  # Legend and Styling
  guides(fill = guide_legend(override.aes = list(size = 3, alpha = 0.4))) +
  theme_bw(base_size = 12) + 
  theme(legend.position  = "top",
        strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold")) +
  labs(x = "polyA tail size", y = "frequency") +
  facet_grid(read_class ~ Genotype)   # <- rows = read_class, cols = Genotype

# Display Overlay Plot
weighted

# Display Faceted Plot
weighted + facet_wrap(~Genotype)

# Save Final Output
ggsave(filename = "flepseq2/fragment_analysis/polyA_distribution_across_tags.pdf", weighted, width = 11, height = 9, dpi = 1200)


# peak functions
# detect the first peak
find_first_peak <- function(x, y) {
  # Find local maxima
  peaks <- which(diff(sign(diff(y))) == -2) + 1
  
  # If no peaks are found, return NA
  if (length(peaks) == 0) {
    return(list(x = NA_real_, y = NA_real_))
  }
  
  # Select the first peak based on its position in the density (smallest x value)
  first_peak <- peaks[1]
  
  # Return the x and y value of the first peak
  return(list(x = x[first_peak], y = y[first_peak]))
}

# detect the most prominent peak in polyA size densities
find_single_peak <- function(x, y) {
  # Identify peaks
  peaks <- which(diff(sign(diff(y))) == -2) + 1
  if (length(peaks) == 0) {
    return(list(peak.x = NA_real_))
  }
  
  # Calculate prominences for each peak
  prominences <- y[peaks] - sapply(peaks, function(peak) {
    min(c(min(y[1:peak]), min(y[(peak + 1):length(y)])))
  })
  
  # Return x position of the most prominent peak as named list
  return(list(peak_x = x[peaks[which.max(prominences)]]))
}

# Function to calculate smoothed values
precalc_smooth <- function(data, x_col, y_col, group_cols, span = span) {
  smoothed_data <- C50[, {
    fit <- loess(get(y_col) ~ get(x_col), 
                 data= .SD, 
                 span = span, 
                 control = loess.control(surface = "direct"))  # Add this!
    smooth <- predict(fit)
    .(x = get(x_col), y = smooth)
  }, by = group_cols]
  return(smoothed_data)
}

# detect the most prom peak between 35 and 60
find_second_peak <- function(x, y) {
  # Restrict to the desired range
  in_range <- which(x >= 35 & x <= 60)
  if (length(in_range) < 3) {
    return(list(peak_x = NA_real_))
  }
  
  x_sub <- x[in_range]
  y_sub <- y[in_range]
  
  # Identify peaks: local maxima with a slope change
  peaks <- which(diff(sign(diff(y_sub))) == -2) + 1
  if (length(peaks) == 0) {
    return(list(peak_2 = NA_real_))
  }
  
  # Calculate prominences for each peak
  prominences <- sapply(peaks, function(peak) {
    left_min <- min(y_sub[1:peak])
    right_min <- min(y_sub[(peak + 1):length(y_sub)])
    y_sub[peak] - min(c(left_min, right_min))
  })
  
  # Find the most prominent peak
  most_prominent_index <- which.max(prominences)
  peak_2_value <- x_sub[peaks[most_prominent_index]]
  
  return(list(peak_2 =peak_2_value))
}


## peaks per replicate r density by Agi and rep
DT_density <- C50[, .(density_data = list(density(polya_size, bw = "sj", adjust = 1))), by = .(AGI, Genotype, rep)]

# Extract x and y from density_data list column
density_table <- DT_density[, {
  x <- density_data[[1]]$x
  y <- density_data[[1]]$y
  .(x, y)
}, by = .(AGI, Genotype, rep)]

# calculate density and retrieve the results  to a table
all_1_pooled_peaks<- density_table[, {maxima <- find_first_peak(x, y)}, by = .(AGI, Genotype, rep)]
all_1_pooled_peaks[, peak_1:= round(x, 0)]
all_1_pooled_peaks[,x := NULL]

all_p_pooled_peaks <- density_table[, {maxima <- find_single_peak(x, y)}, by = .(AGI, Genotype, rep)]
all_p_pooled_peaks[, peak_p := round(peak_x, 0)]
all_p_pooled_peaks[, peak_x := NULL]

Table_peaks <- merge(all_1_pooled_peaks, all_p_pooled_peaks, by = c("AGI", "Genotype","rep"))
Table_peaks <- Table_peaks[, .(AGI, Genotype, rep, peak_1, peak_p)]
rm(all_1_pooled_peaks, all_p_pooled_peaks, DT_density, density_table)

## profile classification
Table_peaks[, profile_1 := fcase(
  peak_p < 35, "short",
  peak_p >= 35 & peak_p <= 60, "inter",
  peak_p > 60, "long"
)]
Table_peaks[, profile_1 := factor(profile_1, levels  = c("short", "inter", "long"))]


#################################################################################

