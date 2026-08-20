# ==============================================================================
# 21/04/2026 - FLEPSEQ FRAGMENTATION ANALYSIS
# ==============================================================================
setwd("~/AMHLABO/")

library(data.table)
library(arrow)     # For Out-of-Memory processing
library(dplyr)     # For Arrow "remote control"
library(pbapply)   # For progress bars
library(ggplot2)

setDTthreads(0)    # Use all 12 threads

# --- 1. SETTINGS, THEME & COLORS ---
parquet_dir   <- "flepseq2/parquet_dataset"
metadata_path <- "flepseq2/bamparsed/barcode_correspondance.tsv"
SampleTable <- fread("flepseq2/bamparsed/barcode_correspondance.tsv")
ds <- open_dataset(parquet_dir)

# Factors and Colors
actual_tags <- c("Intact", "5'_truncated", "3'_truncated", "Both_truncated")
tag_colors  <- c("Intact"="#666666", "5'_truncated"="#D95F02", "3'_truncated"="#1F78B4", "Both_truncated"="#009E73")
genotype_colors <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9")

# Updated Theme with Grids
bio_journal_theme <- theme_minimal(base_size = 12) +
  theme(
    # Remove the distracting tiny minor lines
    panel.grid.minor = element_blank(),
    
    # Add light gray lines for both X and Y axes
    panel.grid.major.x = element_line(color = "black", linetype = "dotted", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "black", linetype = "dotted", linewidth = 0.3),
    
    # Keep the rest of your clean styling
    axis.line = element_line(color = "black", linewidth = 0.5),
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "top",   
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
  )

# Get total reads per barcode from your full dataset
library_sizes <- ds %>%
  group_by(Barcode) %>%
  summarize(libsize = n()) %>%
  collect() %>% setDT()

# Merge this into your sample table
SampleTable <- merge(SampleTable, library_sizes, by = "Barcode", all.x = TRUE)

# --- 2. COVERAGE FUNCTION ---
calculate_coverage_dt <- function(dt, bin_width = 10) { #10
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt_c <- dt[!is.na(Read.start) & !is.na(Read.end)]
  # binning
  dt_c[, b_left  := floor(Read.start / bin_width) * bin_width]
  dt_c[, b_right := floor(Read.end / bin_width) * bin_width]
  # indexing rows
  dt_c[, row_id := .I]
  
  coverage <- dt_c[, .(
    pos = seq(b_left, b_right, by = bin_width)
  ), by = .(Barcode, Genotype, Feature.name, Feature.strand, Degradation_Tag, row_id)] #% add Rep

# .N returns the coverage depth (in number of reads) at the position pos (can be 110, 120 etc) grouped by Metadata 
  return(coverage[, .(depth = .N), by = .(Barcode, Genotype, Feature.name, Feature.strand, Degradation_Tag, pos)]) #% add Rep
}

# --- 3. DATA INITIALIZATION ---

dne1_ids <- c("AT1G13245", "AT1G78080", "AT2G18160", "AT2G22430", "AT2G47400", 
              "AT3G22380", "AT4G34138", "AT5G19120", "AT4G64260", "AT2G43020")
dne1_targets_full <- paste0(dne1_ids, ".1")

# Lists of differentiolly uridilated in xrn4_dne1 vs xrn4; 
# Possibility of visualization of profiles of fragmentts in different genotypes and see if there any differnces in distribution of polyA tails
MUG <- fread("flepseq2/provided/134_genes_MUG.tsv", header=TRUE, sep="\t")
LUG <- fread("flepseq2/provided/33_genes_LUG.tsv", header=TRUE, sep="\t")

trunc3_Stat <- fread("flepseq2/diff_frag_helene_trunc3_Stat.txt", header=TRUE, sep="\t")

#-------------------------------------------------------------------------------------------------------
#  Col0-xrn4 ANALYSIS 5' truncated  ##########################################################################################################
Col0_xrn4_trunc5 <- fread("flepseq2/fragment_analysis/diff_frag_5_truncated_Col0_xrn4_Stat.txt",header=TRUE, sep="\t")
target_barcodes <- SampleTable[Genotype %in% c("Col0", "xrn4") & Condition == "Seedling", Barcode] #Col0

# PART A: PROPORTION ANALYSIS (Example Gene)
df_prop <- ds %>%
  filter(Feature.name %in% Col0_xrn4_trunc5$Feature.name,
         Barcode %in% target_barcodes)%>%
  dplyr::select(Barcode, Feature.name, Degradation_Tag) %>%
  collect() %>% setDT()

df_prop <- merge(df_prop, SampleTable, by = "Barcode", all.x = TRUE)
df_prop[, Degradation_Tag := factor(trimws(Degradation_Tag), levels = actual_tags)]

prop_summary <- df_prop[, .(n = .N), by = .(Feature.name, Genotype, Rep, Degradation_Tag)]
prop_summary[, perc := 100 * n / sum(n), by = .(Feature.name, Genotype, Rep)]
prop_mean    <- prop_summary[, .(mean = mean(perc)), by = .(Feature.name, Genotype, Degradation_Tag)]

p1 <- ggplot(prop_summary, aes(x = Degradation_Tag, fill = Degradation_Tag)) + 
  geom_bar(data = prop_mean, aes(y = mean), stat = "identity", alpha = 0.7, position = position_dodge(width = 0.9)) +
  geom_point(aes(y = perc, group = Genotype), 
             position = position_dodge(width = 0.9), size = 1.2, shape = 21, fill = "white") +
  
  facet_grid(Genotype ~ Feature.name) +
  scale_fill_manual(values =tag_colors) +
  scale_x_discrete(drop = FALSE) +
  bio_journal_theme + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + # Vertical text for narrow columns
  labs(title = "Proportions: Differentially fragmented Col0-xrn4 5' truncated ", y = "Reads (%)")

print(p1)

ggsave("flepseq2/plots/Prop_Col0-xrn4_5p.pdf", p1, width = 12, height = 5)

# PART B: COVERAGE ANALYSIS (DNE1 TARGETS)
df_targets <- ds %>%
  filter(Feature.name %in% Col0_xrn4_trunc5$Feature.name,           #sig_rows$groupID[6:10], #MUG$Feature.name[73:100],#dne1_targets_full, ################################## sig_rows$groupID[0:11]
         Barcode %in% target_barcodes) %>%
  dplyr::select(Barcode, Feature.name, Feature.strand, Read.start, Read.end, Degradation_Tag) %>%
  collect() %>% setDT()

df_targets <- df_targets[!is.na(Feature.name)] #%
# Merge metadata and fix factors
df_targets <- merge(unique(df_targets), 
                    SampleTable[, .(Barcode, Genotype, Rep, libsize)], #% Add Rep
                    by = "Barcode", all.x = TRUE)
df_targets[, Degradation_Tag := factor(trimws(Degradation_Tag), levels = actual_tags)]


coverage_data <- calculate_coverage_dt(df_targets, bin_width = 10)
# Check the depth magnitude (Should be much less than 1 million now)
print(paste("Maximum Raw Depth:", max(coverage_data$depth, na.rm = TRUE)))

# Attach the Library Sizes (libsize) from our SampleTable
coverage_norm <- merge(coverage_data, 
                       SampleTable[, .(Barcode, libsize)], 
                       by = "Barcode", all.x = TRUE)

# Calculate CPM (Counts Per Million) Formula: (Reads at pos / Total reads in library) * 1,000,000
coverage_norm[, cpm := (depth / libsize) * 1e6]


final_coverage <- coverage_norm[, .(mean_cpm = mean(cpm)), 
                                by = .(Genotype, Feature.name, Feature.strand, Degradation_Tag, pos)] #% Add Rep

p2 <- ggplot(final_coverage, aes(x = pos, y = mean_cpm, color = Degradation_Tag)) +
  geom_line(linewidth = 1) +
  # Use "free" to allow both X and Y scales to adapt to each gene
  facet_grid(Genotype ~ Feature.name + Feature.strand, scales = "free_x") +  #% Add + Rep
  scale_color_manual(values = tag_colors, drop = FALSE) + 
  bio_journal_theme +
  # Use pretty_breaks so each unique scale gets appropriate labels
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), 
                     expand = expansion(mult = c(0, 0.1))) + 
  labs(title = "Coverage: Differentially fragmented Col0-xrn4 5' truncated", 
       x = "Chromosome Coordinate", 
       y = "Normalized Coverage (CPM)")
print(p2)

ggsave("flepseq2/plots/Coverage_Col_xnr4_5p.pdf", p2, width = 12, height = 5)

library(patchwork)
# Combine plots: / means top/bottom, | means side-by-side
combined_plot <- (p1 | p2) + 
  plot_annotation(title = "DNE1 Target Analysis Summary",
                  tag_levels = 'A') # Automatically adds labels A, B, etc.

#ggsave("flepseq2/Combined_DNE1_Analysis.pdf", combined_plot, width = 15, height = 20)


#  Col0-xrn4 ANALYSIS 3' truncated  ##########################################################################################################
Col0_xrn4_trunc3 <- fread("flepseq2/fragment_analysis/diff_frag_3_truncated_Col0_xrn4_Stat.txt",header=TRUE, sep="\t")
target_barcodes <- SampleTable[Genotype %in% c("Col0", "xrn4") & Condition == "Seedling", Barcode] #Col0

# PART A: PROPORTION ANALYSIS (Example Gene)
df_prop <- ds %>%
  filter(Feature.name %in% Col0_xrn4_trunc3$Feature.name[0:12],
         Barcode %in% target_barcodes)%>%
  dplyr::select(Barcode, Feature.name, Degradation_Tag) %>%
  collect() %>% setDT()

df_prop <- merge(df_prop, SampleTable, by = "Barcode", all.x = TRUE)
df_prop[, Degradation_Tag := factor(trimws(Degradation_Tag), levels = actual_tags)]

prop_summary <- df_prop[, .(n = .N), by = .(Feature.name, Genotype, Rep, Degradation_Tag)]
prop_summary[, perc := 100 * n / sum(n), by = .(Feature.name, Genotype, Rep)]
prop_mean    <- prop_summary[, .(mean = mean(perc)), by = .(Feature.name, Genotype, Degradation_Tag)]

p1 <- ggplot(prop_summary, aes(x = Degradation_Tag, fill = Degradation_Tag)) + 
  geom_bar(data = prop_mean, aes(y = mean), stat = "identity", alpha = 0.7, position = position_dodge(width = 0.9)) +
  geom_point(aes(y = perc, group = Genotype), 
             position = position_dodge(width = 0.9), size = 1.2, shape = 21, fill = "white") +
  
  facet_grid(Genotype ~ Feature.name) +
  scale_fill_manual(values =tag_colors) +
  scale_x_discrete(drop = FALSE) +
  bio_journal_theme + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + # Vertical text for narrow columns
  labs(title = "Proportions: Differentially fragmented Col0-xrn4 3' truncated ", y = "Reads (%)")

print(p1)

ggsave("flepseq2/plots/Prop_Col0-xrn4_3p.pdf", p1, width = 12, height = 5)

# PART B: COVERAGE ANALYSIS (DNE1 TARGETS)
df_targets <- ds %>%
  filter(Feature.name %in% Col0_xrn4_trunc3$Feature.name[0:12],           #sig_rows$groupID[6:10], #MUG$Feature.name[73:100],#dne1_targets_full, ################################## sig_rows$groupID[0:11]
         Barcode %in% target_barcodes) %>%
  dplyr::select(Barcode, Feature.name, Feature.strand, Read.start, Read.end, Degradation_Tag) %>%
  collect() %>% setDT()

df_targets <- df_targets[!is.na(Feature.name)] #%
# Merge metadata and fix factors
df_targets <- merge(unique(df_targets), 
                    SampleTable[, .(Barcode, Genotype, Rep, libsize)], #% Add Rep
                    by = "Barcode", all.x = TRUE)
df_targets[, Degradation_Tag := factor(trimws(Degradation_Tag), levels = actual_tags)]


coverage_data <- calculate_coverage_dt(df_targets, bin_width = 10)
# Check the depth magnitude (Should be much less than 1 million now)
print(paste("Maximum Raw Depth:", max(coverage_data$depth, na.rm = TRUE)))

# Attach the Library Sizes (libsize) from our SampleTable
coverage_norm <- merge(coverage_data, 
                       SampleTable[, .(Barcode, libsize)], 
                       by = "Barcode", all.x = TRUE)

# Calculate CPM (Counts Per Million) Formula: (Reads at pos / Total reads in library) * 1,000,000
coverage_norm[, cpm := (depth / libsize) * 1e6]


final_coverage <- coverage_norm[, .(mean_cpm = mean(cpm)), 
                                by = .(Genotype, Feature.name, Feature.strand, Degradation_Tag, pos)] #% Add Rep

p2 <- ggplot(final_coverage, aes(x = pos, y = mean_cpm, color = Degradation_Tag)) +
  geom_line(linewidth = 1) +
  # Use "free" to allow both X and Y scales to adapt to each gene
  facet_grid(Genotype ~ Feature.name + Feature.strand, scales = "free_x") +  #% Add + Rep
  scale_color_manual(values = tag_colors, drop = FALSE) + 
  bio_journal_theme +
  # Use pretty_breaks so each unique scale gets appropriate labels
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), 
                     expand = expansion(mult = c(0, 0.1))) + 
  labs(title = "Coverage: Differentially fragmented Col0-xrn4 3' truncated", 
       x = "Chromosome Coordinate", 
       y = "Normalized Coverage (CPM)")
print(p2)

ggsave("flepseq2/plots/Coverage_Col_xnr4_3p.pdf", p2, width = 12, height = 5)

library(patchwork)
# Combine plots: / means top/bottom, | means side-by-side
combined_plot <- (p1 | p2) + 
  plot_annotation(title = "DNE1 Target Analysis Summary",
                  tag_levels = 'A') # Automatically adds labels A, B, etc.

#ggsave("flepseq2/Combined_DNE1_Analysis.pdf", combined_plot, width = 15, height = 20)

#  xrn4-xrn4_dne1 ANALYSIS 3' truncated  ##########################################################################################################
xrn4_xrn4_dne1_trunc3 <- fread("flepseq2/fragment_analysis/diff_frag_3_truncated_xrn4_xrn4_dne1_Stat.txt",header=TRUE, sep="\t")
target_barcodes <- SampleTable[Genotype %in% c("xrn4", "xrn4_dne1") & Condition == "Seedling", Barcode] #Col0

# PART A: PROPORTION ANALYSIS (Example Gene)
df_prop <- ds %>%
  filter(Feature.name %in% xrn4_xrn4_dne1_trunc3$Feature.name[13:28],
         Barcode %in% target_barcodes)%>%
  dplyr::select(Barcode, Feature.name, Degradation_Tag) %>%
  collect() %>% setDT()

df_prop <- merge(df_prop, SampleTable, by = "Barcode", all.x = TRUE)
df_prop[, Degradation_Tag := factor(trimws(Degradation_Tag), levels = actual_tags)]

prop_summary <- df_prop[, .(n = .N), by = .(Feature.name, Genotype, Rep, Degradation_Tag)]
prop_summary[, perc := 100 * n / sum(n), by = .(Feature.name, Genotype, Rep)]
prop_mean    <- prop_summary[, .(mean = mean(perc)), by = .(Feature.name, Genotype, Degradation_Tag)]

p1 <- ggplot(prop_summary, aes(x = Degradation_Tag, fill = Degradation_Tag)) + 
  geom_bar(data = prop_mean, aes(y = mean), stat = "identity", alpha = 0.7, position = position_dodge(width = 0.9)) +
  geom_point(aes(y = perc, group = Genotype), 
             position = position_dodge(width = 0.9), size = 1.2, shape = 21, fill = "white") +
  
  facet_grid(Genotype ~ Feature.name) +
  scale_fill_manual(values =tag_colors) +
  scale_x_discrete(drop = FALSE) +
  bio_journal_theme + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + # Vertical text for narrow columns
  labs(title = "Proportions: Differentially fragmented xrn4-xrn4_dne1 3' truncated ", y = "Reads (%)")

print(p1)

ggsave("flepseq2/plots/Prop_xrn4-xrn4_dne1_3p_13-28.pdf", p1, width = 12, height = 5)

# PART B: COVERAGE ANALYSIS (DNE1 TARGETS)
df_targets <- ds %>%
  filter(Feature.name %in% xrn4_xrn4_dne1_trunc3$Feature.name[13:28],           #sig_rows$groupID[6:10], #MUG$Feature.name[73:100],#dne1_targets_full, ################################## sig_rows$groupID[0:11]
         Barcode %in% target_barcodes) %>%
  dplyr::select(Barcode, Feature.name, Feature.strand, Read.start, Read.end, Degradation_Tag) %>%
  collect() %>% setDT()

df_targets <- df_targets[!is.na(Feature.name)] #%
# Merge metadata and fix factors
df_targets <- merge(unique(df_targets), 
                    SampleTable[, .(Barcode, Genotype, Rep, libsize)], #% Add Rep
                    by = "Barcode", all.x = TRUE)
df_targets[, Degradation_Tag := factor(trimws(Degradation_Tag), levels = actual_tags)]


coverage_data <- calculate_coverage_dt(df_targets, bin_width = 10)
# Check the depth magnitude (Should be much less than 1 million now)
print(paste("Maximum Raw Depth:", max(coverage_data$depth, na.rm = TRUE)))

# Attach the Library Sizes (libsize) from our SampleTable
coverage_norm <- merge(coverage_data, 
                       SampleTable[, .(Barcode, libsize)], 
                       by = "Barcode", all.x = TRUE)

# Calculate CPM (Counts Per Million) Formula: (Reads at pos / Total reads in library) * 1,000,000
coverage_norm[, cpm := (depth / libsize) * 1e6]


final_coverage <- coverage_norm[, .(mean_cpm = mean(cpm)), 
                                by = .(Genotype, Feature.name, Feature.strand, Degradation_Tag, pos)] #% Add Rep

p2 <- ggplot(final_coverage, aes(x = pos, y = mean_cpm, color = Degradation_Tag)) +
  geom_line(linewidth = 1) +
  # Use "free" to allow both X and Y scales to adapt to each gene
  facet_grid(Genotype ~ Feature.name + Feature.strand, scales = "free_x") +  #% Add + Rep
  scale_color_manual(values = tag_colors, drop = FALSE) + 
  bio_journal_theme +
  # Use pretty_breaks so each unique scale gets appropriate labels
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), 
                     expand = expansion(mult = c(0, 0.1))) + 
  labs(title = "Coverage: Differentially fragmented xrn4-xrn4_dne1 3' truncated", 
       x = "Chromosome Coordinate", 
       y = "Normalized Coverage (CPM)")
print(p2)

ggsave("flepseq2/plots/Coverage_xrn4_xrn4_dne1_3p_13-28.pdf", p2, width = 12, height = 5)

library(patchwork)
# Combine plots: / means top/bottom, | means side-by-side
combined_plot <- (p1 | p2) + 
  plot_annotation(title = "DNE1 Target Analysis Summary",
                  tag_levels = 'A') # Automatically adds labels A, B, etc.

#ggsave("flepseq2/Combined_DNE1_Analysis.pdf", combined_plot, width = 15, height = 20)



#  xrn4-xrn4_dne1 ANALYSIS 5' truncated  ##########################################################################################################
xrn4_xrn4_dne1_trunc5 <- fread("flepseq2/fragment_analysis/diff_frag_5_truncated_xrn4_xrn4_dne1_Stat.txt",header=TRUE, sep="\t")
target_barcodes <- SampleTable[Genotype %in% c("xrn4", "xrn4_dne1") & Condition == "Seedling", Barcode] #Col0

# PART A: PROPORTION ANALYSIS (Example Gene)
df_prop <- ds %>%
  filter(Feature.name %in% xrn4_xrn4_dne1_trunc5$Feature.name,
         Barcode %in% target_barcodes)%>%
  dplyr::select(Barcode, Feature.name, Degradation_Tag) %>%
  collect() %>% setDT()

df_prop <- merge(df_prop, SampleTable, by = "Barcode", all.x = TRUE)
df_prop[, Degradation_Tag := factor(trimws(Degradation_Tag), levels = actual_tags)]

prop_summary <- df_prop[, .(n = .N), by = .(Feature.name, Genotype, Rep, Degradation_Tag)]
prop_summary[, perc := 100 * n / sum(n), by = .(Feature.name, Genotype, Rep)]
prop_mean    <- prop_summary[, .(mean = mean(perc)), by = .(Feature.name, Genotype, Degradation_Tag)]

p1 <- ggplot(prop_summary, aes(x = Degradation_Tag, fill = Degradation_Tag)) + 
  geom_bar(data = prop_mean, aes(y = mean), stat = "identity", alpha = 0.7, position = position_dodge(width = 0.9)) +
  geom_point(aes(y = perc, group = Genotype), 
             position = position_dodge(width = 0.9), size = 1.2, shape = 21, fill = "white") +
  
  facet_grid(Genotype ~ Feature.name) +
  scale_fill_manual(values =tag_colors) +
  scale_x_discrete(drop = FALSE) +
  bio_journal_theme + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + # Vertical text for narrow columns
  labs(title = "Proportions: Differentially fragmented xrn4-xrn4_dne1 5' truncated ", y = "Reads (%)")

print(p1)

ggsave("flepseq2/plots/Prop_xrn4-xrn4_dne1_5p.pdf", p1, width = 12, height = 5)

# PART B: COVERAGE ANALYSIS (DNE1 TARGETS)
df_targets <- ds %>%
  filter(Feature.name %in% xrn4_xrn4_dne1_trunc5$Feature.name,           #sig_rows$groupID[6:10], #MUG$Feature.name[73:100],#dne1_targets_full, ################################## sig_rows$groupID[0:11]
         Barcode %in% target_barcodes) %>%
  dplyr::select(Barcode, Feature.name, Feature.strand, Read.start, Read.end, Degradation_Tag) %>%
  collect() %>% setDT()

df_targets <- df_targets[!is.na(Feature.name)] #%
# Merge metadata and fix factors
df_targets <- merge(unique(df_targets), 
                    SampleTable[, .(Barcode, Genotype, Rep, libsize)], #% Add Rep
                    by = "Barcode", all.x = TRUE)
df_targets[, Degradation_Tag := factor(trimws(Degradation_Tag), levels = actual_tags)]


coverage_data <- calculate_coverage_dt(df_targets, bin_width = 10)
# Check the depth magnitude (Should be much less than 1 million now)
print(paste("Maximum Raw Depth:", max(coverage_data$depth, na.rm = TRUE)))

# Attach the Library Sizes (libsize) from our SampleTable
coverage_norm <- merge(coverage_data, 
                       SampleTable[, .(Barcode, libsize)], 
                       by = "Barcode", all.x = TRUE)

# Calculate CPM (Counts Per Million) Formula: (Reads at pos / Total reads in library) * 1,000,000
coverage_norm[, cpm := (depth / libsize) * 1e6]


final_coverage <- coverage_norm[, .(mean_cpm = mean(cpm)), 
                                by = .(Genotype, Feature.name, Feature.strand, Degradation_Tag, pos)] #% Add Rep

p2 <- ggplot(final_coverage, aes(x = pos, y = mean_cpm, color = Degradation_Tag)) +
  geom_line(linewidth = 1) +
  # Use "free" to allow both X and Y scales to adapt to each gene
  facet_grid(Genotype ~ Feature.name + Feature.strand, scales = "free_x") +  #% Add + Rep
  scale_color_manual(values = tag_colors, drop = FALSE) + 
  bio_journal_theme +
  # Use pretty_breaks so each unique scale gets appropriate labels
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), 
                     expand = expansion(mult = c(0, 0.1))) + 
  labs(title = "Coverage: Differentially fragmented xrn4-xrn4_dne1 5' truncated", 
       x = "Chromosome Coordinate", 
       y = "Normalized Coverage (CPM)")
print(p2)

ggsave("flepseq2/plots/Coverage_xrn4_xrn4_dne1_5p.pdf", p2, width = 12, height = 5)

library(patchwork)
# Combine plots: / means top/bottom, | means side-by-side
combined_plot <- (p1 | p2) + 
  plot_annotation(title = "DNE1 Target Analysis Summary",
                  tag_levels = 'A') # Automatically adds labels A, B, etc.

#ggsave("flepseq2/Combined_DNE1_Analysis.pdf", combined_plot, width = 15, height = 20)


#  xrn4-xrn4_dne1 ANALYSIS BOTH truncated  ##########################################################################################################
xrn4_xrn4_dne1_both <- fread("flepseq2/fragment_analysis/diff_frag_Both_truncated_xrn4_xrn4_dne1_Stat.txt",header=TRUE, sep="\t")
target_barcodes <- SampleTable[Genotype %in% c("xrn4", "xrn4_dne1") & Condition == "Seedling", Barcode] #Col0

# PART A: PROPORTION ANALYSIS (Example Gene)
df_prop <- ds %>%
  filter(Feature.name %in% xrn4_xrn4_dne1_both$Feature.name,
         Barcode %in% target_barcodes)%>%
  dplyr::select(Barcode, Feature.name, Degradation_Tag) %>%
  collect() %>% setDT()

df_prop <- merge(df_prop, SampleTable, by = "Barcode", all.x = TRUE)
df_prop[, Degradation_Tag := factor(trimws(Degradation_Tag), levels = actual_tags)]

prop_summary <- df_prop[, .(n = .N), by = .(Feature.name, Genotype, Rep, Degradation_Tag)]
prop_summary[, perc := 100 * n / sum(n), by = .(Feature.name, Genotype, Rep)]
prop_mean    <- prop_summary[, .(mean = mean(perc)), by = .(Feature.name, Genotype, Degradation_Tag)]

p1 <- ggplot(prop_summary, aes(x = Degradation_Tag, fill = Degradation_Tag)) + 
  geom_bar(data = prop_mean, aes(y = mean), stat = "identity", alpha = 0.7, position = position_dodge(width = 0.9)) +
  geom_point(aes(y = perc, group = Genotype), 
             position = position_dodge(width = 0.9), size = 1.2, shape = 21, fill = "white") +
  
  facet_grid(Genotype ~ Feature.name) +
  scale_fill_manual(values =tag_colors) +
  scale_x_discrete(drop = FALSE) +
  bio_journal_theme + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + # Vertical text for narrow columns
  labs(title = "Proportions: Differentially fragmented xrn4-xrn4_dne1 Both truncated ", y = "Reads (%)")

print(p1)

ggsave("flepseq2/plots/Prop_xrn4-xrn4_dne1_both.pdf", p1, width = 12, height = 12)

# PART B: COVERAGE ANALYSIS (DNE1 TARGETS)
df_targets <- ds %>%
  filter(Feature.name %in% xrn4_xrn4_dne1_both$Feature.name,           #sig_rows$groupID[6:10], #MUG$Feature.name[73:100],#dne1_targets_full, ################################## sig_rows$groupID[0:11]
         Barcode %in% target_barcodes) %>%
  dplyr::select(Barcode, Feature.name, Feature.strand, Read.start, Read.end, Degradation_Tag) %>%
  collect() %>% setDT()

df_targets <- df_targets[!is.na(Feature.name)] #%
# Merge metadata and fix factors
df_targets <- merge(unique(df_targets), 
                    SampleTable[, .(Barcode, Genotype, Rep, libsize)], #% Add Rep
                    by = "Barcode", all.x = TRUE)
df_targets[, Degradation_Tag := factor(trimws(Degradation_Tag), levels = actual_tags)]


coverage_data <- calculate_coverage_dt(df_targets, bin_width = 10)
# Check the depth magnitude (Should be much less than 1 million now)
print(paste("Maximum Raw Depth:", max(coverage_data$depth, na.rm = TRUE)))

# Attach the Library Sizes (libsize) from our SampleTable
coverage_norm <- merge(coverage_data, 
                       SampleTable[, .(Barcode, libsize)], 
                       by = "Barcode", all.x = TRUE)

# Calculate CPM (Counts Per Million) Formula: (Reads at pos / Total reads in library) * 1,000,000
coverage_norm[, cpm := (depth / libsize) * 1e6]


final_coverage <- coverage_norm[, .(mean_cpm = mean(cpm)), 
                                by = .(Genotype, Feature.name, Feature.strand, Degradation_Tag, pos)] #% Add Rep

p2 <- ggplot(final_coverage, aes(x = pos, y = mean_cpm, color = Degradation_Tag)) +
  geom_line(linewidth = 1) +
  # Use "free" to allow both X and Y scales to adapt to each gene
  facet_grid(Genotype ~ Feature.name + Feature.strand, scales = "free_x") +  #% Add + Rep
  scale_color_manual(values = tag_colors, drop = FALSE) + 
  bio_journal_theme +
  # Use pretty_breaks so each unique scale gets appropriate labels
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), 
                     expand = expansion(mult = c(0, 0.1))) + 
  labs(title = "Coverage: Differentially fragmented xrn4-xrn4_dne1 Both truncated", 
       x = "Chromosome Coordinate", 
       y = "Normalized Coverage (CPM)")
print(p2)

ggsave("flepseq2/plots/Coverage_xrn4_xrn4_dne1_both.pdf", p2, width = 7, height = 5)

library(patchwork)
# Combine plots: / means top/bottom, | means side-by-side
combined_plot <- (p1 | p2) + 
  plot_annotation(title = "DNE1 Target Analysis Summary",
                  tag_levels = 'A') # Automatically adds labels A, B, etc.

#ggsave("flepseq2/Combined_DNE1_Analysis.pdf", combined_plot, width = 15, height = 20)




