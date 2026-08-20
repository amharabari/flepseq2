



setwd("~/AMHLABO/")

library(data.table)
library(arrow)     # For Out-of-Memory processing
library(dplyr)     # For Arrow "remote control"
library(pbapply)   # For progress bars
library(ggplot2)
library(stringr)

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
  ), by = .(Barcode, Genotype, Rep, Condition, Sample_name, Feature.name, Feature.strand, Degradation_Tag, row_id)] #% add Rep
  
  # .N returns the coverage depth (in number of reads) at the position pos (can be 110, 120 etc) grouped by Metadata 
  return(coverage[, .(depth = .N), by = .(Barcode, Genotype, Rep, Condition, Sample_name, Feature.name, Feature.strand, Degradation_Tag, pos)]) #% add Rep
}

# --- 3. DATA INITIALIZATION ---

new <- "AT3G60630"


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
#Col0_xrn4_trunc5 <- fread("flepseq2/fragment_analysis/diff_frag_5_truncated_Col0_xrn4_Stat.txt",header=TRUE, sep="\t")


target_barcodes <- SampleTable[Genotype %in% c("Col0", "xrn4") & Condition == "Seedling", Barcode] #Col0

barcodes <- SampleTable[, Barcode] 

# PART B: COVERAGE ANALYSIS (DNE1 TARGETS)
df_targets <- ds %>%
  dplyr::filter(
    Barcode %in% barcodes,
    sub("\\..*$", "", Feature.name) == new
  ) %>%
  dplyr::select(Barcode, Feature.name, Feature.strand, Read.start, Read.end, Degradation_Tag, RA5_tag) %>%
  collect() %>%
  as.data.table()

df_targets <- df_targets[!is.na(Feature.name)] #%
# Merge metadata and fix factors
df_targets <- merge(unique(df_targets), 
                    SampleTable[, .(Barcode, Genotype, Rep, Condition, Sample_name, libsize)], #% Add Rep
                    by = "Barcode", all.x = TRUE)

df_targets[, Degradation_Tag := factor(trimws(Degradation_Tag), levels = actual_tags)]


coverage_data <- calculate_coverage_dt(df_targets, bin_width = 10)

# Check the depth magnitude (Should be much less than 1 million now)
print(paste("Maximum Raw Depth:", max(coverage_data$depth, na.rm = TRUE)))

final_coverage <- coverage_data[
  ,
  .(depth = depth),  # no averaging
  by = .(Genotype, Rep, Condition, Sample_name, Barcode, Feature.name, Feature.strand, Degradation_Tag, pos)
]
final_coverage[, RepID := paste(Genotype, Condition, Rep, Sample_name, sep = "_")]


p2 <- ggplot(final_coverage, aes(x = pos, y = depth, color = Degradation_Tag)) +
  geom_line(linewidth = 1) +
  # Use "free" to allow both X and Y scales to adapt to each gene
  facet_grid(Genotype + RepID ~ Feature.name + Feature.strand, scales = "free_x") +  #% Add + Rep
  scale_color_manual(values = tag_colors, drop = FALSE) + 
  bio_journal_theme +
  # Use pretty_breaks so each unique scale gets appropriate labels
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), 
                     expand = expansion(mult = c(0, 0.1))) + 
  labs(title = "Coverage: ", 
       x = "Chromosome Coordinate", 
       y = "coverage_raw")
print(p2)

#ggsave("flepseq2/plots/LOM2 AT3G60630 miRNA target.png", p2, width = 12, height = 12)




