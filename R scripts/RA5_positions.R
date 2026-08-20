library(dplyr)
library(data.table)
library(ggplot2)

dir.create("r_plots/objects", showWarnings = FALSE, recursive = TRUE)

SampleTable <- fread("provided/barcode_correspondance.tsv")
SampleTable <- SampleTable[, .(Barcode, Sample_name, Genotype, Condition, Rep)]

filter_A_for_all_BC <- function(df) {
  # Extraire toutes les combinaisons distinctes de Genotype et Rep
  all_BC <- df %>%
    distinct(Genotype, Rep)
  
  # Pour chaque Feature.name, compter combien de combinaisons (B,C) il couvre
  valid_A <- df %>%
    distinct(Genotype, Rep, Feature.name) %>%
    group_by(Feature.name) %>%
    summarise(n_bc = n()) %>%
    ungroup()
  
  # Nombre total de combinaisons Genotype, Rep dans la table
  total_bc <- nrow(all_BC)
  
  # Garder les Feature.name qui ont toutes les combinaisons Genotype, Rep
  valid_A <- valid_A %>%
    filter(n_bc == total_bc) %>%
    pull(Feature.name)
  
  # Filtrer le df original pour garder uniquement les lignes avec ces A valides
  df %>% filter(Feature.name %in% valid_A)
}


## 1. Parameters -----------------------------------------------------------

target_barcodes <- SampleTable[Genotype %in% c("Col0", "xrn4", "xrn4_dne1") & Condition == "Seedling", Barcode]
bin_width <- 0.02                      # relative-position bin size
x_min <- -0.25                           # plotting window
x_max <- 1.5

## 2. Prepare CDS lengths per transcript ----------------------------------

# cds_lengths <- ds %>%
#   dplyr::select(Feature.name, CDS.start_codon, CDS.stop_codon) %>%
#   distinct() %>%
#   collect() %>%
#   as.data.table()
# 
# # Positive CDS length in transcript coordinates
# cds_lengths[, cds_length := abs(CDS.stop_codon - CDS.start_codon)]
# cds_lengths <- cds_lengths[, .(Feature.name, cds_length)]

#fwrite(cds_lengths, "r_tables/cds_lengths_RA5.txt", row.names = FALSE )
cds_lengths <- fread("r_tables/cds_lengths_RA5.txt")


## 3. Extract relevant columns and join CDS lengths -----------------------

# rel_positions <- ds %>%
#   filter(
#     Barcode %in% target_barcodes,
#     !is.na(CDS_dist_5prime),
#     !is.na(CDS_dist_3prime)
#   ) %>%
#   dplyr::select(
#     Barcode,
#     Feature.name,
#     RA5_tag,
#     CDS_dist_5prime,
#     CDS_dist_3prime,
#     Feature.strand
#   ) %>%
#   collect() %>%
#   setDT()

#fwrite(rel_positions, "r_tables/rel_positions_RA5.txt", row.names = FALSE )
rel_positions <- fread("r_tables/rel_positions_RA5.txt")


# Join sample metadata
rel_positions <- merge(rel_positions,SampleTable,by = "Barcode")

RA5_genelist <-unique(rel_positions$Feature.name)
cat("Unique genes:", uniqueN(rel_positions$Feature.name), "\n")

rel_positions <- filter_A_for_all_BC(rel_positions)
cat("Unique genes:", uniqueN(rel_positions$Feature.name), "\n")

# Join CDS lengths
rel_positions[cds_lengths,
              cds_length := i.cds_length,
              on = .(Feature.name)]

# Drop any rows without CDS length
rel_positions <- rel_positions[!is.na(cds_length) & cds_length > 0]

## 4. Define relative CDS coordinate --------------------------------------

# Relative distance from CDS start: 0 = start codon, ~1 = stop codon
rel_positions[, rel_5prime := CDS_dist_5prime / cds_length]

# Bin the relative coordinate
rel_positions[, bin := round(rel_5prime / bin_width) * bin_width]

# Optional: restrict to a window around CDS
rel_positions <- rel_positions[bin >= x_min & bin <= x_max]

## 5. Aggregate counts per barcode / tag / bin ----------------------------

# Count reads in each bin per sample and RA5_tag
rel_positions_counts <- rel_positions[,.(n = .N), by = .(Barcode, Genotype, Rep, RA5_tag, bin)]


# Compute within-sample proportions (per RA5_tag)
rel_positions_counts[,total := sum(n), by = .(Barcode, Genotype, Rep, RA5_tag)]

rel_positions_counts[
  ,
  prop := n / total
]

## 6. Average across replicates / genotypes -------------------------------

rel_positions_summary <- rel_positions_counts[
  ,
  .(mean_prop = mean(prop)),
  by = .(Genotype, Rep, RA5_tag, bin)
]

## 7. Plot CDS-relative profiles ------------------------------------------

p <- ggplot(
  rel_positions_summary,
  aes(x = bin, y = mean_prop, color = RA5_tag, fill = RA5_tag)
) +
  geom_area(alpha = 0.3, position = "identity") +
  geom_line(linewidth = 0.6) +
  # 0 = CDS start, 1 ≈ CDS stop
  geom_vline(xintercept = c(0, 1), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("RA5" = "#DA27F5", "NO_RA5" = "#696969")) +
  scale_fill_manual(values  = c("RA5" = "#DA27F5", "NO_RA5" = "#696969")) +
  scale_x_continuous(breaks = seq(floor(min(rel_positions_summary$bin)),
                                  ceiling(max(rel_positions_summary$bin)),
                                  by = 0.2)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
  facet_grid(Rep ~ Genotype) +
  theme_bw(base_size = 15) +
  theme(
    panel.spacing.x = unit(1.2, "lines")) +
  labs(
    title = "Distribution of read 5′ ends relative to the CDS",
    x = "Relative position (CDS start to stop)",
    y = "Proportion of reads by position"
  )

print(p)

saveRDS(p, "r_plots/objects/rel_position_RA5v2.rds")
ggsave("r_plots/rel_position_RA5v2.svg", p, width = 10, height = 6, scale = 1.3)
ggsave("r_plots/rel_position_RA5v2.png", p, width = 10, height = 6, scale = 1.3)

#---------------------------------------------------------------------------

