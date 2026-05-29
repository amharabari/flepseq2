# RA5 GLOBAL DIFFERENTIAL ANALYSIS OF FRAGMENTS IN Col0-xrn4
setwd("~/AMHLABO/")

library(data.table)
library(arrow)     # For Out-of-Memory processing
library(dplyr)     # For Arrow "remote control"
library(pbapply)   # For progress bars
library(ggplot2)
library(multcomp)
library(UpSetR)
library(patchwork)
library(ComplexUpset)

setDTthreads(0)    # Use all 12 threads


table_glht <- function(x) {
  pq     <- summary(x)$test
  mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
  colnames(mtests) <- c("Estimate", "Std.Error", "z.value", "p.value")
  return(mtests)
}

# ==============================================================================
# 1. INITIALIZATION & METADATA
# ==============================================================================
parquet_dir <- "flepseq2/parquet_dataset"

metadata_path <- "flepseq2/bamparsed/barcode_correspondance.tsv"
actual_tags   <- c("NO_RA5", "RA5")

ds <- open_dataset(parquet_dir)
SampleTable <- fread(metadata_path)

library_sizes <- ds %>%
  group_by(Barcode) %>%
  summarize(libsize = n()) %>%
  collect() %>% setDT()

SampleTable <- merge(SampleTable, library_sizes, by = "Barcode", all.x = TRUE)


# ==============================================================================
# 3. MAIN COMPARISON FUNCTION
# ==============================================================================

target_barcodes <- SampleTable[
  Genotype %in% c("Col0", "xrn4") & Condition == "Seedling", Barcode]

#target_barcodes <- SampleTable[
#  Genotype %in% c("Col0", "xrn4") & Condition == "Seedling", Barcode] #%


global_tag_raw <- ds %>%
  dplyr::select(Barcode, Feature.name, RA5_tag) %>%
  collect() %>%
  as.data.table()

# --- 3. Row filters ---
global_tag_raw <- global_tag_raw[Barcode %chin% target_barcodes]
global_tag_raw <- global_tag_raw[!grepl("^ATCG|^ATMG", Feature.name)]
global_tag_raw <- global_tag_raw[RA5_tag != "No_CDS"]
global_tag_raw <- global_tag_raw[RA5_tag %in% actual_tags]

# --- 4. Merge metadata  ---
global_tag_raw <- merge(
  global_tag_raw,
  SampleTable[, .(Barcode, Sample_name, Genotype, Rep, Condition)],
  by = "Barcode", all.x = TRUE
)

# --- 5. Gene filter (≥150 reads in all 3 reps per genotype) ---
gene_totals <- global_tag_raw[, .(Total = .N),
                              by = .(Feature.name, Sample_name, Genotype, Rep)]

list_targets <- gene_totals[, .(
  passed = sum(Total >= 150) == 3
), by = .(Feature.name, Genotype)]

list_targets <- list_targets[, .(both_pass = all(passed)), 
                             by = Feature.name][both_pass == TRUE, Feature.name]

cat("Genes passing flexible filter:", length(list_targets), "\n")
global_tag_raw <- global_tag_raw[Feature.name %in% list_targets]


# --- 7. Aggregate: lib totals and per-tag counts ---
lib_totals <- global_tag_raw[, .(lib_total = .N), by = Barcode]
tag_counts <- global_tag_raw[, .(raw_count = .N), by = .(Barcode, Sample_name, RA5_tag)]
tag_counts <- merge(tag_counts, lib_totals, by = "Barcode")

# --- 8. Merge Genotype and Rep onto tag_counts for GLM ---
tag_counts <- merge(
  tag_counts,
  SampleTable[, .(Barcode, Genotype, Rep, Condition)],
  by = "Barcode", all.x = TRUE
)

# --- 9. Derived columns ---
tag_counts[, NonTag_reads     := lib_total - raw_count]
tag_counts[, raw_count_adj    := raw_count    + 0.5]
tag_counts[, NonTag_reads_adj := NonTag_reads + 0.5]
tag_counts[, variable         := factor(Genotype, levels = c("Col0", "xrn4"))] #% 

# --- 10. GLM per tag ---
results <- lapply(actual_tags, function(tag) {
  tag_df <- tag_counts[RA5_tag == tag]
  if (nrow(tag_df) <= 2) {
    cat("  WARNING: Not enough observations for tag:", tag, "\n")
    return(NULL)
  }
  tryCatch({
    fit      <- glm(cbind(raw_count_adj, NonTag_reads_adj) ~ variable,
                    family = quasibinomial, data = tag_df)
    glht_obj <- glht(fit, mcp(variable = "Tukey"))
    Tab      <- as.data.frame(table_glht(summary(glht_obj)))
    Tab$comparaison     <- rownames(Tab)
    Tab$RA5_tag <- tag
    rownames(Tab)       <- NULL
    as.data.table(Tab)
  }, error = function(e) {
    cat("  ERROR for tag", tag, ":", conditionMessage(e), "\n")
    NULL
  })
})

names(results)  <- actual_tags
result_table    <- rbindlist(results, fill = TRUE)
result_table[, padj := p.adjust(p.value, "BH")]

# --- 11. Mean proportions per genotype per tag ---
tag_counts[, cpm := (raw_count / lib_total) * 1e6]

cpm_wide <- dcast(
  tag_counts,
  RA5_tag ~ Genotype + Sample_name,
  value.var = "cpm"
)

mean_cpm_summary <- tag_counts[, .(
  mean_cpm = mean(cpm, na.rm = TRUE),
  sd_cpm   = sd(cpm,   na.rm = TRUE)
), by = .(RA5_tag, Genotype)]

mean_cpm_wide <- dcast(mean_cpm_summary, RA5_tag ~ Genotype, value.var = "mean_cpm")
setnames(mean_cpm_wide,
         old = c("Col0", "xrn4"), #%
         new = c("mean_cpm_Col0", "mean_cpm_xrn4")) #%
mean_cpm_wide[, fold_change := mean_cpm_xrn4 / mean_cpm_Col0] #%


result_table <- merge(result_table, mean_cpm_wide, by = "RA5_tag", all.x = TRUE)
result_table <- merge(result_table, cpm_wide, by = "RA5_tag", all.x = TRUE)

result_table[, sig := ifelse(padj <= 0.05, "significant", "non_significant")]

cat("\n--- GLM Results ---\n")
print(result_table[, .(RA5_tag, Estimate, Std.Error, z.value, p.value, padj, sig,
                       mean_cpm_Col0, mean_cpm_xrn4, fold_change)]) #%

# --- 12. Boxplot (CPM, significant tags annotated) ---

# Per-replicate CPM for boxplot (one value per Barcode × tag)
box_df <- copy(tag_counts)
box_df[, RA5_tag := factor(RA5_tag, levels = actual_tags)]
box_df[, Genotype        := factor(Genotype,        levels = c("Col0", "xrn4"))] #
box_df[, Rep             := factor(Rep)]

# Bracket y position for significance stars
bracket_y <- box_df[, .(bracket_y = max(cpm, na.rm = TRUE) * 1.08),
                    by = RA5_tag]

sig_labels <- result_table[, .(RA5_tag, sig, padj)]
bracket_y  <- merge(bracket_y, sig_labels, by = "RA5_tag", all.x = TRUE)
p <- ggplot(box_df, aes(x = RA5_tag, y = cpm, fill = Genotype)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.65,
               outlier.shape = 19, outlier.size = 1.5, alpha = 0.85) +
  # Significance star above significant tags only
  geom_text(
    data = bracket_y[sig == "significant"],
    aes(x = RA5_tag, y = bracket_y, label = "*"),
    inherit.aes = FALSE, size = 7, fontface = "bold"
  ) +
  scale_y_log10( breaks = c(10000, 20000, 50000, 100000, 200000, 500000, 1000000),
                 labels = scales::comma,
                 limits = c(10000, 1000000)) +
  labs(
    title    = paste("Global Tag CPM:", "Col0", "vs", "xrn4"), #%
    subtitle = "log10 CPM per replicate | * padj ≤ 0.05",
    x        = "Degradation Tag",
    y        = "CPM (log10 scale)",
    fill     = "Genotype"
  ) +
  
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p)

#ggsave("flepseq2/fragment_analysis/global_frag_analysis_Col0_xrn4.pdf", p, width = 5, height = 5)


mean_lines <- box_df[, .(mean_cpm = mean(cpm, na.rm = TRUE)),
                     by = .(RA5_tag, Genotype)]

points <- ggplot(
  box_df[RA5_tag == "RA5"],
  aes(x = Genotype, y = cpm, color = Genotype)) +
  
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.4, linewidth = 0.7, show.legend = FALSE) +
  geom_point(aes(shape = Rep),
             position = position_jitter(width = 0.08, seed = 42),
             size = 4, alpha = 0.9) +
  geom_text(
    data = bracket_y[sig == "significant" & RA5_tag == "RA5"],
    aes(x = 1.5, y = bracket_y, label = "*"),
    inherit.aes = FALSE, size = 8, fontface = "bold", color = "black"
  ) +
  scale_shape_manual(values = c(15, 16, 17), name = "Replicate") +
  scale_color_manual(values = c("Col0" = "#4A90D9", "xrn4" = "#E87040")) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    title    = "RA5-tagged reads: Col0 vs xrn4",
    subtitle = "CPM per replicate | * padj ≤ 0.05 | crossbar = mean",
    x        = NULL, y = "CPM", color = "Genotype"
  ) +
  theme_bw(base_size = 13)

print(points)


ggsave("flepseq2/fragment_analysis/global_RA5_Col0_xrn4.pdf", points, width = 5, height = 5)

