# GLOBAL DIFFERENTIAL ANALYSIS OF FRAGMENTS IN Col0-xnr4

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
actual_tags   <- c("Intact", "5'_truncated", "3'_truncated", "Both_truncated")

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
  Genotype %in% c("xrn4", "xrn4_dne1") & Condition == "Seedling", Barcode]

#target_barcodes <- SampleTable[
#  Genotype %in% c("Col0", "xrn4") & Condition == "Seedling", Barcode] #%


global_tag_raw <- ds %>%
  dplyr::select(Barcode, Feature.name, Degradation_Tag) %>%
  collect() %>%
  as.data.table()

# --- 3. Row filters (while Feature.name still exists) ---
global_tag_raw <- global_tag_raw[Barcode %chin% target_barcodes]
global_tag_raw <- global_tag_raw[!grepl("^ATCG|^ATMG", Feature.name)]
global_tag_raw <- global_tag_raw[Degradation_Tag != "No_CDS"]
global_tag_raw <- global_tag_raw[Degradation_Tag %in% actual_tags]

# --- 4. Merge metadata (needed by the gene filter below) ---
global_tag_raw <- merge(
  global_tag_raw,
  SampleTable[, .(Barcode, Sample_name, Genotype, Rep, Condition)],
  by = "Barcode", all.x = TRUE
)

# --- 5. Flexible gene filter (≥150 reads in all 3 reps per genotype) ---
# Each row is still one read here, so use .N
gene_totals <- global_tag_raw[, .(Total = .N),
                              by = .(Feature.name, Sample_name, Genotype, Rep)]

list_targets <- gene_totals[, .(
  passed = sum(Total >= 150) == 3
), by = .(Feature.name, Genotype)]

list_targets <- list_targets[, .(both_pass = all(passed)), 
                             by = Feature.name][both_pass == TRUE, Feature.name]

cat("Genes passing flexible filter:", length(list_targets), "\n")
global_tag_raw <- global_tag_raw[Feature.name %in% list_targets]

# --- 6. Drop Feature.name — no longer needed ---
global_tag_raw[, Feature.name := NULL]

# --- 7. Aggregate: lib totals and per-tag counts (each row = 1 read → .N is correct) ---
lib_totals <- global_tag_raw[, .(lib_total = .N), by = Barcode]
tag_counts <- global_tag_raw[, .(raw_count = .N), by = .(Barcode, Sample_name, Degradation_Tag)]
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
tag_counts[, variable         := factor(Genotype, levels = c("xrn4", "xrn4_dne1"))] #% 

# --- 10. GLM per tag ---
results <- lapply(actual_tags, function(tag) {
  tag_df <- tag_counts[Degradation_Tag == tag]
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
    Tab$Degradation_Tag <- tag
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

mean_cpm_summary <- tag_counts[, .(
  mean_cpm = mean(cpm, na.rm = TRUE),
  sd_cpm   = sd(cpm,   na.rm = TRUE)
), by = .(Degradation_Tag, Genotype)]

mean_cpm_wide <- dcast(mean_cpm_summary, Degradation_Tag ~ Genotype, value.var = "mean_cpm")
setnames(mean_cpm_wide,
         old = c("xrn4", "xrn4_dne1"), #%
         new = c("mean_cpm_xrn4", "mean_cpm_xrn4_dne1")) #%
mean_cpm_wide[, fold_change := mean_cpm_xrn4_dne1 / mean_cpm_xrn4] #%

result_table <- merge(result_table, mean_cpm_wide, by = "Degradation_Tag", all.x = TRUE)
result_table[, sig := ifelse(padj <= 0.05, "significant", "non_significant")]

cat("\n--- GLM Results ---\n")
print(result_table[, .(Degradation_Tag, Estimate, Std.Error, z.value, p.value, padj, sig,
                       mean_cpm_xrn4, mean_cpm_xrn4_dne1, fold_change)]) #%

# --- 12. Boxplot (CPM, significant tags annotated) ---

# Per-replicate CPM for boxplot (one value per Barcode × tag)
box_df <- copy(tag_counts)
box_df[, Degradation_Tag := factor(Degradation_Tag, levels = actual_tags)]
box_df[, Genotype        := factor(Genotype,        levels = c("xrn4", "xrn4_dne1"))] #%

# Bracket y position for significance stars
bracket_y <- box_df[, .(bracket_y = max(cpm, na.rm = TRUE) * 1.08),
                    by = Degradation_Tag]

sig_labels <- result_table[, .(Degradation_Tag, sig, padj)]
bracket_y  <- merge(bracket_y, sig_labels, by = "Degradation_Tag", all.x = TRUE)

p <- ggplot(box_df, aes(x = Degradation_Tag, y = cpm,
                        fill = Genotype)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.65,
               outlier.shape = 19, outlier.size = 1.5, alpha = 0.85) +
  # Significance star above significant tags only
  geom_text(
    data = bracket_y[sig == "significant"],
    aes(x = Degradation_Tag, y = bracket_y, label = "*"),
    inherit.aes = FALSE, size = 7, fontface = "bold"
  ) +
  scale_y_log10( breaks = c(10000, 20000, 50000, 100000, 200000, 500000, 1000000),
                 labels = scales::comma,
                 limits = c(10000, 1000000)) +
  labs(
    title    = paste("Global Tag CPM:", "xrn4", "vs", "xrn4_dne1"), #%
    subtitle = "log10 CPM per replicate | * padj ≤ 0.05",
    x        = "Degradation Tag",
    y        = "CPM (log10 scale)",
    fill     = "Genotype"
  ) +

  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p)

#ggsave("flepseq2/fragment_analysis/global_frag_analysis_Col0_xrn4.pdf", p, width = 5, height = 5)


#p <- ggplot(box_df, aes(x = Degradation_Tag, y = cpm,
#                        fill = Genotype)) +
#  geom_violin(position = position_dodge(0.8), width = 0.65,
#              quantile.linetype = F, alpha = 0.85) +
#  
#  # Significance star above significant tags only
#  geom_text(
#    data = bracket_y[sig == "significant"],
#    aes(x = Degradation_Tag, y = bracket_y, label = "*"),
#    inherit.aes = FALSE, size = 7, fontface = "bold"
#  ) +
#  scale_y_log10( breaks = c(10000, 20000, 50000, 100000, 200000, 500000, 1000000),
#                 labels = scales::comma,
#                 limits = c(10000, 1000000)) +
#  labs(
#    title    = paste("Global Tag CPM:", "Col0", "vs", "xrn4"), #%
#    subtitle = "log10 CPM per replicate | * padj ≤ 0.05",
#    x        = "Degradation Tag",
#    y        = "CPM (log10 scale)",
#    fill     = "Genotype"
#  ) +
#  
#  theme_bw(base_size = 13) +
#  theme(axis.text.x = element_text(angle = 30, hjust = 1))
#print(p)