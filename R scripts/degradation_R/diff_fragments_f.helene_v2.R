
library(data.table)
library(arrow)     # For Out-of-Memory processing
library(dplyr)     # For Arrow "remote control"
library(pbapply)   # For progress bars
library(ggplot2)
library(multcomp)
library(UpSetR)
library(patchwork)
setDTthreads(0)    # Use all 12 threads
# ==============================================================================
# 1. INITIALIZATION & METADATA
# ==============================================================================
actual_tags   <- c("Intact", "5'_truncated", "3'_truncated", "Both_truncated")
SampleTable <- fread("1_Mapping_Flepseq2/barcode_correspondance.tsv")

bam_dir <- "2_Output"
bam_files <- list.files(
  path    = bam_dir,
  pattern = "_bamparsed\\.tsv$",
  full.names = TRUE
)

bamparsed_all <- rbindlist(
  lapply(bam_files, function(f) {
    dt <- fread(f)
    dt[, Barcode := sub(".*_(barcode[0-9]+)_bamparsed\\.tsv$", "\\1", basename(f))]
    dt
  }),
  use.names = TRUE,
  fill      = TRUE
)
library_sizes <- bamparsed_all %>%
  group_by(Barcode) %>%
  summarize(libsize = n()) %>%
  collect() %>% setDT()

SampleTable <- merge(SampleTable, library_sizes, by = "Barcode", all.x = TRUE)

# ==============================================================================
# 2. HELPER FUNCTIONS 
# ==============================================================================

my_prepare <- function(x, level1, level2) {
  x <- as.data.table(x)  
  x$variable <- x$Genotype
  x$variable  <- factor(x$variable, levels = c(level1, level2))
  x$Tag_reads    <- x$raw_count
  x$NonTag_reads <- x$gene_raw_total - x$raw_count
  x$Tag_reads_adj    <- x$Tag_reads    + 0.5
  x$NonTag_reads_adj <- x$NonTag_reads + 0.5
  x$prop_tag <- x$Tag_reads_adj / (x$Tag_reads_adj + x$NonTag_reads_adj)
  x$Percent_tag_co <- x$prop_tag * 100
  return(x)
}

filter_A_for_all_BC <- function(df) {
  # Force to data.table to avoid dplyr/Arrow src conflict
  df <- as.data.table(df)
  all_BC <- unique(df[, .(Genotype, Rep)])
  total_bc <- nrow(all_BC)
  # Count how many distinct (Genotype, Rep) combinations each gene covers
  coverage <- unique(df[, .(Genotype, Rep, Feature.name)])[
    , .(n_bc = .N), by = Feature.name]
  valid_genes <- coverage[n_bc == total_bc, Feature.name]
  df[Feature.name %in% valid_genes]
}

table_glht <- function(x) {
  pq     <- summary(x)$test
  mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
  colnames(mtests) <- c("Estimate", "Std.Error", "z.value", "p.value")
  return(mtests)
}

my_glm_test <- function(x, level1, level2) {
  x <- as.data.table(x)
  results_by_condition <- list()
  
  for (t in unique(x$Condition)) {
    Selec_condition <- as.data.table(x[Condition == t])
    Selec_condition <- filter_A_for_all_BC(Selec_condition)
    gene_list <- unique(Selec_condition$Feature.name)
    
    tab_list <- vector("list", length(gene_list))
    
    for (i in seq_along(gene_list)) {
      if (i %% 500 == 0) cat("  Progress:", i, "/", length(gene_list), "\n")
      gene   <- gene_list[i]
      perTag <- Selec_condition[Feature.name == gene]
      if (nrow(perTag) <= 2) next
      tryCatch({
        Frag     <- glm(cbind(Tag_reads_adj, NonTag_reads_adj) ~ variable,
                        family = quasibinomial, data = perTag)
        glht_obj <- summary(glht(Frag, mcp(variable = "Tukey")))  # ← missing line
        Tab      <- as.data.frame(table_glht(glht_obj))
        Tab$comparaison  <- row.names(Tab)
        Tab$Feature.name <- gene
        rownames(Tab)    <- NULL
        tab_list[[i]]    <- as.data.table(Tab)
      }, error = function(e) NULL)
    }   # ← closes for(i in seq_along)
    
    Result_table <- rbindlist(tab_list, fill = TRUE, use.names = TRUE)
    Result_table[, padj      := p.adjust(p.value, "BH")]
    Result_table[, Condition := t]
    results_by_condition[[t]] <- Result_table
  }   # ← closes for(t in ...)
  
  rbindlist(results_by_condition, fill = TRUE)
}   # ← closes my_glm_test


strip_isoform <- function(x) sub("\\.[0-9]+$", "", x)

# ==============================================================================
# 3. MAIN COMPARISON FUNCTION
# ==============================================================================

level1 <- "Col0"
level2 <- "xrn4"
comp_name <- paste0(level1, "_", level2)

# 3.1 target barcodes
target_barcodes <- SampleTable[
  Genotype %in% c(level1, level2) & Condition == "Seedling", Barcode
]

# 3.2 global counts with polyA columns
global_counts <- bamparsed_all %>%
  dplyr::select(
    Barcode,
    Feature.name,
    Degradation_Tag,
    PolyA,
    PolyA_length_signal,
    Add_tail
  ) %>%
  collect() %>%
  as.data.table()

global_counts <- global_counts[Barcode %chin% target_barcodes]
global_counts <- global_counts[!grepl("^ATCG|^ATMG", Feature.name)]
global_counts <- global_counts[Degradation_Tag != "No_CDS"]

table_before <- table(global_counts$Degradation_Tag)


# enforce: 3'-truncated reads should NOT have polyA
global_counts[, polya_len := as.numeric(PolyA_length_signal)]

global_counts <- global_counts[
  !(
    Degradation_Tag == "3'_truncated" & (
      (!is.na(PolyA)    & PolyA    != "") |
        (!is.na(Add_tail) & Add_tail != "") |
        (!is.na(polya_len) & polya_len >= 9)
    )
  )
]
table_after  <- table(global_counts$Degradation_Tag)

table_before
table_after
# aggregate counts per tag
global_counts <- global_counts[, .(raw_count = .N),
                               by = .(Barcode, Feature.name, Degradation_Tag)]

# 3.3 add metadata
global_counts <- merge(
  global_counts,
  SampleTable[, .(Barcode, Sample_name, Genotype, Rep, libsize)],
  by = "Barcode", all.x = TRUE
)

# 4.1 flexible filter ≥150 reads in all 3 reps per genotype
gene_totals <- global_counts[, .(Total = sum(raw_count)),
                             by = .(Feature.name, Sample_name, Genotype, Rep)]

list_targets <- gene_totals[, .(
  passed = sum(Total >= 50) == 3
), by = .(Feature.name, Genotype)][passed == TRUE, unique(Feature.name)]

cat("Genes passing flexible filter:", length(list_targets), "\n")

# 4.2 final stats table
final_stats <- global_counts[Feature.name %in% list_targets]
final_stats[, gene_raw_total := sum(raw_count), by = .(Sample_name, Feature.name)]
final_stats[, prop_tag       := raw_count / pmax(gene_raw_total, 1)]
final_stats[, cpm            := (raw_count / libsize) * 1e6]
final_stats[, Condition      := "Seedling"]

delta_threshold <- 0.

results_list <- lapply(actual_tags, function(tag) {
  tag_df <- my_prepare(final_stats[Degradation_Tag == tag], level1, level2)
  Stat   <- my_glm_test(tag_df)
  
  if (nrow(Stat) == 0) {
    cat("WARNING: No results for tag:", tag, "\n")
    return(NULL)
  }
  
  mean_props <- tag_df[, .(mean_prop = mean(prop_tag, na.rm = TRUE)),
                       by = .(Feature.name, variable)]
  mean_props <- dcast(mean_props, Feature.name ~ variable, value.var = "mean_prop")
  setnames(mean_props, c("Feature.name", level1, level2))
  mean_props[, delta_prop := abs(get(level1) - get(level2))]
  
  Stat <- as.data.table(Stat)
  Stat <- merge(Stat,
                mean_props[, .(Feature.name, get(level1), get(level2), delta_prop)],
                by = "Feature.name", all.x = TRUE)
  
  Stat$diff   <- "non_significant"
  Stat[padj <= 0.05 & delta_prop >= delta_threshold, diff := "significant"]
  Stat$Change <- "no-change"
  Stat[Estimate < 0, Change := "decrease"]
  Stat[Estimate > 0, Change := "increase"]
  Stat$diff            <- factor(Stat$diff)
  Stat$Change          <- factor(Stat$Change)
  Stat$Degradation_Tag <- tag
  
  Count_tag <- Stat[, .(number = .N),
                    by = .(comparaison, Condition, diff, Change)]
  Count_tag[, total   := sum(number), by = .(comparaison, Condition)]
  Count_tag[, Percent := 100 * number / total]
  Count_tag[, Degradation_Tag := tag]
  
  Stat_sig <- Stat[diff == "significant"]
  cat("Significant genes for", tag, ":", nrow(Stat_sig), "\n")
  
  list(Stat = Stat, Stat_sig = Stat_sig, Count = Count_tag)
})

names(results_list) <- actual_tags

all_sig    <- rbindlist(lapply(results_list, `[[`, "Stat_sig"), fill = TRUE)
all_counts <- rbindlist(lapply(results_list, `[[`, "Count"),    fill = TRUE)

cat("\nSignificant genes per tag:\n")
print(all_sig[, .N, by = Degradation_Tag])
print(all_counts)