# ==============================================================================
# RA5 GENE BY GENE
# ==============================================================================
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

# Folders
raw_path    <- "flepseq2/bamparsed/"
parquet_dir <- "flepseq2/parquet_dataset"
if(!dir.exists(parquet_dir)) dir.create(parquet_dir, recursive = TRUE)

cols_to_keep <- unique(c(
  "Read.name", "Read.start", "Read.end",
  "Feature.name", "Feature.start", "Feature.end", "Feature.strand",
  "CDS.start_codon", "CDS.stop_codon", "CDS_dist_5prime", "CDS_dist_3prime", 
  "Degradation_Tag", "RA5_tag", "PolyA", "Add_tail", "PolyA_length_basecall", 
  "PolyA_length_signal", "Add_tail_length", "Dist_from_5prime", 
  "Dist_from_3prime", "Tail_tag"
))

# ==============================================================================
# IMPORTING BAMPARSED:
# PHASE 1: CONVERSION 
# ==============================================================================

file_paths <- list.files(path = raw_path, pattern = "\\.gz$", full.names = TRUE)
start_total <- Sys.time()

# Read just the header of one raw file
header <- fread(file_paths[1], sep = "\t", nrows = 0)
colnames(header)

cat("\n--- PHASE 1: Converting TSV.GZ to Parquet (Memory-Safe) ---\n")

for (i in seq_along(file_paths)) {
  f <- file_paths[i]
  bc_name <- regmatches(f, regexpr("barcode\\d+", f))
  out_path <- file.path(parquet_dir, paste0(bc_name, ".parquet"))
  
  if(!file.exists(out_path)) {
    t1 <- Sys.time()
    
    # Read one file at a time
    temp_dt <- fread(f, sep = "\t", select = cols_to_keep, showProgress = FALSE)
    temp_dt[, Barcode := bc_name]
    
    # Save as Parquet (Columnar format)
    write_parquet(temp_dt, out_path, compression = "snappy")
    
    # Reporting
    t2 <- Sys.time()
    obj_gb  <- round(object.size(temp_dt) / 1024^3, 2)
    sys_gb  <- round(sum(gc()[, 2]) / 1024, 2)
    elapsed <- round(as.numeric(t2 - t1), 1)
    
    cat(sprintf("[%d/%d] %s | Object: %s GB | Sys RAM: %s GB | Time: %s s\n", 
                i, length(file_paths), bc_name, obj_gb, sys_gb, elapsed))
    
    # Forced clean up
    rm(temp_dt); gc(verbose = FALSE)
  } else {
    cat(sprintf("[%d/%d] Skipping %s (Already converted)\n", i, length(file_paths), bc_name))
  }
}

# ==============================================================================
# PHASE 1b: PATCH — add RA5_tag to existing Parquet files
# ==============================================================================
cat("\n--- PHASE 1b: Patching RA5_tag into Parquet files ---\n")

for (i in seq_along(file_paths)) {
  f        <- file_paths[i]
  bc_name  <- regmatches(f, regexpr("barcode\\d+", f))
  out_path <- file.path(parquet_dir, paste0(bc_name, ".parquet"))
  
  existing_schema <- schema(open_dataset(out_path))
  
  if ("RA5_tag" %in% names(existing_schema)) {
    cat(sprintf("[%d/%d] Skipping %s (RA5_tag already present)\n",
                i, length(file_paths), bc_name))
    next
  }
  
  patch    <- fread(f, sep = "\t", select = c("Read.name", "RA5_tag"), showProgress = FALSE)
  patch    <- unique(patch, by = "Read.name")
  existing <- read_parquet(out_path)
  updated  <- merge(existing, patch, by = "Read.name", all.x = TRUE)
  write_parquet(updated, out_path, compression = "snappy")
  
  rm(patch, existing, updated); gc(verbose = FALSE)
  cat(sprintf("[%d/%d] Patched %s\n", i, length(file_paths), bc_name))
}

ds <- open_dataset(parquet_dir)
cat("Columns:", paste(names(schema(ds)), collapse = ", "), "\n")

# ==============================================================================
# PHASE 2: STACKING DATASETS + SAMPLE TABLE MAPPING
# ==============================================================================

cat("\n--- PHASE 2: Stacking and Mapping Datasets to Sample table---\n")

# Connect to all 12 parquet files at once. This uses NO memory.
ds <- open_dataset(parquet_dir)

cat("✅ Dataset mapped! Total rows available:", format(nrow(ds), big.mark=","), "\n")

# ==============================================================================
# 1. INITIALIZATION & METADATA
# ==============================================================================
metadata_path <- "flepseq2/bamparsed/barcode_correspondance.tsv"
actual_tags   <- c("NO_RA5", "RA5")

ds <- open_dataset(parquet_dir)
colnames(ds)
SampleTable <- fread(metadata_path)

library_sizes <- ds %>%
  group_by(Barcode) %>%
  summarize(libsize = n()) %>%
  collect() %>% setDT()

SampleTable <- merge(SampleTable, library_sizes, by = "Barcode", all.x = TRUE)


# ==============================================================================
# 2. HELPER FUNCTIONS (defined once)
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
    
    # ── GUARD: skip BH adjustment if no results came back ──────────────────
    if (nrow(Result_table) == 0 || !"p.value" %in% colnames(Result_table)) {
      cat("  WARNING: No usable results for condition:", t, "\n")
      next
    }
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
target_barcodes <- SampleTable[Genotype %in% c("Col0", "xrn4", "xrn4_dne1") & Condition == "Seedling", Barcode]

global_counts <- ds %>%
  dplyr::select(Barcode, Feature.name, RA5_tag) %>%
  collect() %>%
  as.data.table()

#global_counts[RA5_tag == "RA5", .N, by = Barcode]

#global_counts[RA5_tag == "RA5", .N, by = .(Genotype, RA5_tag)]

# What do raw RA5 proportions look like per gene?
#global_counts[RA5_tag == "RA5", .(total_RA5 = sum(raw_count)), 
#              by = Feature.name][order(-total_RA5)][1:20]

global_counts <- global_counts[Barcode %chin% target_barcodes]
global_counts <- global_counts[!grepl("^ATCG|^ATMG", Feature.name)]
global_counts <- global_counts[RA5_tag != "No_CDS"]
global_counts <- global_counts[, .(raw_count = .N),
                               by = .(Barcode, Feature.name, RA5_tag)]

# ==============================================================================
# 4. METADATA MERGE
# ==============================================================================
global_counts <- merge(global_counts,
                       SampleTable[, .(Barcode, Sample_name, Genotype, Rep, libsize)],
                       by = "Barcode", all.x = TRUE)

# ==============================================================================
# 5. FLEXIBLE FILTER (≥150 reads in all 3 reps per genotype)
# ==============================================================================
gene_totals <- global_counts[, .(Total = sum(raw_count)),
                             by = .(Feature.name, Sample_name, Genotype, Rep)]

list_targets <- gene_totals[, .(
  passed = sum(Total >= 50) == 3
), by = .(Feature.name, Genotype)][passed == TRUE, unique(Feature.name)]

cat("Genes passing flexible filter:", length(list_targets), "\n")

# ==============================================================================
# 6. FINAL STATS
# ==============================================================================
final_stats <- global_counts[Feature.name %in% list_targets]
final_stats[, gene_raw_total := sum(raw_count), by = .(Sample_name, Feature.name)]
final_stats[, prop_tag       := raw_count / pmax(gene_raw_total, 1)]
final_stats[, cpm            := (raw_count / libsize) * 1e6]
final_stats[, Condition      := "Seedling"]

# ==============================================================================
# 7. GLM LOOP OVER TAGS
# ==============================================================================
results_list <- lapply(actual_tags, function(tag) {
  
  tag_df <- my_prepare(final_stats[RA5_tag == tag], level1, level2)
  Stat   <- my_glm_test(tag_df, level1, level2)
  
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
                mean_props[, c("Feature.name", level1, level2, "delta_prop"), with = FALSE],
                by = "Feature.name", all.x = TRUE)
  
  # Run this BEFORE the diff assignment, inside the lapply, after the merge:
  cat("Columns in Stat:", paste(colnames(Stat), collapse=", "), "\n")
  cat("padj range:", range(Stat$padj, na.rm=TRUE), "\n")
  cat("delta_prop range:", range(Stat$delta_prop, na.rm=TRUE), "\n")
  cat("Rows with padj <= 0.05:", sum(Stat$padj <= 0.05, na.rm=TRUE), "\n")
  cat("Rows with delta_prop >=", delta_threshold, ":", sum(Stat$delta_prop >= delta_threshold, na.rm=TRUE), "\n")
  
  Stat$diff   <- "non_significant"
  Stat[padj <= 0.05 & delta_prop >= delta_threshold, diff := "significant"]
  Stat$Change <- "no-change"
  Stat[Estimate < 0, Change := "decrease"]
  Stat[Estimate > 0, Change := "increase"]
  Stat$diff    <- as.factor(Stat$diff)
  Stat$Change  <- as.factor(Stat$Change)
  Stat$RA5_tag <- tag
  
  Count_tag <- Stat[, .(number = .N), by = .(comparaison, Condition, diff, Change)]
  Count_tag[, total   := sum(number), by = .(comparaison, Condition)]
  Count_tag[, Percent := 100 * number / total]
  Count_tag[, RA5_tag := tag]
  
  Stat_sig <- Stat[diff == "significant"]
  
  safe_tag <- gsub("'", "", tag)
 # write.table(Stat_sig,
  #            paste0("flepseq2/fragment_analysis/diff_frag_", safe_tag, "_Col0_xrn4_Stat.txt"),
   #           sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Significant genes for", tag, ":", nrow(Stat_sig), "\n")
  
  list(Stat = Stat, Stat_sig = Stat_sig, Count = Count_tag)
})

names(results_list) <- actual_tags

# ==============================================================================
# 8. COMBINED TABLES
# ==============================================================================
all_sig    <- rbindlist(lapply(results_list, `[[`, "Stat_sig"), fill = TRUE)
all_counts <- rbindlist(lapply(results_list, `[[`, "Count"),    fill = TRUE)

cat("\nSignificant genes per tag:\n")
print(all_sig[, .N, by = RA5_tag])


final_stats[RA5_tag == "RA5", .N, by = Feature.name][order(-N)][1:20]
# How many genes survive the filter per tag?
final_stats[, .N, by = RA5_tag]

