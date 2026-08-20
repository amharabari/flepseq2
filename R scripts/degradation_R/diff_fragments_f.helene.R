# ==============================================================================
# 21/04/2026 - FLEPSEQ FRAGMENTATION ANALYSIS
# ==============================================================================

library(data.table)
library(arrow)     # For Out-of-Memory processing
library(dplyr)     # For Arrow "remote control"
library(pbapply)   # For progress bars
library(ggplot2)
library(multcomp)
library(UpSetR)
library(patchwork)


# Folders
parquet_dir <- "parquet_dataset"

# if(!dir.exists(parquet_dir)) dir.create(parquet_dir, recursive = TRUE)
# 
# cols_to_keep <- unique(c(
#   "Read.name", "Read.start", "Read.end",
#   "Feature.name", "Feature.start", "Feature.end", "Feature.strand",
#   "CDS.start_codon", "CDS.stop_codon", "CDS_dist_5prime", "CDS_dist_3prime", 
#   "Degradation_Tag", "RA5_tag", "PolyA", "Add_tail", "PolyA_length_basecall", 
#   "PolyA_length_signal", "Add_tail_length", "Dist_from_5prime", 
#   "Dist_from_3prime", "Tail_tag"
# ))
# 
# # ==============================================================================
# # IMPORTING BAMPARSED:
# # PHASE 1: CONVERSION 
# # ==============================================================================
# 
# file_paths <- list.files(path = raw_path, pattern = "\\.gz$", full.names = TRUE)
# start_total <- Sys.time()
# 
# # Read just the header of one raw file
# header <- fread(file_paths[1], sep = "\t", nrows = 0)
# colnames(header)
# 
# cat("\n--- PHASE 1: Converting TSV.GZ to Parquet (Memory-Safe) ---\n")
# 
# for (i in seq_along(file_paths)) {
#   f <- file_paths[i]
#   bc_name <- regmatches(f, regexpr("barcode\\d+", f))
#   out_path <- file.path(parquet_dir, paste0(bc_name, ".parquet"))
#   
#   if(!file.exists(out_path)) {
#     t1 <- Sys.time()
#     
#     # Read one file at a time
#     temp_dt <- fread(f, sep = "\t", select = cols_to_keep, showProgress = FALSE)
#     temp_dt[, Barcode := bc_name]
#     
#     # Save as Parquet (Columnar format)
#     write_parquet(temp_dt, out_path, compression = "snappy")
#     
#     # Reporting
#     t2 <- Sys.time()
#     obj_gb  <- round(object.size(temp_dt) / 1024^3, 2)
#     sys_gb  <- round(sum(gc()[, 2]) / 1024, 2)
#     elapsed <- round(as.numeric(t2 - t1), 1)
#     
#     cat(sprintf("[%d/%d] %s | Object: %s GB | Sys RAM: %s GB | Time: %s s\n", 
#                 i, length(file_paths), bc_name, obj_gb, sys_gb, elapsed))
#     
#     # Forced clean up
#     rm(temp_dt); gc(verbose = FALSE)
#   } else {
#     cat(sprintf("[%d/%d] Skipping %s (Already converted)\n", i, length(file_paths), bc_name))
#   }
# }
# 
# # ==============================================================================
# # PHASE 1b: PATCH — add RA5_tag to existing Parquet files
# # ==============================================================================
# cat("\n--- PHASE 1b: Patching RA5_tag into Parquet files ---\n")
# 
# for (i in seq_along(file_paths)) {
#   f        <- file_paths[i]
#   bc_name  <- regmatches(f, regexpr("barcode\\d+", f))
#   out_path <- file.path(parquet_dir, paste0(bc_name, ".parquet"))
#   
#   existing_schema <- schema(open_dataset(out_path))
#   
#   if ("RA5_tag" %in% names(existing_schema)) {
#     cat(sprintf("[%d/%d] Skipping %s (RA5_tag already present)\n",
#                 i, length(file_paths), bc_name))
#     next
#   }
#   
#   patch    <- fread(f, sep = "\t", select = c("Read.name", "RA5_tag"), showProgress = FALSE)
#   patch    <- unique(patch, by = "Read.name")
#   existing <- read_parquet(out_path)
#   updated  <- merge(existing, patch, by = "Read.name", all.x = TRUE)
#   write_parquet(updated, out_path, compression = "snappy")
#   
#   rm(patch, existing, updated); gc(verbose = FALSE)
#   cat(sprintf("[%d/%d] Patched %s\n", i, length(file_paths), bc_name))
# }
# 
# ds <- open_dataset(parquet_dir)
# cat("Columns:", paste(names(schema(ds)), collapse = ", "), "\n")
# 
# # ==============================================================================
# # PHASE 2: STACKING DATASETS + SAMPLE TABLE MAPPING
# # ==============================================================================
# 
# cat("\n--- PHASE 2: Stacking and Mapping Datasets to Sample table---\n")
# 
# # Connect to all 12 parquet files at once. This uses NO memory.
# ds <- open_dataset(parquet_dir)
# 
# cat("✅ Dataset mapped! Total rows available:", format(nrow(ds), big.mark=","), "\n")

# ==============================================================================
# 1. INITIALIZATION & METADATA
# ==============================================================================
actual_tags   <- c("Intact", "5'_truncated", "3'_truncated", "Both_truncated")

ds <- open_dataset(parquet_dir)
colnames(ds)
SampleTable <- fread("1_Mapping_Flepseq2/barcode_correspondance.tsv")

library_sizes <- ds %>%
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

run_comparison <- function(level1, level2, SampleTable, ds,
                           actual_tags, delta_threshold = 0) {
  comp_name <- paste0(level1, "_", level2)
  
  # --- 3.1 Target barcodes ---
  target_barcodes <- SampleTable[Genotype %in% c(level1, level2) & Condition == "Seedling", Barcode]
  
  # --- 3.2 Global counts ---
  global_counts <- ds %>%
    dplyr::select(Barcode, Feature.name, Degradation_Tag) %>%   # pull only needed columns
    collect() %>%
    as.data.table()
  
  global_counts <- global_counts[Barcode %chin% target_barcodes]
  global_counts <- global_counts[!grepl("^ATCG|^ATMG", Feature.name)]
  global_counts <- global_counts[Degradation_Tag != "No_CDS"]
  global_counts <- global_counts[, .(raw_count = .N),
                                 by = .(Barcode, Feature.name, Degradation_Tag)]
  # --- 3.3 Metadata merge ---
  global_counts <- merge(global_counts,
                         SampleTable[, .(Barcode, Sample_name, Genotype, Rep, libsize)],
                         by = "Barcode", all.x = TRUE)
  
  # --- 3.4 Flexible filter (≥150 reads in all 3 reps per genotype) ---
  gene_totals <- global_counts[, .(Total = sum(raw_count)),
                               by = .(Feature.name, Sample_name, Genotype, Rep)]
  
  list_targets <- gene_totals[, .(
    passed = sum(Total >= 50) == 3
  ), by = .(Feature.name, Genotype)][passed == TRUE, unique(Feature.name)]
  
  cat("Genes passing flexible filter:", length(list_targets), "\n")
  
  # --- 3.5 Final stats ---
  final_stats <- global_counts[Feature.name %in% list_targets]
  final_stats[, gene_raw_total := sum(raw_count), by = .(Sample_name, Feature.name)]
  final_stats[, prop_tag       := raw_count / pmax(gene_raw_total, 1)]
  final_stats[, cpm            := (raw_count / libsize) * 1e6]
  final_stats[, Condition      := "Seedling"]
  
  # --- 3.6 GLM loop over all tags ---
  results_list <- lapply(actual_tags, function(tag) {
    
    tag_df <- my_prepare(final_stats[Degradation_Tag == tag], level1, level2)
    Stat <- my_glm_test(tag_df, level1, level2)    
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
    Stat <- merge(Stat, mean_props[, c("Feature.name", level1, level2, "delta_prop"), with = FALSE],
                  by = "Feature.name", all.x = TRUE)
    
    
    Stat$diff   <- "non_significant"
    Stat[Stat$padj <= 0.05 & Stat$delta_prop >= delta_threshold, ]$diff <- "significant"
    Stat$Change <- "no-change"
    Stat[Stat$Estimate < 0, ]$Change <- "decrease"
    Stat[Stat$Estimate > 0, ]$Change <- "increase"
    Stat$diff            <- as.factor(Stat$diff)
    Stat$Change          <- as.factor(Stat$Change)
    Stat$Degradation_Tag <- tag
    
    Count_tag <- as.data.table(Stat)[, .(number = .N),
                                     by = .(comparaison, Condition, diff, Change)]
    Count_tag[, total   := sum(number), by = .(comparaison, Condition)]
    Count_tag[, Percent := 100 * number / total]
    Count_tag[, Degradation_Tag := tag]
    
    Stat_sig <- Stat[diff == "significant"]
    
    safe_tag <- gsub("'", "", tag)
    write.table(Stat_sig,
                paste0("diff_frag_", safe_tag, "_", comp_name, "_Stat.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    cat("Significant genes for", tag, ":", nrow(Stat_sig), "\n")
    
    list(Stat = as.data.table(Stat), Stat_sig = Stat_sig, Count = as.data.table(Count_tag))
  })
  
  names(results_list) <- actual_tags
  
  # --- 3.7 Combined tables ---
  all_sig    <- rbindlist(lapply(results_list, function(x) x$Stat_sig), fill = TRUE)
  all_counts <- rbindlist(lapply(results_list, function(x) x$Count),    fill = TRUE)
  
  cat("\nSignificant genes per tag:\n")
  print(all_sig[, .N, by = Degradation_Tag]) #!!
  
  write.table(all_sig,
              paste0("diff_frag_", comp_name, ".txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  write.table(all_sig[, .N, by = Degradation_Tag], paste0("diff_frag_Nb_", comp_name, ".txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  

# ==============================================================================
# 9. UPSETR
# ==============================================================================

# # Extract gene lists 
# # Extract both directions per tag
# intact_increase   <- results_list[["Intact"]]$Stat_sig[Change == "increase", Feature.name]
# intact_decrease   <- results_list[["Intact"]]$Stat_sig[Change == "decrease", Feature.name]
# 
# trunc5_increase   <- results_list[["5'_truncated"]]$Stat_sig[Change == "increase", Feature.name]
# trunc5_decrease   <- results_list[["5'_truncated"]]$Stat_sig[Change == "decrease", Feature.name]
# 
# trunc3_increase   <- results_list[["3'_truncated"]]$Stat_sig[Change == "increase", Feature.name]
# trunc3_decrease   <- results_list[["3'_truncated"]]$Stat_sig[Change == "decrease", Feature.name]
# 
# both_increase     <- results_list[["Both_truncated"]]$Stat_sig[Change == "increase", Feature.name]
# both_decrease     <- results_list[["Both_truncated"]]$Stat_sig[Change == "decrease", Feature.name]
# 
# 
# dexseq_frag <- fread("sig_rows.csv", header=TRUE)
# MUG <- fread("provided/134_genes_MUG.tsv", header=TRUE, sep="\t")
# LUG <- fread("provided/33_genes_LUG.tsv", header=TRUE, sep="\t")
# GMUCT_DOWN <- fread("Genelists/GMUCT_DOWN.txt", header=F, sep="\t")
# HIGH_CONFIDENCE_TARGETS <- fread("Genelists/HIGH_CONFIDENCE_TARGETS.txt", header=F, sep="\t")
# NAGARAJAN_DNE1_TARGETS <- fread("Genelists/NAGARAJAN_DNE1_TARGETS.txt", header=F, sep="\t")
# TRIBE_D153N <- fread("Genelists/TRIBE_D153N.txt", header=F, sep="\t")
# TRIBE_DNE1 <- fread("Genelists/TRIBE_DNE1.txt", header=F, sep="\t")
# dne1_ids <- c("AT1G13245", "AT1G78080", "AT2G18160", "AT2G22430", "AT2G47400", 
#               "AT3G22380", "AT4G34138", "AT5G19120", "AT4G64260", "AT2G43020")
# 
# gene_sets <- list(
#   "Intact"     = union(strip_isoform(intact_increase),  strip_isoform(intact_decrease)),
#   "5p_trunc"   = union(strip_isoform(trunc5_increase),  strip_isoform(trunc5_decrease)),
#   "3p_trunc"   = union(strip_isoform(trunc3_increase),  strip_isoform(trunc3_decrease)),
#   "Both_trunc" = union(strip_isoform(both_increase),    strip_isoform(both_decrease)),
#   "dexseq_frag"             = strip_isoform(dexseq_frag$groupID),  # your existing lists
#   "MUG"                     = strip_isoform(MUG$Feature.name),
#   "LUG"                     = strip_isoform(LUG$Feature.name),
#   "GMUCT_DOWN"              = GMUCT_DOWN$V1,
#   "HIGH_CONFIDENCE" = HIGH_CONFIDENCE_TARGETS$V1,
#   "NAGARAJAN_DNE1"  = NAGARAJAN_DNE1_TARGETS$V1,
#   "TRIBE_D153N"             = TRIBE_D153N$V1,
#   "TRIBE_DNE1"              = TRIBE_DNE1$V1,
#   "dne1_ids"                = dne1_ids
# )
# 
# 
# elements    <- unique(unlist(gene_sets))
# upset_input <- fromList(gene_sets)
# rownames(upset_input) <- elements
# 
# png(paste0("upset_gene_lists_", comp_name, ".png"),
#     width = 1800, height = 1050, res = 150)
# print(UpSetR::upset(UpSetR::fromList(gene_sets),
#                     nsets      = length(gene_sets),
#                     order.by   = "freq",
#                     text.scale = 1.3
#                     )) # keep.order = TRUE
# 
# dev.off()
# 
# # Intersection table
# all_combos <- unique(upset_input)
# results_upset <- data.frame()
# for (i in 1:nrow(all_combos)) {
#   active   <- names(all_combos)[all_combos[i, ] == 1]
#   if (length(active) == 0) next
#   is_match <- apply(upset_input, 1, function(x) all(x == all_combos[i, ]))
#   genes    <- rownames(upset_input)[is_match]
#   n        <- length(genes)
#   if (n > 0) {
#     results_upset <- rbind(results_upset, data.frame(
#       count = n,
#       sets  = paste(active, collapse = " & "),
#       genes = paste(genes, collapse = ", ")
#     ))
#   }
# }
# 
# cat("\nTop intersections:\n")
# print(head(results_upset[order(-results_upset$count), ], 50))
# 
# 
# write.table(results_upset[order(-results_upset$count), ],
#            paste0("results_upset_", comp_name, ".txt"),
#            sep = "\t", row.names = FALSE, quote = FALSE)
# 
# 
# # Return everything for downstream access
# invisible(list(
#   comp_name    = comp_name,
#   results_list = results_list,
#   all_sig      = all_sig,
#   all_counts   = all_counts,
#   upset_input  = upset_input,
#   intersections = results_upset,
#   gene_sets     = gene_sets  
# )) 

multi_hit <- all_sig[, .N, by = Feature.name][N >= 2, Feature.name]
heat_data <- all_sig[Feature.name %in% multi_hit]
heat_data[, signed_est := ifelse(Change == "increase", abs(Estimate), -abs(Estimate))]

p_heat <- ggplot(heat_data, aes(x = Degradation_Tag, y = Feature.name, fill = signed_est)) +
  geom_tile() +
  scale_fill_gradient2(low = "#56B4E9", mid = "white", high = "#E69F00", midpoint = 0) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 6)) +
  labs(x = "Fragment type", y = "Gene", fill = "Log-odds")

}  
# ==============================================================================
# 4. RUN BOTH COMPARISONS INDEPENDENTLY
# ==============================================================================

comp_Col0_xrn4      <- run_comparison("Col0",  "xrn4",      SampleTable, ds, actual_tags)
comp_xrn4_xrn4_dne1 <- run_comparison("xrn4",  "xrn4_dne1", SampleTable, ds, actual_tags)


