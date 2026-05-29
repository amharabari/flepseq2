library(stageR)
library(DRIMSeq)

# ==============================================================================
# 1. INITIALIZATION & METADATA
# ==============================================================================
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
# 2. HELPER FUNCTIONS
# ==============================================================================

strip_isoform <- function(x) sub("\\.[0-9]+$", "", x)
no_na         <- function(x) ifelse(is.na(x), 1, x)

# ==============================================================================
# 3. STATIC GENE LISTS
# ==============================================================================

# dexseq_frag             <- fread("flepseq2/sig_rows.csv", header = TRUE)
# MUG                     <- fread("flepseq2/provided/134_genes_MUG.tsv",            header = TRUE,  sep = "\t")
# LUG                     <- fread("flepseq2/provided/33_genes_LUG.tsv",             header = TRUE,  sep = "\t")
# GMUCT_DOWN              <- fread("flepseq2/Genelists/GMUCT_DOWN.txt",              header = FALSE, sep = "\t")
# HIGH_CONFIDENCE_TARGETS <- fread("flepseq2/Genelists/HIGH_CONFIDENCE_TARGETS.txt", header = FALSE, sep = "\t")
# NAGARAJAN_DNE1_TARGETS  <- fread("flepseq2/Genelists/NAGARAJAN_DNE1_TARGETS.txt",  header = FALSE, sep = "\t")
# TRIBE_D153N             <- fread("flepseq2/Genelists/TRIBE_D153N.txt",             header = FALSE, sep = "\t")
# TRIBE_DNE1              <- fread("flepseq2/Genelists/TRIBE_DNE1.txt",              header = FALSE, sep = "\t")
# dne1_ids                <- c("AT1G13245", "AT1G78080", "AT2G18160", "AT2G22430", "AT2G47400",
#                              "AT3G22380", "AT4G34138", "AT5G19120", "AT4G64260", "AT2G43020")
# 
# static_gene_sets <- list(
#   dexseq_frag     = strip_isoform(dexseq_frag$groupID),
#   MUG             = strip_isoform(MUG$Feature.name),
#   LUG             = strip_isoform(LUG$Feature.name),
#   GMUCT_DOWN      = GMUCT_DOWN$V1,
#   HIGH_CONFIDENCE = HIGH_CONFIDENCE_TARGETS$V1,
#   NAGARAJAN_DNE1  = NAGARAJAN_DNE1_TARGETS$V1,
#   TRIBE_D153N     = TRIBE_D153N$V1,
#   TRIBE_DNE1      = TRIBE_DNE1$V1,
#   dne1_ids        = dne1_ids
# )

# ==============================================================================
# 4. COMPARISON 1: Col0 vs xrn4
# ==============================================================================
level1    <- "Col0"
level2    <- "xrn4"
comp_name <- "Col0_xrn4"

target_barcodes <- SampleTable[Genotype %in% c(level1, level2) & 
                                 Condition == "Seedling", Barcode]

global_counts <- ds %>%
  dplyr::select(Barcode, Feature.name, Degradation_Tag) %>%
  collect() %>% as.data.table()
global_counts <- global_counts[Barcode %chin% target_barcodes]
global_counts <- global_counts[!grepl("^ATCG|^ATMG", Feature.name)]
global_counts <- global_counts[Degradation_Tag != "No_CDS"]
global_counts <- global_counts[, .(raw_count = .N), 
                               by = .(Barcode, Feature.name, Degradation_Tag)]
global_counts <- merge(global_counts,
                       SampleTable[, .(Barcode, Sample_name, Genotype, Rep, libsize)],
                       by = "Barcode", all.x = TRUE)

gene_totals  <- global_counts[, .(Total = sum(raw_count)), 
                              by = .(Feature.name, Sample_name, Genotype, Rep)]
list_targets <- gene_totals[, .(passed = sum(Total >= 150) == 3),
                            by = .(Feature.name, Genotype)][passed == TRUE, 
                                                            unique(Feature.name)]
cat("Genes passing flexible filter:", length(list_targets), "\n")

final_stats <- global_counts[Feature.name %in% list_targets]
final_stats[, gene_raw_total := sum(raw_count), by = .(Sample_name, Feature.name)]
final_stats[, prop_tag       := raw_count / pmax(gene_raw_total, 1)]
final_stats[, cpm            := (raw_count / libsize) * 1e6]
final_stats[, Condition      := "Seedling"]

# 1) Summarize counts per sample (Sample_name) per gene per tag
counts_long <- final_stats[,
                           .(raw_count = sum(raw_count)),
                           by = .(Feature.name, Degradation_Tag, Sample_name)
]

# 2) Create gene_id and feature_id
counts_long[, gene_id    := Feature.name]
counts_long[, feature_id := paste0(Feature.name, ":", Degradation_Tag)]

counts_wide <- dcast(
  counts_long,                              # ← use counts_long, NOT final_stats
  gene_id + feature_id ~ Sample_name,
  value.var = "raw_count",
  fill = 0
)

setcolorder(counts_wide, c("gene_id", "feature_id"))

samps <- data.frame(
  sample_id = SampleTable[Genotype %in% c(level1, level2) & Condition == "Seedling", Sample_name],
  group     = SampleTable[Genotype %in% c(level1, level2) & Condition == "Seedling", Genotype]
)
samps$group <- relevel(factor(samps$group), ref = level1)
samps <- samps[match(colnames(counts_wide)[-(1:2)], samps$sample_id), ]
stopifnot(all(colnames(counts_wide)[-(1:2)] == samps$sample_id))

# ← batch defined HERE, before dmDSdata
samps$batch <- factor(ifelse(samps$sample_id %in% c("Col0_Rep2", "xrn4_Rep3"), "B", "A"))
print(table(samps$group, samps$batch))

n       <- nrow(samps)
n.small <- min(table(samps$group))

d <- dmDSdata(counts = as.data.frame(counts_wide), samples = samps)
d <- dmFilter(d, min_samps_feature_expr = n.small, min_feature_expr = 10,
              min_samps_feature_prop = n.small, min_feature_prop = 0.05,
              min_samps_gene_expr = n, min_gene_expr = 150)
cat("Genes after DRIMSeq filter:", length(levels(factor(counts(d)$gene_id))), "\n")

# ← single design, batch-corrected
design <- model.matrix(~ batch + group, data = DRIMSeq::samples(d))
d <- dmPrecision(d, design = design)
d <- dmFit(d,      design = design)
d <- dmTest(d,     coef   = paste0("group", level2))

# ← results extracted FIRST
res_gene <- as.data.table(results(d));                     res_gene[, pvalue := no_na(pvalue)]
res_feat <- as.data.table(results(d, level = "feature"));  res_feat[, pvalue := no_na(pvalue)]

# ← diagnostics AFTER
hist(res_gene$pvalue, breaks = 50, main = "Gene-level raw p-values (batch-corrected)")
hist(res_feat$pvalue, breaks = 50, main = "Feature-level raw p-values (batch-corrected)")
summary(res_gene$pvalue)
cat("Genes with p < 0.05:", sum(res_gene$pvalue < 0.05, na.rm = TRUE), "\n")
# Critical check — features per gene after filtering
filtered_counts <- counts(d)
print(table(table(filtered_counts$gene_id)))


res_gene <- as.data.table(results(d));        res_gene[, pvalue := no_na(pvalue)]
res_feat <- as.data.table(results(d, level = "feature")); res_feat[, pvalue := no_na(pvalue)]

res_gene[, padj_gene := p.adjust(pvalue, "BH")]
res_gene[, padj_gene := p.adjust(pvalue, "BH")]

feat_annotated <- merge(
  res_feat[, .(gene_id, feature_id, pvalue_feat = pvalue)],
  res_gene[, .(gene_id, pvalue_gene = pvalue, padj_gene)],
  by = "gene_id"
)
feat_annotated[, Degradation_Tag := sub("^.*:", "", feature_id)]
feat_annotated <- merge(feat_annotated,
                        mean_props[, .(gene_id, feature_id, mean_prop_level1, mean_prop_level2, delta_prop, Change)],
                        by = c("gene_id", "feature_id"), all.x = TRUE)

# FDR sweep — gene-level BH + raw feature p < 0.05
for (fdr in c(0.05, 0.10, 0.20)) {
  n_sig <- feat_annotated[padj_gene <= fdr & pvalue_feat <= 0.05, uniqueN(gene_id)]
  cat("FDR", fdr, "→", n_sig, "genes\n")
}

Stat_sig_C0xrn4 <- feat_annotated[padj_gene <= 0.10 & pvalue_feat <= 0.05]
cat("Significant gene×tag combinations:", nrow(Stat_sig_C0xrn4), "\n")
print(Stat_sig_C0xrn4[, .N, by = Degradation_Tag])
print(Stat_sig_C0xrn4[, .N, by = .(Degradation_Tag, Change)])


# ==============================================================================
# 5. COMPARISON 2: xrn4 vs xrn4_dne1
level1    <- "xrn4"
level2    <- "xrn4_dne1"
comp_name <- "xrn4_xrn4_dne1"

target_barcodes <- SampleTable[Genotype %in% c(level1, level2) & Condition == "Seedling", Barcode]

global_counts <- ds %>%
  dplyr::select(Barcode, Feature.name, Degradation_Tag) %>%
  collect() %>% as.data.table()
global_counts <- global_counts[Barcode %chin% target_barcodes]
global_counts <- global_counts[!grepl("^ATCG|^ATMG", Feature.name)]
global_counts <- global_counts[Degradation_Tag != "No_CDS"]
global_counts <- global_counts[, .(raw_count = .N), by = .(Barcode, Feature.name, Degradation_Tag)]
global_counts <- merge(global_counts,
                       SampleTable[, .(Barcode, Sample_name, Genotype, Rep, libsize)],
                       by = "Barcode", all.x = TRUE)

gene_totals  <- global_counts[, .(Total = sum(raw_count)), by = .(Feature.name, Sample_name, Genotype, Rep)]
list_targets <- gene_totals[, .(passed = sum(Total >= 150) == 3),
                            by = .(Feature.name, Genotype)][passed == TRUE, unique(Feature.name)]
cat("Genes passing flexible filter:", length(list_targets), "\n")

final_stats <- global_counts[Feature.name %in% list_targets]
final_stats[, gene_raw_total := sum(raw_count), by = .(Sample_name, Feature.name)]
final_stats[, prop_tag       := raw_count / pmax(gene_raw_total, 1)]
final_stats[, cpm            := (raw_count / libsize) * 1e6]
final_stats[, Condition      := "Seedling"]

# 1) Summarize counts per sample (Sample_name) per gene per tag
counts_long <- final_stats[,
                           .(raw_count = sum(raw_count)),
                           by = .(Feature.name, Degradation_Tag, Sample_name)
]

# 2) Create gene_id and feature_id
counts_long[, gene_id    := Feature.name]
counts_long[, feature_id := paste0(Feature.name, ":", Degradation_Tag)]


counts_wide <- dcast(
  counts_long,
  gene_id + feature_id ~ Sample_name,
  value.var = "raw_count",
  fill = 0)


samps <- data.frame(
  sample_id = SampleTable[Genotype %in% c(level1, level2) & Condition == "Seedling", Sample_name],
  group     = SampleTable[Genotype %in% c(level1, level2) & Condition == "Seedling", Genotype]
)
samps$group <- relevel(factor(samps$group), ref = level1)
samps <- samps[match(colnames(counts_wide)[-(1:2)], samps$sample_id), ]
stopifnot(all(colnames(counts_wide)[-(1:2)] == samps$sample_id))
n           <- nrow(samps)
n.small     <- min(table(samps$group))

d <- dmDSdata(counts = as.data.frame(counts_wide), samples = samps)
d <- dmFilter(d, 
              min_samps_feature_expr = n.small, 
              min_feature_expr = 5,
              min_samps_feature_prop = 1, 
              min_feature_prop = 0.05,
              min_samps_gene_expr = n.small, 
              min_gene_expr = 50)

cat("Genes after DRIMSeq filter:", length(levels(factor(counts(d)$gene_id))), "\n")

filtered_counts <- counts(d)
table(table(filtered_counts$gene_id))  # distribution of feature counts per gene

design <- model.matrix(~ group, data = DRIMSeq::samples(d))
d <- dmPrecision(d, design = design)
d <- dmFit(d,      design = design)
d <- dmTest(d,     coef   = paste0("group", level2))

res_gene <- as.data.table(results(d));                    res_gene[, pvalue := no_na(pvalue)]
res_feat <- as.data.table(results(d, level = "feature")); res_feat[, pvalue := no_na(pvalue)]

pScreen            <- res_gene$pvalue    
names(pScreen)     <- res_gene$gene_id
pConfirm           <- matrix(res_feat$pvalue, ncol = 1)
rownames(pConfirm) <- res_feat$feature_id
tx2gene <- as.data.frame(res_feat[, .(feature_id, gene_id)])

stageRObj <- stageRTx(pScreen, pConfirm, tx2gene = as.data.frame(res_feat[, .(feature_id, gene_id)]))
stageRObj <- stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)

adj <- as.data.table(
  getAdjustedPValues(stageRObj,
                     order = FALSE,
                     onlySignificantGenes = FALSE)
)
# adj has columns geneID (gene-level padj), txID (feature-level padj)
setnames(adj, c("geneID", "txID", "gene", "transcript"),
         c("gene_id", "feature_id", "p_gene", "p_transcript"))
adj[, gene_id    := as.character(gene_id)]
adj[, feature_id := as.character(feature_id)]


feat_full <- merge(
  res_feat[, .(gene_id, feature_id, pvalue_raw = pvalue)],
  adj[,   .(gene_id, feature_id, p_gene, p_transcript)],
  by = c("gene_id", "feature_id"),
  all.x = TRUE
)

feat_full[, Degradation_Tag := sub("^.*:", "", feature_id)]
mean_props <- final_stats[
  Feature.name %in% unique(feat_full$gene_id),
  .(mean_prop = mean(prop_tag, na.rm = TRUE)),
  by = .(Feature.name, Degradation_Tag, Genotype)
]

mean_props <- dcast(
  mean_props,
  Feature.name + Degradation_Tag ~ Genotype,
  value.var = "mean_prop"
)
setnames(mean_props, c("Feature.name", "Degradation_Tag"),
         c("gene_id", "Degradation_Tag"))

mean_props[, feature_id := paste0(gene_id, ":", Degradation_Tag)]
setnames(mean_props, c(level1, level2),
         c("mean_prop_level1", "mean_prop_level2"))
mean_props[, delta_prop := abs(mean_prop_level1 - mean_prop_level2)]
mean_props[, Change := fcase(
  mean_prop_level2 >  mean_prop_level1, "increase",
  mean_prop_level2 <  mean_prop_level1, "decrease",
  default = "no-change"
)]

drim_padj_xrn4dne1 <- merge(
  feat_full,
  mean_props[, .(gene_id, feature_id,
                 mean_prop_level1, mean_prop_level2,
                 delta_prop, Change)],
  by = c("gene_id", "feature_id"),
  all.x = TRUE
)

drim_padj_xrn4dne1[,
                   diff := ifelse(!is.na(p_gene)       & p_gene       <= 0.05 &
                                    !is.na(p_transcript) & p_transcript <= 0.05 &
                                    !is.na(delta_prop)   & delta_prop   >= 0.05,
                                  "significant", "non_significant")
]


Stat_sig_xrn4dne1 <- drim_padj_xrn4dne1[diff == "significant"]
cat("Significant gene×tag combinations:", nrow(Stat_sig_xrn4dne1), "\n")
print(Stat_sig_xrn4dne1[, .N, by = Degradation_Tag])

for (tag in actual_tags) {
  safe_tag <- gsub("'", "", tag)
#  write.table(Stat_sig_xrn4dne1[Degradation_Tag == tag],
#              paste0("flepseq2/fragment_analysis/diff_frag_", safe_tag, "_", comp_name, "_Stat.txt"),
#              sep = "\t", row.names = FALSE, quote = FALSE)
}
#write.table(Stat_sig_xrn4dne1,
#            paste0("flepseq2/fragment_analysis/drim_diff_frag_", comp_name, ".txt"),
 #           sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(Stat_sig_xrn4dne1[, .N, by = Degradation_Tag],
#            paste0("flepseq2/fragment_analysis/drim_diff_frag_Nb_", comp_name, ".txt"),
#            sep = "\t", row.names = FALSE, quote = FALSE)

dynamic_gene_sets_xrn4dne1 <- setNames(
  lapply(actual_tags, function(tag) strip_isoform(Stat_sig_xrn4dne1[Degradation_Tag == tag, gene_id])),
  gsub("'", "p", gsub("_truncated", "_trunc", actual_tags))
)
gene_sets_xrn4dne1   <- c(dynamic_gene_sets_xrn4dne1, static_gene_sets)
upset_input_xrn4dne1 <- fromList(gene_sets_xrn4dne1)
rownames(upset_input_xrn4dne1) <- unique(unlist(gene_sets_xrn4dne1))

#png(paste0("flepseq2/fragment_analysis/drim_upset_gene_lists_", comp_name, ".png"), width = 1800, height = 1050, res = 150)
#print(UpSetR::upset(UpSetR::fromList(gene_sets_xrn4dne1), nsets = length(gene_sets_xrn4dne1), order.by = "freq", text.scale = 1.3))
#dev.off()

results_upset_xrn4dne1 <- rbindlist(lapply(seq_len(nrow(unique(upset_input_xrn4dne1))), function(i) {
  combo  <- unique(upset_input_xrn4dne1)[i, ]
  active <- names(combo)[combo == 1]
  if (length(active) == 0) return(NULL)
  genes  <- rownames(upset_input_xrn4dne1)[apply(upset_input_xrn4dne1, 1, function(x) all(x == combo))]
  if (length(genes) == 0) return(NULL)
  data.table(count = length(genes), sets = paste(active, collapse = " & "), genes = paste(genes, collapse = ", "))
}))[order(-count)]
cat("\nTop intersections:\n"); print(head(results_upset_xrn4dne1, 50))
#write.table(results_upset_xrn4dne1,
#            paste0("flepseq2/fragment_analysis/results_upset_", comp_name, ".txt"),
#            sep = "\t", row.names = FALSE, quote = FALSE)