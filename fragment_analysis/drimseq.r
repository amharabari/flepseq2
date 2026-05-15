library(DRIMSeq)

# Reshape from your existing final_stats data.table
# final_stats has: Feature.name, Degradation_Tag, Barcode, raw_count, Sample_name

counts_wide <- dcast(
  final_stats,
  Feature.name + Degradation_Tag ~ Sample_name,
  value.var = "raw_count",
  fill = 0
)

# Rename to match DRIMSeq expectations
setnames(counts_wide,
         old = c("Feature.name", "Degradation_Tag"),
         new = c("gene_id", "feature_id"))

# Sample table for DRIMSeq
samps <- data.frame(
  sample_id = SampleTable[Genotype %in% c("Col0","xrn4") & Condition == "Seedling", Sample_name],
  group     = SampleTable[Genotype %in% c("Col0","xrn4") & Condition == "Seedling", Genotype]
)

# Build DRIMSeq object
d <- dmDSdata(counts = as.data.frame(counts_wide), samples = samps)

# Filter (equivalent to your ≥150 reads filter)
n       <- nrow(samps)
n.small <- min(table(samps$group))

d <- dmFilter(d,
              min_samps_feature_expr  = n.small,
              min_feature_expr        = 10,
              min_samps_feature_prop  = n.small,
              min_feature_prop        = 0.05,
              min_samps_gene_expr     = n,
              min_gene_expr           = 150)   # your existing threshold

# Design matrix
design <- model.matrix(~ group, data = DRIMSeq::samples(d))

# Fit
d <- dmPrecision(d, design = design)
d <- dmFit(d,      design = design)
d <- dmTest(d,     coef   = "groupxrn4")

# Results
res_gene <- results(d)          # one p-value per gene (any tag changed?)
res_feat <- results(d, level = "feature")  # one p-value per gene × tag

library(stageR)

pScreen    <- res_gene$pvalue
names(pScreen) <- res_gene$gene_id

pConfirm   <- matrix(res_feat$pvalue, ncol = 1)
rownames(pConfirm) <- res_feat$feature_id

stageRObj  <- stageRTx(pScreen, pConfirm, allowNA = TRUE, tx2gene = res_feat[, c("feature_id","gene_id")])
stageRObj  <- stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)
drim.padj  <- getAdjustedPValues(stageRObj, order = TRUE, onlySignificantGenes = FALSE)