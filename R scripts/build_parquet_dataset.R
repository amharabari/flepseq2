# Converts raw FLEP-seq2 bamparsed TSV.GZ files into a Parquet dataset,
# selecting relevant columns and tagging each read with its barcode.

library(data.table)
library(arrow)

raw_path    <- "flepseq2/bamparsed/"
parquet_dir <- "flepseq2/parquet_dataset"

cols_to_keep <- c(
  "Read.name", "Read.start", "Read.end",
  "Feature.name", "Feature.start", "Feature.end", "Feature.strand",
  "CDS.start_codon", "CDS.stop_codon", "CDS_dist_5prime", "CDS_dist_3prime",
  "Degradation_Tag", "RA5_tag", "PolyA", "Add_tail", "PolyA_length_basecall",
  "PolyA_length_signal", "Add_tail_length", "Dist_from_5prime",
  "Dist_from_3prime", "Tail_tag")

if (!dir.exists(parquet_dir)) dir.create(parquet_dir, recursive = TRUE)

convert_one_file <- function(f, out_dir, cols) {
  bc_name  <- regmatches(f, regexpr("barcode\\d+", f))
  out_path <- file.path(out_dir, paste0(bc_name, ".parquet"))
  
  if (file.exists(out_path)) {
    cat(sprintf("Skipping %s (already converted)\n", bc_name))
    return(invisible(NULL))
  }
  
  t1 <- Sys.time()
  dt <- fread(f, sep = "\t", select = cols, showProgress = FALSE)
  dt[, Barcode := bc_name]
  write_parquet(dt, out_path, compression = "snappy")
  t2 <- Sys.time()
  
  cat(sprintf("[%s] converted | rows: %s | time: %.1fs\n",
              bc_name, format(nrow(dt), big.mark = ","),
              as.numeric(t2 - t1)))
  
  rm(dt); gc(verbose = FALSE)
  invisible(NULL)
}

# ---- Run conversion for all files ----------------------------------------
file_paths <- list.files(raw_path, pattern = "\\.gz$", full.names = TRUE)

cat("--- Converting", length(file_paths), "files to Parquet ---\n")
invisible(lapply(file_paths, convert_one_file, out_dir = parquet_dir, cols = cols_to_keep))

ds <- open_dataset(parquet_dir)
cat("Dataset ready. Columns:", paste(names(schema(ds)), collapse = ", "), "\n")
cat("Total rows:", format(nrow(ds), big.mark = ","), "\n")