library (grid)
library(tidyverse)
library(data.table)
library(ggplot2)

# 1. Classifications
# 2. 
# 3. 


setwd("/Users/hyu/Documents/FlepSeq/xrn4/Plots_frag")
Tab <- fread("/Users/hyu/Documents/FlepSeq/xrn4/DataTables/RUN_Remy_seedlings_xrn4_cleaned.txt")

# add read_start and read_end
# Tab <- Tab %>% separate(read_core_id, into = c("ID", "Chr", "read_start", "read_end"), sep = ",")

Tab$read_start <- sub(".*,.*,(.*),.*", "\\1", Tab$read_core_id)
Tab$read_end <- sub(".*,.*,(.*),(.*)", "\\2", Tab$read_core_id)

test <- slice(Tab, 1:3000)
test
rm(test)



# ==============================================================================
# 1. Classifications 
# ==============================================================================


## (1a) Adding of cds end position on AGIs

annotation <- fread("/Users/hyu/Documents/FlepSeq/analysis2408/R7 frag exon/Araport11_GFF3_genes_transposons.current.gff", encoding="UTF-8")
colnames(annotation) <- c("Chr", "source", "feature", "cds_start", "cds_end", "?", "orientation", "??", "info")
setDT(annotation)[orientation == "-", cds_end := cds_start]

annotation[, AGI := gsub(".*ID=(AT[1-5]G\\d{5}).*","\\1", info)]

annotation1 <- annotation %>%
  filter(feature == "protein") %>%
  select(AGI, cds_end, orientation)

annotation1a <- filter(annotation1, orientation == "+")
annotation1b <- filter(annotation1, orientation == "-")

annotation2a <- annotation1a %>%
  group_by(AGI) %>%
  arrange(desc(cds_end)) %>%
  slice(1)

annotation2b <- annotation1b %>%
  group_by(AGI) %>%
  arrange(cds_end) %>%
  slice(1)

annotation_cds_end <- bind_rows(annotation2a, annotation2b)
annotation_cds_end <- arrange(annotation_cds_end, AGI)

setDT(annotation_cds_end)
Tab_annotat_pre <- annotation_cds_end[Tab, on = "AGI", nomatch = 0]

head(Tab_annotat_pre)

rm(list=setdiff(ls(), "Tab_annotat_pre"))
#

## (1b) Adding of cds start position on AGIs

annotation <- fread("/Users/hyu/Documents/FlepSeq/analysis2408/R7 frag exon/Araport11_GFF3_genes_transposons.current.gff", encoding="UTF-8")
colnames(annotation) <- c("Chr", "source", "feature", "cds_start", "cds_end", "?", "orientation", "??", "info")
setDT(annotation)[orientation == "-", cds_start := cds_end]

annotation[, AGI := gsub(".*ID=(AT[1-5]G\\d{5}).*","\\1", info)]

annotation1 <- annotation %>%
  filter(feature == "protein") %>%
  select(AGI, cds_start, orientation)

annotation1a <- filter(annotation1, orientation == "+")
annotation1b <- filter(annotation1, orientation == "-")

annotation2a <- annotation1a %>%
  group_by(AGI) %>%
  arrange(cds_start) %>%
  slice(1)

annotation2b <- annotation1b %>%
  group_by(AGI) %>%
  arrange(desc(cds_start)) %>%
  slice(1)

annotation_cds_start <- bind_rows(annotation2a, annotation2b)
annotation_cds_start <- arrange(annotation_cds_start, AGI)
annotation_cds_start <- select(annotation_cds_start, AGI, cds_start)

setDT(annotation_cds_start)
Tab_annotat <- annotation_cds_start[Tab_annotat_pre, on = "AGI", nomatch = 0]

head(Tab_annotat)

# calculate 5dist_cds_start & 3dist_cds_end

rm(list=setdiff(ls(), "Tab_annotat"))

setwd("/Users/hyu/Documents/FlepSeq/xrn4/DataTables")
fwrite(Tab_annotat, "RUN_Remy_seedlings_xrn4_cleaned_frag_pre.csv")

Tab_annotat <- fread("/Users/hyu/Documents/FlepSeq/xrn4/DataTables/RUN_Remy_seedlings_xrn4_cleaned_frag_pre.csv")

rm(Tab_annotat)

prot_annotat_a <- filter(Tab_annotat, orientation == "+")
prot_annotat_a$five_dist_cds <- as.numeric(prot_annotat_a$cds_start) - as.numeric(prot_annotat_a$read_start)
prot_annotat_a$three_dist_cds <- as.numeric(prot_annotat_a$read_end) - as.numeric(prot_annotat_a$cds_end)

setwd("/Users/hyu/Documents/FlepSeq/xrn4/DataTables")
fwrite(prot_annotat_a, "prot_annotat_a.csv")

prot_annotat_b <- filter(Tab_annotat, orientation == "-")
prot_annotat_b$five_dist_cds <- as.numeric(prot_annotat_b$read_end) - as.numeric(prot_annotat_b$cds_start)
prot_annotat_b$three_dist_cds <- as.numeric(prot_annotat_b$cds_end) - as.numeric(prot_annotat_b$read_start)

#test <- slice(prot_annotat_b, 1:3000)
#test
#rm(test)

rm(Tab_annotat)

Tab_annotat2 <- bind_rows(prot_annotat_a, prot_annotat_b)



library(data.table)
prot_annotat_a <- as.data.table(prot_annotat_a)
prot_annotat_b <- as.data.table(prot_annotat_b)
Tab_annotat2 <- rbindlist(list(prot_annotat_a, prot_annotat_b))


setwd("/Users/hyu/Documents/FlepSeq/xrn4/DataTables")
fwrite(Tab_annotat2, "RUN_Remy_seedlings_xrn4_cleaned_frag_pre.csv")


## (1c) Reads classification

Tab <- fread("/Users/hyu/Documents/FlepSeq/xrn4/DataTables/RUN_Remy_seedlings_xrn4_cleaned_frag_pre.csv")

#head(Tab)

Tab <- as.data.table(Tab)

Tab$info_cds <- "intact"
Tab$info_polyA <- ">60"

Tab[as.numeric(three_dist_cds) <= 0, info_cds := "Three_in_cds"] 
Tab[as.numeric(five_dist_cds) <= 0, info_cds := "Five_in_cds"] 
Tab[as.numeric(three_dist_cds) <= 0 & as.numeric(five_dist_cds) <= 0, info_cds := "Both_in_cds"] 


Tab[polya_round <= 10, info_polyA := "<=10"] 
Tab[polya_round > 10 & polya_round <= 30, info_polyA := ">10_30"] 
Tab[polya_round > 30 & polya_round <= 60, info_polyA := ">30_60"] 

test <- slice(Tab, 1:5000)
view(test)

rm(test)

fwrite(Tab, "RUN_Remy_seedlings_xrn4_cleaned_frag.txt")

Tab <- fread("/Users/hyu/Documents/FlepSeq/xrn4/DataTables/RUN_Remy_seedlings_xrn4_cleaned_frag.txt")

