# Adapting Heike's scripts for gene by gene PolyA analysis
#rm(list=ls())
gc()

library(data.table)
library(stringr)
library(arrow)   # For Out-of-Memory processing
library(dplyr)   # For Arrow "remote control"
library(ggplot2)
library(tidyr)
library(ggplot2)
library(scales)
library(viridis)
library(viridisLite)
library(eulerr)
library(grid)
library(patchwork)
library(forcats)
setDTthreads(0)

genotype_colors <- c(
  "Col0"      = "#696969", 
  "xrn4"      = "#FF1493", 
  "xrn4_dne1" = "#00BFFF")

#ticks and breaks for polyA profiles
custom_breaks <- c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200)
every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}

# Imporiting the dataset
ds_full <- open_dataset("parquet_dataset")
SampleTable <- fread("bamparsed/barcode_correspondance.tsv")
SampleTable <- SampleTable[, .(Barcode, Sample_name, Genotype, Condition, Rep)]
# Identify Target Samples
target_barcodes <- SampleTable[
  Genotype %in% c("Col0", "xrn4", "xrn4_dne1") & Condition == "Seedling", Barcode]

gc()

# Extract and Filter Data using Arrow Pipeline
ds_full %>%
  group_by(Barcode) %>%
  write_dataset("parquet_partitioned", format = "parquet")

rm(ds_full); gc()

ds_part <- open_dataset("parquet_partitioned")

df_targets <- ds_part %>%
  filter(Barcode %in% target_barcodes) %>%
  select(Barcode, Feature.name, PolyA, Add_tail, PolyA_length_signal) %>%
  collect() %>%
  setDT()
gc()

df_targets <- df_targets[!grepl("^ATCG|^ATMG", Feature.name)]


setnames(df_targets,
         old = c("Feature.name", "PolyA", "PolyA_length_signal", "Add_tail"),
         new = c("AGI", "polytail", "polya_length", "additional_tail"))



df_targets[, polytail := fifelse(polytail == "", NA_character_, polytail)]
df_targets[, additional_tail := fifelse(additional_tail == "", NA_character_, additional_tail)]
df_targets[, additional_tail := fifelse(str_detect(polytail, "(T+)$"), 
                                 str_extract(polytail, "(T+)$"), 
                                 additional_tail)]

# count the letters in additional tail
df_targets[, `:=`(
  add_tail_A = stringi::stri_count_fixed(additional_tail, "A"),
  add_tail_T = stringi::stri_count_fixed(additional_tail, "T")
)]

# Calculate the percentage of A, T, C, and G
df_targets[, `:=`(
  add_tail_pct_A = add_tail_A / nchar(additional_tail) * 100,
  add_tail_pct_T = add_tail_T / nchar(additional_tail) * 100
)]

# tail sorting
df_targets[, tail := "other"]
df_targets[is.na(additional_tail),  tail := "A-Tails"]
df_targets[add_tail_pct_T > 74, tail := "uridylated A-Tails"]
df_targets[(is.na(polytail)), tail := "none"] 

# remove all 3' truncated fragments / untailed reads if you want to analyse tail size. 
df_targets <- df_targets[polya_size > 9] 
df_targets <- df_targets[!tail == "none"] 

# tails with mismatch in polytail / add. tail has a long A stretch and ends with A or U 
df_targets[tail =="other" & str_detect(additional_tail, "(?=.*A{10,})A$"), tail := "A-Tails"]
df_targets[tail =="other" & str_detect(additional_tail, "(?=.*A{10,})T$"), tail := "uridylated A-Tails"]

df_targets[, polya_size := floor(polya_length)]

# remove rows with NA or emty string in polya_size
df_targets <- df_targets[!is.na(polya_size)]
df_targets <- df_targets[!polytail ==""]

df_targets <- merge(SampleTable, df_targets, by = "Barcode", df_targets.y = TRUE)
setnames(df_targets, "Rep", "rep")
rm(SampleTable); gc()  


grouped_counts <- df_targets[, .N, by = .(AGI, Genotype, rep)]
filtered_counts <- grouped_counts[N >= 50]
unique_combinations <- uniqueN(grouped_counts, by = c("rep", "Genotype"))
list_C50 <- unique(filtered_counts[, if (.N == unique_combinations) .SD, by = AGI]$AGI)
print(length(list_C50))
rm(filtered_counts, grouped_counts, unique_combinations)

C50 <- df_targets[AGI %in% list_C50]
C50[, weight := 1/.N, by = .(AGI, rep, Genotype)]
rm(df_targets, list_C50); gc()

colnames(C50)

# intergenic
# "intergenic" using median across all AGIs, limits effects of very abundant AGIs. The mean of the median is calculated across reps for each polyA size. 
## the cross-joint table merge is the data.table equivalent to the Drop=FALSE; 
## It ensures that the non-existing values are not omitted and can be filled with 0.

all_combinations <- CJ( 
  Genotype = unique(C50$Genotype),
  rep = unique(C50$rep),
  polya_size = unique(C50$polya_size),
  AGI = unique(C50$AGI))

distri <- C50[ , .(n = .N), by = .(Genotype, polya_size, rep, AGI)] # total Nr of reads for each polyA size in each library by AGIs

distri <- merge(all_combinations, distri, by = c("AGI", "Genotype", "rep", "polya_size"), all.x = TRUE) 
distri[is.na(n), n := 0] # fills the not-detected values with 0
distri[, `:=` (total = sum(n), perc = 100*(n) / sum(n)), by = .(Genotype, AGI, rep)]

distri[, median := median(perc, na.rm = T),  by = .(Genotype, polya_size, rep)]
distri[, `:=` (mean = mean(median, na.rm = T), SD = sd(median, na.rm = T)), by = .(Genotype, polya_size)] # mean of the median across reps but not AGIs. The result is that high abundant AGis have still more impact than the low abundant ones
distri <- distri[, .(polya_size, Genotype, SD, mean)] # reduce to necessary columns, optional

distri <- unique(distri, by = c("polya_size", "Genotype", "SD", "mean")) # collaps to unique entries for each genotype and polya bin

rm(all_combinations)

# two-step distri invented by hanzhang. Same weight for each AGI
# The mean of the median is calculated across all AGI-rep-pairs for each polyA size. This gives each AGI-rep pair the same weight. A second mean operation is neede to caclulate the final distribution. This gives exactly the same result as the "weighted" calculation, but it takes more processing time.

all_combinations <- CJ(
  Genotype   = unique(C50$Genotype),
  rep        = unique(C50$rep),
  AGI        = unique(C50$AGI),
  polya_size = unique(C50$polya_size))

distribution_full <- C50[, .(n = .N), by = .(Genotype, rep, AGI, polya_size)]
distribution_full <- distribution_full[ all_combinations, on = .(Genotype, rep, AGI, polya_size)] # same as merge, just the data.table way

distribution_full[is.na(n), n := 0] 
distribution_full[, total_reads := sum(n), by = .(Genotype, rep, AGI)]
distribution_full[, perc := 100 * n / total_reads]

distribution_full_mean <- distribution_full[total_reads > 0, .(median_perc = median(perc)), by = .(Genotype, AGI, polya_size, rep)]

distribution_b <- distribution_full_mean[ ,.( mean = mean(median_perc) ), by = .(Genotype, polya_size)]
distribution_b <- unique(distribution_b, by = c("polya_size", "Genotype", "mean"))

rm(all_combinations)

# two-step by hanzhang in dplyr syntax.
# fills the missing values only for 10:200, the rest becomes NA and is dropped. same wight for each AGIn.

distribution_full <- C50 %>% #  mean distri for each AGI
  group_by(Genotype, rep, AGI, polya_size) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(nesting(Genotype, rep, AGI), polya_size = 10:200, fill = list(n = 0)) %>% # this effectively limits the calculation to 10-200 nt
  group_by(Genotype, rep, AGI) %>%
  mutate(
    total_reads = sum(n), 
    perc = ifelse(total_reads == 0, 0, 100 * n / total_reads)
  ) %>%
  ungroup()

distribution_full_mean <- distribution_full %>%  # mean distri for each AGI
  group_by(Genotype, AGI, polya_size) %>%
  summarise(
    mean_perc = mean(perc, na.rm = TRUE),
    .groups = "drop"
  )

distribution_b  <- distribution_full_mean %>% # mean across all mean_AGI_distributions
  group_by(Genotype, polya_size) %>%
  summarise(
    mean = mean(mean_perc, na.rm = TRUE), 
    .groups = "drop"
  )
distribution_b <- as.data.table(distribution_b) # needed because the smooth precalc function wants a data.table


# weighted distri. single step. same weight for each AGI
## same weight for each AGI. instead of summing up the number of reads at a given poly A size, it sums up the weights (are defined as  1/.N) at a given polyA size. There is no need to calculate individual AGI distributions and than the distribution of the means or medians. Faster and less heavy as compared to the code above. the CJ table is smaller because it does not need to be calculated for each AGI.

all_combinations <- CJ(
  Genotype = unique(C50$Genotype),
  rep = unique(C50$rep),
  polya_size = unique(C50$polya_size)
)
weighted_data <- C50[, .(weighted_n = sum(weight)), by = .(Genotype, polya_size, rep)]
weighted_data <- merge(all_combinations, weighted_data, by = c("Genotype", "rep", "polya_size"), all.x = TRUE)
weighted_data <- weighted_data[is.na(weighted_n), weighted_n := 0]
weighted_data[, `:=` (total_weighted = sum(weighted_n), perc = 100 * (weighted_n) / sum(weighted_n)), by = .(Genotype, rep)]

weighted_data[, mean := mean(perc, na.rm = TRUE), by = .(Genotype, polya_size)] # the na.rm= TRUE is needed when the groups are not complete. Example, one Genotype has only one rep. 

# bulk. profile is dominated by abundant AGIs
bulk <- C50[, .(n = .N), by = .(Genotype, polya_size, rep)] 
bulk <- merge(all_combinations, bulk, by = c("Genotype", "rep", "polya_size"), all.x = TRUE) 
bulk <- bulk[is.na(n), n := 0]
bulk[, `:=` (total = sum(n), perc = 100*(n) / sum(n)), by = .(Genotype, rep)] 
bulk[, mean := mean(perc, na.rm = T),  by = .(Genotype, polya_size)]

rm (all_combinations)

#each of these way of calculating can be adapted to include additional groups. For example, if you want to split your AGI by whatsoever category.

# precalc smooth function precalculate smooth for each calculation}
# this can be done inside ggplot, but doing it before is less heavy and faster. works better for me.
precalc_smooth <- function(data, x_col, y_col, group_cols, span = 0.02) {
  smoothed_data <- data[, {
    fit <- loess(get(y_col) ~ get(x_col), data = .SD, span = span)
    smooth <- predict(fit)
    .(x = get(x_col), y = smooth)
  }, by = group_cols]
  return(smoothed_data)
}
grouping_vars <- c("Genotype")

distri_smooth <- precalc_smooth(distri, "polya_size", "mean", grouping_vars)
weighted_data_smooth <- precalc_smooth(weighted_data, "polya_size", "mean", grouping_vars)
bulk_smooth <- precalc_smooth(bulk, "polya_size", "mean", grouping_vars)
b_smooth <- precalc_smooth(distribution_b, "polya_size", "mean", grouping_vars)
rm(distri, weighted_data, bulk, distribution_b); gc()

# plot
plot <- ggplot() +
  
  # intergenic
  geom_line(data = distri_smooth, aes(x = x, y = y, color = Genotype, group = Genotype), alpha = 0.8) +
  # weighted
  geom_line(data = weighted_data_smooth, aes(x = x, y = y, color = Genotype, group = Genotype), linetype = "dotted", alpha = 0.8) +
  # bulk
  geom_line(data = bulk_smooth, aes(x = x, y = y, color = Genotype, group = Genotype), linetype = "dashed", alpha = 0.8) +
  
  # hanzhang
  #   geom_area(data = HY_smooth, aes(x = x, y = y, fill = Genotype, group = Genotype), linetype = "longdash", alpha = 0.3) +
  # b
  geom_area(data = b_smooth, aes(x = x, y = y, fill = Genotype, group = Genotype), linetype = "longdash", alpha = 0.5) +
  
  geom_vline(data = peaks_distri,
             aes(xintercept = peak_x),
             linetype = "dotted",
             color = "black",
             linewidth = 0.5) +
  
  facet_wrap(~ Genotype, ncol = 2) +
  scale_fill_manual(values = genotype_colors, name = NULL) +
  scale_color_manual(values = genotype_colors, guide = "none") +
  scale_x_continuous(
    limits = c(0, 300),
    breaks = custom_breaks,
    labels = every_nth(custom_breaks, 6, inverse = TRUE)
  ) +
  
  theme_classic(base_size = 15) +
  theme(legend.position = "top", legend.justification = "left") +
  
  labs(
    title = " profile calculation methods",
    subtitle = "line = intergenic, dotted = weighted, dashed = bulk, area = HY",
    x = "polyA tail size",
    y = "frequency"
  )

plot

# peak functions
# detect the first peak
find_first_peak <- function(x, y) {
  # Find local maxima
  peaks <- which(diff(sign(diff(y))) == -2) + 1
  
  # If no peaks are found, return NA
  if (length(peaks) == 0) {
    return(list(x = NA_real_, y = NA_real_))
  }
  
  # Select the first peak based on its position in the density (smallest x value)
  first_peak <- peaks[1]
  
  # Return the x and y value of the first peak
  return(list(x = x[first_peak], y = y[first_peak]))
}

peaks_distri <- distri_smooth[, {
  p <- find_first_peak(x, y)
  .(peak_x = p$x, peak_y = p$y)
}, by = Genotype]

peaks_weighted <- weighted_data_smooth[, {
  p <- find_first_peak(x, y)
  .(peak_x = p$x, peak_y = p$y)
}, by = Genotype]

peaks_bulk <- bulk_smooth[, {
  p <- find_first_peak(x, y)
  .(peak_x = p$x, peak_y = p$y)
}, by = Genotype]

peaks_two_step <- b_smooth[, {
  p <- find_first_peak(x, y)
  .(peak_x = p$x, peak_y = p$y)
}, by = Genotype]


peaks_distri[, method := "intergenic"]
peaks_weighted[, method := "weighted"]
peaks_bulk[, method := "bulk"]
peaks_two_step[, method := "two_step"]

peaks_all <- rbindlist(list(peaks_distri, peaks_weighted, peaks_bulk, peaks_two_step))


find_single_peak <- function(x, y) {
  # Identify local maxima
  peaks <- which(diff(sign(diff(y))) == -2) + 1
  
  if (length(peaks) == 0) {
    return(list(peak_x = NA_real_, peak_y = NA_real_))
  }
  
  # Calculate prominences for each peak
  prominences <- y[peaks] - sapply(peaks, function(peak) {
    min(c(min(y[1:peak]), min(y[(peak + 1):length(y)])))
  })
  
  # Pick the most prominent peak
  best <- peaks[which.max(prominences)]
  
  list(peak_x = x[best], peak_y = y[best])
}

peaks_distri <- distri_smooth[, {
  p <- find_single_peak(x, y)
  .(peak_x = p$peak_x, peak_y = p$peak_y)
}, by = Genotype]


# detect the most prom peak between 35 and 60
find_second_peak <- function(x, y) {
  # Restrict to the desired range
  in_range <- which(x >= 35 & x <= 60)
  if (length(in_range) < 3) {
    return(list(peak_x = NA_real_))
  }
  
  x_sub <- x[in_range]
  y_sub <- y[in_range]
  
  # Identify peaks: local maxima with a slope change
  peaks <- which(diff(sign(diff(y_sub))) == -2) + 1
  if (length(peaks) == 0) {
    return(list(peak_2 = NA_real_))
  }
  
  # Calculate prominences for each peak
  prominences <- sapply(peaks, function(peak) {
    left_min <- min(y_sub[1:peak])
    right_min <- min(y_sub[(peak + 1):length(y_sub)])
    y_sub[peak] - min(c(left_min, right_min))
  })
  
  # Find the most prominent peak
  most_prominent_index <- which.max(prominences)
  peak_2_value <- x_sub[peaks[most_prominent_index]]
  
  return(list(peak_2 =peak_2_value))
}

peaks_distri_2 <- bulk_smooth[, {
  p <- find_second_peak(x, y)
  .(peak_2 = p$peak_2)
}, by = Genotype]


## percent short reads by AGI
Table_short <- C50[, .(total_reads = .N, short_reads = sum(polya_size < 30, na.rm = TRUE), perc_short = 100 * sum(polya_size < 30, na.rm = TRUE) / .N
), by = .(AGI, rep, Genotype)]
Table_short[, mean_short := mean(perc_short, na.rm = TRUE), by = .(AGI, Genotype)]
mean_Table_short <- unique(Table_short[, .(AGI, mean_short, Genotype)], by = c("AGI", "Genotype"))
mean_Table_short$mean_short <- round(mean_Table_short$mean_short, 2)
Table_short <- Table_short[, .(AGI, Genotype, rep, perc_short)]

# short count
short_count <- C50[, .(total_reads = .N, short_reads = sum(polya_size < 30 , na.rm = TRUE), perc_short = 100 * sum(polya_size < 30, na.rm = TRUE) / .N
), by = .(rep, Genotype)]


# barplot short
plot <- ggplot()+
  geom_bar (data = short_count, aes(x = rep, y = perc_short, fill= Genotype), stat = "identity", position = position_dodge(width = 0.9)) +
  
  facet_wrap(~Genotype, nrow = 1)+
  
  scale_fill_manual(values = genotype_colors, name = NULL)+
  scale_y_continuous(limits = c(0,30), expand = c(0, 0))+
  guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1)))+
  # theme_no_grid()+
  theme_classic(base_size = 18)+
  
  theme (legend.position="none")+
  theme(axis.text.x = element_blank())+
  labs(
    
    x = "",
    y= "% short tails")

plot
#ggsave(filename ="pdf_R10_rosette_bargraph_short.pdf", plot , width = 6, height = 6, dpi = 1200)

# box short

box <- ggplot(Table_short, aes(x = rep, y = perc_short, color = Genotype)) +
  
  geom_boxplot(position = position_dodge(0.9), outlier.shape = NA,  fill = "white",  width=0.5, linewidth = 1, alpha=0.4, notch=F)+
  facet_wrap(~ Genotype, nrow = 1)+
  scale_y_continuous(limits = c(0, 50))+
  scale_color_manual(values = genotype_colors, name = NULL)+
  scale_fill_manual(values = genotype_colors, name = NULL)+
  # theme_no_grid()+
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())+
  theme (legend.position="none")+
  
  labs(
    #x="",
    y= "% short")
box
#ggsave(filename ="pdf_R10_rosette_boxshort.pdf", box, width = 6, height = 6, dpi = 1200)


## percent medium reads 
Table_medium <- C50[, .(total_reads = .N, medium_reads = sum(polya_size >= 35 & polya_size <= 60 , na.rm = TRUE), perc_medium = 100 * sum(polya_size >= 35 & polya_size <= 60, na.rm = TRUE) / .N
), by = .(AGI, rep, Genotype)]
Table_medium[, mean_medium := mean(perc_medium, na.rm = TRUE), by = .(AGI, Genotype)]
mean_Table_medium <- unique(Table_medium[, .(AGI, mean_medium, Genotype)], by = c("AGI", "Genotype"))
mean_Table_medium$mean_medium <- round(mean_Table_medium$mean_medium, 2)
Table_medium <- Table_medium[, .(AGI, Genotype, rep, perc_medium)]

# medium count
medium_count <- C50[, .(total_reads = .N, medium_reads = sum(polya_size >= 35 & polya_size <= 60 , na.rm = TRUE), perc_medium = 100 * sum(polya_size >= 35 & polya_size <= 60, na.rm = TRUE) / .N
), by = .(rep, Genotype)]


# barlot medium plot
plot <- ggplot()+
  geom_bar (data= medium_count, aes(x = rep, y = perc_medium, fill= Genotype), stat = "identity", position = position_dodge(width = 0.9)) +
  
  facet_wrap(~Genotype, nrow = 1)+
  
  scale_fill_manual(values = genotype_colors, name = NULL)+
  scale_y_continuous(limits = c(0,32), expand = c(0, 0))+
  guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1)))+
  # theme_no_grid()+
  theme_classic(base_size = 18)+
  
  theme (legend.position="none")+
  theme(axis.text.x = element_blank())+
  
  theme (legend.position="none")+
  
  labs(
    
    x = "",
    y= "% medium reads")

plot
#ggsave(filename ="pdf_R10_rosette_bar_medium.pdf", plot , width = 6, height = 6, dpi = 1200)


## percent long reads
Table_long <- C50[, .(total_reads = .N, long_reads = sum(polya_size > 60 , na.rm = TRUE), perc_long = 100 * sum(polya_size > 60, na.rm = TRUE) / .N
), by = .(AGI, rep, Genotype)]
Table_long[, mean_long := mean(perc_long, na.rm = TRUE), by = .(AGI, Genotype)]
mean_Table_long <- unique(Table_long[, .(AGI, mean_long, Genotype)], by = c("AGI", "Genotype"))
mean_Table_long$mean_long <- round(mean_Table_long$mean_long, 2)
Table_long <- Table_long[, .(AGI, Genotype, rep, perc_long)]

# long count
long_count <- C50[, .(total_reads = .N, long_reads = sum(polya_size > 60 , na.rm = TRUE), perc_long = 100 * sum(polya_size > 60, na.rm = TRUE) / .N
), by = .(rep, Genotype)]

# barlot long plot}
plot <- ggplot()+
  geom_bar (data = long_count, aes(x = rep, y = perc_long, fill= Genotype), stat = "identity", position = position_dodge(width = 0.9)) +
  
  facet_wrap(~Genotype, nrow = 1)+
  
  scale_fill_manual(values = genotype_colors, name = NULL)+
  scale_y_continuous(limits = c(0,80), expand = c(0, 0))+
  guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1)))+
  # theme_no_grid()+
  theme_classic(base_size = 18)+
  
  theme (legend.position="none")+
  theme(axis.text.x = element_blank())+
  
  theme (legend.position="none")+
  
  labs(
    
    x = "",
    y= "% long reads")

plot
#ggsave(filename ="pdf_R10_degradation_dcp5_bargraph_long.pdf", plot , width = 6, height = 6, dpi = 1200)


## RPM
# Table_reads
Table_reads<- C50[, .(AGI_reads = .N),  by = .(AGI, rep, Genotype)]
total <- C50[, .(total= .N), by = .(rep, Genotype)]
Table_reads <- merge(Table_reads, total, by = c("rep", "Genotype"))
Table_reads[, RPM := 1000000 *AGI_reads/total]
Table_reads[, mean_RPM := mean(RPM, na.rm = TRUE), by = .(AGI, Genotype)]
mean_Table_reads <- unique(Table_reads[, .(AGI, mean_RPM, Genotype)], by = c("AGI", "Genotype"))
mean_Table_reads$mean_RPM <- round(mean_Table_reads$mean_RPM, 2)
Table_reads <- Table_reads[, .(AGI, Genotype, rep, AGI_reads, RPM)]


## peaks per replicate
# density by Agi and rep
DT_density <- C50[, .(density_C50 = list(density(polya_size, bw = "sj", adjust = 1))), by = .(AGI, Genotype, rep)]

# Extract x and y from density_C50 list column
density_table <- DT_density[, {
  x <- density_C50[[1]]$x
  y <- density_C50[[1]]$y
  .(x, y)
}, by = .(AGI, Genotype, rep)]

# calculate density and retrieve the results  to a table

all_1_pooled_peaks<- density_table[, {maxima <- find_first_peak(x, y)}, by = .(AGI, Genotype, rep)]
all_1_pooled_peaks[, peak_1:= round(x, 0)]
all_1_pooled_peaks[,x := NULL]

all_p_pooled_peaks <- density_table[, {maxima <- find_single_peak(x, y)}, by = .(AGI, Genotype, rep)]
all_p_pooled_peaks[, peak_p := round(peak_x, 0)]
all_p_pooled_peaks[, peak_x := NULL]

Table_peaks <- merge(all_1_pooled_peaks, all_p_pooled_peaks, by = c("AGI", "Genotype","rep"))
Table_peaks <- Table_peaks[, .(AGI, Genotype, rep, peak_1, peak_p)]
rm(all_1_pooled_peaks, all_p_pooled_peaks, DT_density, density_table)

## profile classification
#profile_1
Table_peaks[, profile_1 := fcase(
  peak_p < 35, "short",
  peak_p >= 35 & peak_p <= 60, "inter",
  peak_p > 60, "long"
)]
Table_peaks[, profile_1 := factor(profile_1, levels  = c("short", "inter", "long"))]

table_AGI <- merge(Table_short, Table_medium, by = c("AGI", "Genotype", "rep"))
table_AGI <- merge(table_AGI, Table_long, by = c("AGI", "Genotype", "rep"))
table_AGI <- merge(table_AGI, Table_reads, by = c("AGI", "Genotype", "rep"))
#table_AGI <- merge(table_AGI, Table_median, by = c("AGI", "Genotype", "rep"))
table_AGI <- merge(table_AGI, Table_peaks, by = c("AGI", "Genotype", "rep"))
table_AGI <-table_AGI[, .(AGI, 
                          Genotype, 
                          rep, 
                         # median, 
                          perc_short, 
                          perc_medium, 
                          perc_long, 
                          perc_U, 
                          RPM, 
                          peak_1, 
                          peak_p,  
                          profile_1)]