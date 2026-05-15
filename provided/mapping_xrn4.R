library (grid)
library(tidyverse)
library(data.table)
library(ggplot2)


Tab <- fread("/Users/hyu/Documents/FlepSeq/xrn4/DataTables/RUN_Remy_seedlings_xrn4_cleaned_frag.txt")

my_theme <- theme(
  panel.spacing = unit(0.4, "lines"),
  plot.title = element_text(hjust = 0.5),
  panel.background = element_blank(),
  panel.grid.major.y = element_line(linewidth = 0.2, linetype = 'dashed', colour = "gray70"),
  panel.border = element_rect(colour="black",fill=NA,linewidth=0.5),
  #legend.position = "none",
  axis.text.x = element_text(angle = 45, hjust = 1),
  #strip.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank())

TabX <- filter(Tab, AGI == "AT1G13245")

setwd("/Users/hyu/Documents/FlepSeq/xrn4/mapping2")


Proportion_tab <- TabX %>% 
  group_by(Genotype, Rep, info_cds) %>% 
  summarise(n=n()) %>%
  group_by(Genotype, Rep) %>% 
  mutate(total=sum(n), perc = 100*n/sum(n))

Proportion_tab_mean <- Proportion_tab %>%
  group_by(Genotype, info_cds) %>%
  summarise_at(vars(perc),list(mean = mean, sd=sd))

Proportion_tab$info_cds <- factor(Proportion_tab$info_cds, levels=c("intact", "Five_in_cds", "Three_in_cds", "Both_in_cds"))
Proportion_tab_mean$info_cds <- factor(Proportion_tab_mean$info_cds, levels=c("intact", "Five_in_cds", "Three_in_cds", "Both_in_cds"))

my_colors <- c("#ABABAB", "#984cff", "#ff4cbb", "#d1d1d1", "#c9a2ff", "#ff9fda")
p <- ggplot(Proportion_tab) + 
  geom_bar(data=Proportion_tab_mean, aes(x=info_cds, y=mean, fill=Genotype), stat="identity", alpha=0.6, position = "dodge") +
  geom_point(data=Proportion_tab, aes(x=info_cds, y=perc, fill=Genotype, color=Genotype), alpha=1, size=0.8, shape=21) +
  geom_text(data=Proportion_tab_mean, aes(x=info_cds, y=mean, label=round(mean, 1)), 
            vjust=-0.5, color="black", size=3) +
  #facet_grid(rows = vars(Genotype)) +
  facet_wrap(~Genotype, nrow=1)+
  #scale_y_continuous(breaks=seq(0,80,by=10)) +
  #coord_cartesian(ylim = c(0, 50)) +
  scale_fill_manual(values = my_colors) +
  my_theme +
  scale_fill_manual(name="Genotype", values = my_colors) +
  scale_color_manual(name="Genotype", values = my_colors) 
p

ggsave(filename ="AT1G13245_info_cds_proportion.pdf", p, width = 10, height = 5, dpi = 300)



# ------------------------------------------------------------------------------ 
# Coverage Plot

library(dplyr)
library(ggplot2)

# df should contain columns: AGI, read_start, read_end, mRNA_start, mRNA_end

my_colors <- c("#ABABAB", "#984cff", "#ff4cbb", "#d1d1d1", "#c9a2ff", "#ff9fda")

df2 <- select(Tab, Genotype, Rep, AGI, mRNA_start, mRNA_end, read_start, read_end, info_cds, info_polyA, Tag)
rm(Tab)

df2$Genotype <- factor(df2$Genotype, levels = c("Col0_root", "xrn4_root", "XRN4dCTRD_root", "Col0_shoot", "xrn4_shoot", "XRN4dCTRD_shoot"))

count <- df3 %>% 
 group_by(AGI, Genotype, Rep) %>% 
  summarise(n=n())%>%
  arrange(desc(n))

most_expressed_100 <- head(unique(count$AGI), 100)


count <- df2 %>% 
  group_by(Genotype) %>% 
  summarise(n=n())


#===============================================================================

df <- filter(df2, AGI == "AT2G43020") 
bin_width <- 5
df <- df %>%
  rowwise() %>%
  mutate(bin_start = floor(read_start / bin_width) * bin_width,
         bin_end = floor(read_end / bin_width) * bin_width)

calculate_coverage <- function(df) {
  coverage_list <- df %>%
    rowwise() %>%
    mutate(positions = list(bin_start:bin_end)) %>%
    ungroup() %>%
    select(Genotype ,positions, info_cds)
  coverage_data <- coverage_list %>%
    group_by(Genotype, info_cds) %>%
    unnest(cols = c(positions)) %>%  # Unnest the positions
    group_by(Genotype, info_cds, positions) %>%
    summarise(depth = n(), .groups = 'drop')  # Count number of reads covering each position
  return(coverage_data)
}

# Calculate coverage for each position
coverage_data <- calculate_coverage(df)
coverage_data$info_cds <- factor(coverage_data$info_cds, levels=c("intact", "Five_in_cds","Three_in_cds","Both_in_cds"))
my_colors <- c("#ABABAB", "#ff1515", "#1538ff", "#ffb115")

p <- ggplot(coverage_data, aes(x = positions, y = depth, group = info_cds, color = info_cds)) +
  geom_line() +  # Line plot for coverage
  #facet_grid(info_polyA~Genotype, scales="free_y") +
  facet_grid(rows = vars(Genotype)) + #, scales = "free_y"
  theme_minimal() +  
  scale_color_manual(values = my_colors) +
  labs(title = "Coverage AT2G43020",
       x = "Position",
       y = "Coverage Depth") +
  theme(plot.title = element_text(hjust = 0.5))
p

ggsave(filename ="AT2G43020_coverage.pdf", p, width = 5, height = 5, dpi = 300)

# DNE1targets: AT1G13245 AT1G78080 AT2G18160 AT2G22430 AT2G47400 AT3G22380 AT4G34138 AT5G19120 AT4G64260. AT2G43020



#------------------


count <- df2 %>% 
  group_by(AGI, Genotype, Rep) %>% 
  summarise(n=n())%>%
  arrange(desc(n))

list <- unique(count$AGI)
most_expressed_100 <- list[1:10]

#most_expressed_100 <- c("AT1G13245", "AT1G78080", "AT2G18160", "AT2G22430", "AT2G47400", "AT3G22380", "AT4G34138", "AT5G19120")

#df <- filter(df2, AGI == "AT5G19120") 
#df <- filter(df, info_cds == "Five_in_cds")
#df <- filter(df, info_cds == "Three_in_cds")

most_expressed_100 <- filter(df2, AGI %in% most_expressed_100) 

bin_width <- 5

#df <- most_expressed_100
most_expressed_100 <- most_expressed_100 %>%
  rowwise() %>%
  mutate(bin_start = floor(read_start / bin_width) * bin_width,
         bin_end = floor(read_end / bin_width) * bin_width)


calculate_coverage <- function(most_expressed_100) {
  coverage_list <- most_expressed_100 %>%
    rowwise() %>%
    mutate(positions = list(bin_start:bin_end)) %>%
    ungroup() %>%
    select(Genotype, AGI, positions, info_cds)
    coverage_data <- coverage_list %>%
    group_by(Genotype, AGI, info_cds) %>%
    unnest(cols = c(positions)) %>%  # Unnest the positions
    group_by(Genotype, AGI, info_cds, positions) %>%
    summarise(depth = n(), .groups = 'drop')  # Count number of reads covering each position
  
  return(coverage_data)
}

# Calculate coverage for each position
coverage_data <- calculate_coverage(most_expressed_100)

coverage_data$info_cds <- factor(coverage_data$info_cds, levels=c("intact", "Five_in_cds","Three_in_cds","Both_in_cds"))


# Plot the coverage using ggplot2
p <- ggplot(coverage_data, aes(x = positions, y = depth, group = info_cds, color = info_cds)) +
  geom_line() +  # Line plot for coverage
  facet_grid(rows = vars(Genotype), cols = vars(AGI), scales = "free") +  # Allow both x and y scales to be free
  theme_minimal() +  # Clean theme
  scale_color_manual(values = my_colors) +
  labs(title = "Coverage (targets of DNE1)",
       x = "Position",
       y = "Coverage Depth") +
  theme(plot.title = element_text(hjust = 0.5))

p

#ggsave(filename ="AT1G13245_(all readsXinfo_polyA)_coverage.pdf", p, width = 8, height = 5, dpi = 300)
ggsave(filename ="targets of DNE1_coverage.pdf", p, width = 24, height = 8, dpi = 300)
