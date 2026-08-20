# RA5 directly Helene script -> gene by gene analysis

library(data.table)
library(arrow)     # For Out-of-Memory processing
library(dplyr)     # For Arrow "remote control"
#library(pbapply)   # For progress bars
library(ggplot2)
library(multcomp)
#library(UpSetR)
library(patchwork)

dir.create("r_plots/objects", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. INITIALIZATION & METADATA
# ==============================================================================
actual_tags   <- c("NO_RA5", "RA5")

SampleTable <- fread("provided/barcode_correspondance.tsv")
SampleTable <- SampleTable[, .(Barcode, Sample_name, Genotype, Condition, Rep)]
Table_RA5 <- fread("r_tables/Table_RA5.txt")

Table_RA5[, variable := factor(Genotype, levels = c("Col0", "xrn4", "xrn4_dne1"))]



cat("Rows in Table:", nrow(Table_RA5), "\n")
cat("Unique genes:", uniqueN(Table_RA5$Feature.name), "\n")


filter_A_for_all_BC <- function(df) {
  # Extraire toutes les combinaisons distinctes de Genotype et Rep
  all_BC <- df %>%
    distinct(Genotype, Rep)
  
  # Pour chaque Feature.name, compter combien de combinaisons (B,C) il couvre
  valid_A <- df %>%
    distinct(Genotype, Rep, Feature.name) %>%
    group_by(Feature.name) %>%
    summarise(n_bc = n()) %>%
    ungroup()
  
  # Nombre total de combinaisons Genotype, Rep dans la table
  total_bc <- nrow(all_BC)
  
  # Garder les Feature.name qui ont toutes les combinaisons Genotype, Rep
  valid_A <- valid_A %>%
    filter(n_bc == total_bc) %>%
    pull(Feature.name)
  
  # Filtrer le df original pour garder uniquement les lignes avec ces A valides
  df %>% filter(Feature.name %in% valid_A)
}

#Function to apply a glm for proportion
table_glht <- function(x) {
  pq <- summary(x)$test
  mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
  error <- attr(pq$pvalues, "error")
  colnames(mtests) <- c("Estimate", "Std.Error", "z.value", "p.value")
  return(mtests)
}


my_glm_test <- function(x) {
  
  Result_table_exp <- data.frame(Estimate=as.numeric(), Std.Error=as.numeric(), 
                                 z.value=as.numeric(), p.value=as.numeric(), 
                                 comparaison=as.character(), Feature.name=as.character(), 
                                 padj=as.character(), Condition=as.character())
  
  for (t in unique(x$Condition))  {

    Result_table <- data.frame(Estimate=as.numeric(), Std.Error=as.numeric(), 
                               z.value=as.numeric(), p.value=as.numeric(), 
                               comparaison=as.character(), Feature.name=as.character())
    
    Selec_condition <- x[x$Condition==t,]
    Selec_condition <- filter_A_for_all_BC(Selec_condition)
    print
    
    for (i in unique(Selec_condition$Feature.name)) {
      print(i)
      PermRNA_percent_50_tailed_U_gene <- Selec_condition[Selec_condition$Feature.name==i, ]
      if (nrow(PermRNA_percent_50_tailed_U_gene) <= 2) {next}
      Uri1 <- glm(RA5_prop ~ variable, 
                  family=quasibinomial, 
                  weights = total_adj,
                  data = PermRNA_percent_50_tailed_U_gene)
      
      glht_obj <- summary(glht(Uri1 , mcp(variable="Tukey")))
      Tab <- as.data.frame(table_glht(glht_obj))
      Tab$comparaison <- row.names(Tab)
      Tab$Feature.name <- i
      rownames(Tab) <- NULL
      Result_table <- rbind(Result_table, Tab)
    }
    
    
    Result_table$padj <- p.adjust(Result_table$p.value, "BH")
    Result_table$Condition <- t
    Result_table_exp <- rbind(Result_table_exp, Result_table)
    
  }
  
  return(Result_table_exp)
}

# Apply hte function
Stat <- my_glm_test(Table_RA5)

cat("Rows in Table:", nrow(Stat), "\n")
cat("Unique genes:", uniqueN(Stat$Feature.name), "\n")
## count number of differentially expressed genes
Stat$diff <- "no-sign"
Stat[Stat$padj<=0.05,]$diff <- "sign"

Stat$Change <- "no-change"
Stat[Stat$Estimate<0,]$Change <- "decrease"
Stat[Stat$Estimate>0,]$Change <- "increase"


Stat$diff <- as.factor(Stat$diff)
Stat$Change <- as.factor(Stat$Change)

Count_ofDRA5genes <- Stat %>% group_by(comparaison, Condition, diff, Change, .drop=FALSE) %>%
  summarise(number=n())  %>%
  group_by(comparaison, Condition) %>%
  mutate(total=sum(number), Percent=100*number/sum(number))

Count_ofDRA5genes$comparaison <- factor(Count_ofDRA5genes$comparaison, levels=c("xrn4 - Col0", 
                                                                                "xrn4_dne1 - Col0", 
                                                                                "xrn4_dne1 - xrn4"))
Count_ofDRA5genes_diff <- Count_ofDRA5genes[Count_ofDRA5genes$diff=="sign",]

# ==============================================================================
# 8. SUMMARY COUNT TABLE
# ==============================================================================
Count_ofDUgenes <- Stat %>%
  group_by(comparaison, Condition, diff, Change, .drop = FALSE) %>%
  summarise(number = n(), .groups = "drop") %>%
  group_by(comparaison, Condition) %>%
  mutate(total   = sum(number),
         Percent = 100 * number / sum(number))

Count_ofDUgenes_diff <- Count_ofDUgenes[Count_ofDUgenes$diff == "sign", ]

nrow(Stat)
nrow(Table)
unique(Table$RA5_tag)
unique(global_counts$RA5_tag)
table(global_counts$RA5_tag)

#hist(Stat$p.value, breaks=50, main="p-value distribution (RA5 only)")
summary(Stat$p.value)
summary(Stat$padj)

# Get the list of significant genes per comparison
sig_genes <- as.data.table(Stat)[diff == "sign", .(Feature.name, comparaison, Change)]

# For each significant gene, get its Percent_co per sample
sig_data <- merge(
  Table[Feature.name %in% unique(sig_genes$Feature.name)],
  sig_genes,
  by = "Feature.name",
  allow.cartesian = TRUE
)
#View(sig_data)
sig_data$Genotype <- factor(sig_data$Genotype, levels = c("Col0", "xrn4", "xrn4_dne1"))
sig_data$Change   <- factor(sig_data$Change, levels = c("increase", "decrease"))

geno_colors <- c("Col0"      = "#888888",
                 "xrn4"      = "#E6A817",
                 "xrn4_dne1" = "#5B9BD5")

RA5 <- ggplot(sig_data, aes(x = Rep, y = RA5_prop * 100, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6) +
  facet_grid(Change ~ Genotype, scales = "free_x") +
  scale_fill_manual(values = geno_colors) +
  coord_cartesian(ylim = c(0, 5)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "RA5% for differentially fragmented transcripts",
       x = NULL, y = "% RA5 reads per gene") +
  theme_bw(base_size = 13) +
  theme(legend.position  = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank())
print(RA5)

ggsave("r_plots/direct_diff_RA5.svg", RA5, width = 5, height = 5)


#----------------------------------------------------------------------------------------


# Graph
my_theme <- theme(panel.spacing = unit(0.1, "lines"),
                  plot.title = element_text(hjust = 0.5),
                  panel.background = element_blank(),
                  panel.grid.major.y = element_line(linewidth = 0.2, linetype = 'dashed', colour = "gray70"),
                  panel.border = element_rect(colour="black",fill=NA,linewidth=0.5),
                  axis.text.y = element_text(color="black",size=8),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  #axis.title = element_blank(),
                  #legend.position="none",
                  legend.text = element_text(color="black",size=7),
                  legend.title = element_text(color="black",size=5),
                  legend.key.size = unit(0.8,"line"),
                  strip.text = element_blank())


my_colors <- c("#E69F00", "#56B4E9", "#CC79A7", "#009E73")


Count_ofDRA5genes_diff_selec_UP <- Count_ofDRA5genes_diff[Count_ofDRA5genes_diff$Change=="increase", ]

Count_ofDRA5genes_diff_selec_DOWN <- Count_ofDRA5genes_diff[Count_ofDRA5genes_diff$Change=="decrease", ]


p1 <- ggplot(Count_ofDRA5genes_diff_selec_DOWN, aes(x=comparaison, y=Percent, fill=comparaison)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=number), position=position_dodge(width=0.9), vjust=+2, hjust=0.5, size=2) +
  scale_fill_manual(values = my_colors) +
  my_theme +
  scale_y_reverse(limits=c(50,0), breaks = seq(50, 0, by = -10)) +
  
  labs(
    title = "Downregulated RA5 tails",
    x = "Comparison",
    #y = "Percent",
    fill = "Comparison"
  )


p2 <- ggplot(Count_ofDRA5genes_diff_selec_UP, aes(x=comparaison, y=Percent, fill=comparaison)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=number), position=position_dodge(width=0.9), vjust=-0.5, hjust=0.5, size=2) +
  scale_y_continuous(limits=c(0, 50), breaks=seq(0,50,by=10)) +
  scale_fill_manual(values = my_colors) +
  my_theme +
  
  labs(
    title = "Upregulated RA5 RA5 tails",
    x = "Comparison",
    # y = "Percent of genes",
    fill = "Comparison"
  )

library(gridExtra)
p <- grid.arrange(p2, p1, ncol=1, nrow = 2, widths=c(8), heights=c(5,5))