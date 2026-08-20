library(data.table)
library(arrow)     # For Out-of-Memory processing
library(dplyr)     # For Arrow "remote control"
library(pbapply)   # For progress bars
library(ggplot2)
library(multcomp)
library(UpSetR)
library(stringr)
library(patchwork)
#library(ComplexUpset)


# Import files U tails percentages
#Table_U_HZ <- fread("Table_U_HZ.csv")
Utails_per_gene <- 'Uridylation_%'
file_Utails_per_gene <- list.files(Utails_per_gene, full.names = TRUE)
list_Utails_per_gene <- lapply(file_Utails_per_gene, function (x) fread(x))
Table_U_HZ <- do.call("rbind",list_Utails_per_gene)
length(unique(Table_U_HZ$Feature.name))

#colnames(Table_U_HZ)
#unique(Table_U_HZ$Tail_tag)
#unique(Table_U_HZ$Percent)
#sort(unique(Table_U_HZ$Percent))
 
Table_U_HZ <- Table_U_HZ[Genotype %in% c("Col0", "xrn4", "xrn4_dne1") & Condition == "Seedling"]

colnames(Table_U_HZ)

summarize_table <- function(df, name) {
  data.table::data.table(
    source = name,
    n_rows = nrow(df),
    n_transcripts = length(unique(df$Feature.name)),
    n_zero = sum(df$Percent == 0, na.rm = TRUE),
    mean_percent = mean(df$Percent, na.rm = TRUE),
    median_percent = median(df$Percent, na.rm = TRUE)
  )
}

cmp <- rbind(
  summarize_table(Table_U_HZ, "r_scripts/Table_f_HZ_diff_Uridylation_analysis.R"),
  summarize_table(Table_U, "full_pipeline_uridylation_analysis.R")
)


A <- unique(Table_U_HZ[, .(Feature.name, Genotype, Rep, Tail_tag, Percent)])
B <- unique(Table_U[, .(Feature.name, Genotype, Rep, Tail_tag, Percent)])

setkey(A, Feature.name, Genotype, Rep, Tail_tag)
setkey(B, Feature.name, Genotype, Rep, Tail_tag)

cmp_rows <- merge(
  A[, .(Feature.name, Genotype, Rep, Tail_tag, Percent_A = Percent)],
  B[, .(Feature.name, Genotype, Rep, Tail_tag, Percent_B = Percent)],
  by = c("Feature.name", "Genotype", "Rep", "Tail_tag"),
  all = TRUE
)


cmp_summary <- cmp_rows[, .(
  n_in_A = sum(!is.na(Percent_A)),
  n_in_B = sum(!is.na(Percent_B)),
  overlap = sum(!is.na(Percent_A) & !is.na(Percent_B)),
  only_A = sum(!is.na(Percent_A) & is.na(Percent_B)),
  only_B = sum(is.na(Percent_A) & !is.na(Percent_B))
)]

#write.table(Table_U_HZ, "Table_U_HZ.tsv", sep='\t', row.names=FALSE, quote=FALSE)

my_pRepare <- function(x) { 
  
  x$variable <- x$Genotype
  x$variable <- factor(x$variable, levels=c("Col0", "xrn4", "xrn4_dne1"))
  
  x$Percent_co <- ifelse(x$Percent==0, 0.1/x$total*100, x$Percent)
  
  x$Uridylated <- as.numeric(x$Percent_co)
  x$NonUridylated <- as.numeric(100 - x$Percent_co)
  return(x)
}

Table <- my_pRepare(Table_U_HZ)

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

table_glht <- function(x) {
  pq <- summary(x)$test
  mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
  error <- attr(pq$pvalues, "error")
  colnames(mtests) <- c("Estimate", "Std.Error", "z.value", "p.value")
  return(mtests)
}

my_glm_test <- function(x) {
  
  Result_table_exp <- data.frame(Estimate=as.numeric(), Std.Error=as.numeric(), z.value=as.numeric(), p.value=as.numeric(), comparaison=as.character(), Feature.name=as.character(), padj=as.character(), Condition=as.character())
  
  for (t in unique(x$Condition))  {
    print(t)
    
    Result_table <- data.frame(Estimate=as.numeric(), Std.Error=as.numeric(), z.value=as.numeric(), p.value=as.numeric(), comparaison=as.character(), Feature.name=as.character())
    Selec_condition <- x[x$Condition==t,]
    Selec_condition <- filter_A_for_all_BC(Selec_condition)
    
    
    for (i in unique(Selec_condition$Feature.name)) {
      print(i)
      PermRNA_percent_50_tailed_U_gene <- Selec_condition[Selec_condition$Feature.name==i, ]
      if (nrow(PermRNA_percent_50_tailed_U_gene) <= 2) {next}
      tryCatch({     
      Uri1 <- glm(cbind(Uridylated , NonUridylated ) ~ variable, family=quasibinomial, data=PermRNA_percent_50_tailed_U_gene)
      glht_obj <- summary(glht(Uri1 , mcp(variable="Tukey")))
      Tab <- as.data.frame(table_glht(glht_obj))
      Tab$comparaison <- row.names(Tab)
      Tab$Feature.name <- i
      rownames(Tab) <- NULL
      Result_table <- rbind(Result_table, Tab)
    }, error = function(e) {
      message("Skipping ", i, ": ", e$message)
    })
    }
    
    Result_table$padj <- p.adjust(Result_table$p.value, "BH")
    Result_table$Condition <- t
    Result_table_exp <- rbind(Result_table_exp, Result_table)
    
  }
  
  return(Result_table_exp)
}

Stat <- my_glm_test(Table)

Stat$diff <- "no-sign"
Stat[Stat$padj<=0.05,]$diff <- "sign"

Stat$Change <- "no-change"
Stat[Stat$Estimate<0,]$Change <- "decrease"
Stat[Stat$Estimate>0,]$Change <- "increase"

Stat$diff <- as.factor(Stat$diff)
Stat$Change <- as.factor(Stat$Change)

Count_ofDUgenes <- Stat %>% group_by(comparaison, Condition, diff, Change, .drop=FALSE) %>%
  summarise(number=n())  %>%
  group_by(comparaison, Condition) %>%
  mutate(total=sum(number), Percent=100*number/sum(number))

Diff_uri_transcripts <- Stat[Stat[["diff"]] == "sign" & Stat[["Change"]] == "increase", ]

Diff_Table_U_HZ <- Table_U_HZ[Feature.name %in% unique(Diff_uri_transcripts$Feature.name)]


Count_ofDUgenes$comparaison <- factor(Count_ofDUgenes$comparaison, levels=c("xrn4 - Col0", "xrn4_dne1 - Col0", "xrn4_dne1 - xrn4"))
Count_ofDUgenes_diff <- Count_ofDUgenes[Count_ofDUgenes$diff=="sign",]


write.table(Count_ofDUgenes_diff, "r_tables/Countofdifferentially_uridylatedmRNA_HZ.txt", sep='\t', row.names=FALSE, quote=FALSE)

my_theme <- theme(panel.spacing = unit(0.1, "lines"),
                  plot.title = element_text(hjust = 0.5),
                  panel.background = element_blank(),
                  panel.grid.major.y = element_line(linewidth = 0.2, linetype = 'dashed', colour = "gray70"),
                  panel.border = element_rect(colour="black",fill=NA,linewidth=0.5),
                  axis.text.y = element_text(color="black",size=8),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title = element_blank(),
                  #legend.position="none",
                  legend.text = element_text(color="black",size=7),
                  legend.title = element_text(color="black",size=7),
                  legend.key.size = unit(0.8,"line"),
                  strip.text = element_blank())


my_colors <- c("#E69F00", "#56B4E9", "#CC79A7", "#009E73")


Count_ofDUgenes_diff_selec_UP <- Count_ofDUgenes_diff[Count_ofDUgenes_diff$Change=="increase", ]
Count_ofDUgenes_diff_selec_DOWN <- Count_ofDUgenes_diff[Count_ofDUgenes_diff$Change=="decrease", ]

p1 <- ggplot(Count_ofDUgenes_diff_selec_DOWN, aes(x=comparaison, y=Percent, fill=comparaison)) + 
  geom_bar(stat="identity") +
  my_theme +
  theme(axis.title.y = element_text(size = 12)) +
  geom_text(aes(label=number), position=position_dodge(width=0.9), vjust=+2, hjust=0.5, size=2) +
  scale_fill_manual(values = my_colors) +
  scale_y_reverse(limits=c(50,0), breaks = seq(50, 0, by = -10)) +
  
  labs(
    title = "Less uridylatd transcripts",
    x = "Comparison",
   # y = "% of genes",
    fill = "Comparison"
  )


p2 <- ggplot(Count_ofDUgenes_diff_selec_UP, aes(x=comparaison, y=Percent, fill=comparaison)) + 
  geom_bar(stat="identity") +
  my_theme +
  theme(axis.title.y = element_text(size = 15)) +
  geom_text(aes(label=number), position=position_dodge(width=0.9), vjust=-0.5, hjust=0.5, size=2) +
  scale_y_continuous(limits=c(0, 50), breaks=seq(0,50,by=10)) +
  scale_fill_manual(values = my_colors) +
 
  
  labs(
    title = "More uridylated transcripts",
    x = "Comparison",
 #   y = "% of genes",
    fill = "Comparison"
  )



library(gridExtra)
p <- grid.arrange(p2, p1, ncol=1, nrow = 2, widths=c(8), heights=c(5,5))

#ggsave("r_plots/Percentofdifferentiallyuridylated.pdf", p, width = 5, height = 5, dpi = 300)
#svg
#ggsave("r_plots/Percentofdifferentiallyuridylated.svg", p, width = 5, height = 5, dpi = 300)


# box U
box <- ggplot(Table_U_HZ, aes(x = Rep, y = Percent, color = Genotype, fill = Genotype)) +
  
  geom_boxplot(position = position_dodge(0.9), outlier.shape = NA,
               width = 0.5, linewidth = 1, alpha = 0.4, notch = F) +
  scale_y_continuous(limits = c(0, 50)) +
  scale_color_manual(values = genotype_colors, name = "Genotype") +
  scale_fill_manual(values = genotype_colors, name = "Genotype") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust =1),
    axis.title.x  = element_blank(),
    legend.position  = "right",
    strip.background = element_rect(fill = "white", color = "grey80"),
    strip.text       = element_text(face = "bold")                      
  ) +
  labs(
    title = "Percentage of uridylated poly(A) tails",  
    x = "",
    y = "% uridylation"
  )
box
#pdf
#ggsave("r_plots/boxplot_U_per_gene.pdf", box, width = 7, height = 6, dpi = 1200)
#svg
#ggsave("r_plots/boxplot_U_per_gene.svg", box, width = 7, height = 6, dpi = 1200)





# Uridylation in diff transcripts 
box <- ggplot(Diff_Table_U_HZ, aes(x = Genotype, y = Percent, color = Genotype, fill = Genotype)) +
  
  geom_boxplot(position = position_dodge(0.9), outlier.shape = NA,
               width = 0.5, linewidth = 1, alpha = 0.4, notch = F) +
  scale_y_continuous(limits = c(0, 50)) +
  scale_color_manual(values = genotype_colors, name = "Genotype") +
  scale_fill_manual(values = genotype_colors, name = "Genotype") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust =1),
    axis.title.x  = element_blank(),
    legend.position  = "right",
    strip.background = element_rect(fill = "white", color = "grey80"),
    strip.text       = element_text(face = "bold")                      
  ) +
  labs(
    title = "uri",  
    x = "",
    y = "% uridylation"
  )
box

unique(Diff_Table_U_HZ$Feature.name)