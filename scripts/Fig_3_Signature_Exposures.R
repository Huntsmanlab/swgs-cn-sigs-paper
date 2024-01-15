library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggtext)
library(patchwork)

signatures <- read.csv('./sig_exposures/signature_exposures_component_models_britroc_aCNs_component_by_signature_britmodelsvansigs5_aCNs_2023-08-01.csv')

nsigs <- nrow(signatures)
long_data <- tidyr::gather(signatures)
long_data$max_sig <- rep(apply(signatures, 2, function(x) which.max(x)),
                         times = 1,
                         each = nsigs)
long_data$sigs <- rep(1:nsigs,dim(signatures)[2])
colnames(long_data) <- c('X', 'Z', 'max_sig', 'Y')
long_data <- long_data %>% dplyr::arrange(max_sig)
long_data$X <- factor(long_data$X, levels = unique(long_data$X))

g <- ggplot2::ggplot(long_data, ggplot2::aes(X, Y, fill=Z)) +
  ggplot2::geom_tile() +
  viridis::scale_fill_viridis(discrete = FALSE, name = "Signature Exposure") +
  ggplot2::theme(plot.title = element_text(size = 12),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::labs(fill = '', x = "", y = "") +
  ggplot2::scale_y_discrete(limits = paste0('S', 1:dim(signatures)[1])) +
  ggplot2::ggtitle("Signature exposures")

sample <- read.csv("./clinical_data/sample.csv", sep="")

cancer_type <- read.csv("./clinical_data/cancer_type_rescored_histotypes_187.csv")
gene_tracks <- read.csv("./cn_calls/gene_tracks_187.csv")
gene_tracks_top_6 <- read.csv("./cn_calls/gene_tracks_top_6_187.csv")
gene_tracks_targeted_panel_seq <- read.csv("./mutation_calls/gene_Mutation_187_top_5.csv")
survival_tracks <- read.csv("./clinical_data/survival_track_Three_Year_PFS_187_rescored.csv")
ploidy <- read.csv("./cn_calls//ploidy_classification_187.csv")
p53abn_MUT <- read.csv("./mutation_calls/p53abn_mut_188.csv")
HER2_Track <- readxl::read_xlsx("./clinical_data/HER2_Track_simplified.xlsx")
tumour_stage <- read.csv("./clinical_data/tumour_stage_187.csv")

colnames(cancer_type)[1] <- 'sample'
colnames(gene_tracks)[1] <- 'sample'
colnames(gene_tracks_top_6)[1] <- 'sample'
colnames(gene_tracks_targeted_panel_seq)[1] <- 'sample'
colnames(survival_tracks)[1] <- 'sample'
colnames(ploidy)[1] <- 'sample'
colnames(p53abn_MUT)[1] <- 'sample'
colnames(HER2_Track)[1] <- 'sample'
colnames(tumour_stage)[1] <- 'sample'

cancer_type$sample <- gsub("-", "\\.", cancer_type$sample)
gene_tracks$sample <- gsub("-", "\\.", gene_tracks$sample)
gene_tracks_top_6$sample <- gsub("-", "\\.", gene_tracks_top_6$sample)
gene_tracks_targeted_panel_seq$sample <- gsub("-", "\\.", gene_tracks_targeted_panel_seq$sample)
survival_tracks$sample <- gsub("-", "\\.", survival_tracks$sample)
ploidy$sample <- gsub("-", "\\.", ploidy$sample)
p53abn_MUT$sample <- gsub("-", "\\.", p53abn_MUT$sample)
HER2_Track$sample <- gsub("-", "\\.", HER2_Track$sample)
tumour_stage$sample <- gsub("-", "\\.", tumour_stage$sample)

HER2_Track$sample[192] <- "YW.EC002"
HER2_Track <- sample %>% left_join(HER2_Track)
HER2_Track$HER2_IHC_score[which(is.na(HER2_Track$HER2_IHC_score))] <- 'No Data'
HER2_Track$HER2_IHC_score[which(HER2_Track$HER2_IHC_score == 'NA')] <- 'No Data'

p53abn_MUT <- sample %>% left_join(p53abn_MUT)

cancer_type$sample <- factor(cancer_type$sample, levels(g[["data"]][["X"]]))
gene_tracks$sample <- factor(gene_tracks$sample, levels(g[["data"]][["X"]]))
gene_tracks_top_6$sample <- factor(gene_tracks_top_6$sample, levels(g[["data"]][["X"]]))
gene_tracks_targeted_panel_seq$sample <- factor(gene_tracks_targeted_panel_seq$sample, levels(g[["data"]][["X"]]))
ploidy$sample <- factor(ploidy$sample, levels(g[["data"]][["X"]]))
survival_tracks$sample <- factor(survival_tracks$sample, levels(g[["data"]][["X"]]))
p53abn_MUT$sample <- factor(p53abn_MUT$sample, levels(g[["data"]][["X"]]))
HER2_Track$sample <- factor(HER2_Track$sample, levels(g[["data"]][["X"]]))
tumour_stage$sample <- factor(tumour_stage$sample, levels(g[["data"]][["X"]]))

q <- ggplot(HER2_Track) +
  geom_tile(aes(x = sample, y = 30, fill = HER2_IHC_score), width = 1, height = 0.2) +
  scale_fill_manual(values = c("black", "chartreuse3", "lightblue", "darkgoldenrod2", "white"), name = "HER2 IHC") +
  theme(plot.title = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_blank()) +
  labs(fill = '', x = "", y = "", title = "HER2")

r <- ggplot(p53abn_MUT) +
  geom_tile(aes(x = sample, y = 30, fill = p53_mut_type), width = 1, height = 0.2) +
  scale_fill_manual(values = c("red", "chartreuse3", "blue", "darkgoldenrod2", "grey", "white","lightblue", "purple", "darkblue")) +
  theme(plot.title = element_markdown(size = 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_blank()) +
  guides(fill = guide_legend(ncol =  2, title = "*TP53* Mutation", title.theme = element_markdown())) +
  labs(fill = '', x = "", y = "", title = "*TP53*")

t <- ggplot(ploidy) +
  geom_tile(aes(x = sample, y = 30, fill = classification), width = 1, height = 0.2) +
  scale_fill_manual(
    values = c(
      "1" = "red", "2" = "grey", "3" = "lightblue", "4" = "blue", "5+" = "black"
    ),
    name = "Ploidy") +
  theme(plot.title = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_blank()) +
  labs(fill = '', x = "", y = "", title = "Ploidy")

w <- ggplot(cancer_type) +
  geom_tile(aes(x = sample, y = 25, fill = str_wrap(cancer_type, 20)), width = 1, height = 0.2) +
  scale_fill_manual(values = c("black", "chartreuse3", "darkgoldenrod2", "red", "lightblue"), name = "Histotype") +
  theme(plot.title = element_text(size = 12),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size = 12),
      axis.text.y = element_blank()) +
  labs(fill = '', x = "", y = "", title = "Histotype")

y <- ggplot(tumour_stage) +
  geom_tile(aes(x = sample, y = 30, fill = stage_full), width = 1, height = 0.2) +
  scale_fill_manual(values = c("black", "green", "purple", "lightblue", "white"), name = "FIGO Stage") +
  theme(plot.title = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_blank()) +
  labs(fill = '', x = "", y = "", title = "Stage")

gene_plot <- ggplot2::ggplot(gene_tracks, ggplot2::aes(sample, gene, fill=classification)) +
  ggplot2::geom_tile() +
  scale_fill_manual(
    values = c(
      "Gain" = "lightblue", "Amplification" = "blue", "No Data" = "white", "Normal" = "grey", "Loss" = "red", "High Amplification" = "black"
    ), breaks = c('High Amplification', 'Amplification', 'Gain', 'Normal', 'Loss', 'No Data'),
    name = "Copy Number Status"
  ) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 8),
                 axis.title.y = ggplot2::element_blank(),
                 legend.position = "none") +
  ggplot2::ggtitle("Clinically relevant genes")

gene_plot_top_6 <- ggplot2::ggplot(gene_tracks_top_6, ggplot2::aes(sample, gene, fill=classification)) +
  ggplot2::geom_tile() +
  scale_fill_manual(
    values = c(
      "Gain" = "lightblue", "Amplification" = "blue", "No Data" = "white", "Normal" = "grey", "Loss" = "red", "High Amplification" = "black"
    ),
    breaks = c('High Amplification', 'Amplification', 'Gain', 'Normal', 'Loss', 'No Data'),
    name = "Copy Number Status"
  ) +
  ggplot2::theme(plot.title = ggtext::element_markdown(size = 12),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::labs(fill = 'Status', x = "", y = "") +
  ggplot2::ggtitle("Top 5 mutated genes + *TP53*")

gene_targeted_panel_seq_plot <- ggplot2::ggplot(gene_tracks_targeted_panel_seq, ggplot2::aes(sample, gene, fill=classification)) +
  ggplot2::geom_tile() +
  scale_fill_manual(
    values = c(
      "Absent" = "white", "Present" = "black", "No Data" = 'lightgrey'
    ),
    name = "Mutation Status"
  ) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_text(size = 10),
                 axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::labs(fill = '', x = "", y = "") +
  ggplot2::ggtitle("Top 5 mutated genes")

survival_track <- ggplot2::ggplot(survival_tracks, ggplot2::aes(sample, year, fill=classification)) +
  ggplot2::geom_tile() +
  scale_fill_manual(
    values = c(
      "No Event" = "aquamarine3", "No Data" = "white", "Any Event" = "brown1"
    ),
    name = "Survival (3-year PFS)"
  ) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_text(size = 10),
                 axis.text.y = ggplot2::element_blank()) +
  ggplot2::labs(fill = 'Survival Status', x = "", y = "") +
  ggplot2::ggtitle("Survival (3-year PFS)")

merged_plot <- r + gene_targeted_panel_seq_plot + gene_plot_top_6 + gene_plot + t + g + survival_track + w + y + q + guide_area() +
  plot_layout(nrow = 11, ncol = 1,
              heights = c(0.5, 2.5, 3, 2, 0.5, 6, 0.5, 0.5, 0.5, 0.5, 6),
              widths = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4),
              guides = 'collect')
merged_plot <- merged_plot & theme(legend.position = "bottom", legend.direction = "vertical")
merged_plot

ggplot2::ggsave("~/Downloads/vancouver_sigs_all_tracks.png", plot = merged_plot, width = 15, height = 16)
