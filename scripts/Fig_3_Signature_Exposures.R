signatures <- read.csv('~/Local_Documents/CN_Signatures/final_sigs/signature_exposures_component_models_britroc_aCNs_component_by_signature_britmodelsvansigs5_aCNs_2023-08-01.csv')

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
  viridis::scale_fill_viridis(discrete=FALSE) +
  ggplot2::theme(plot.title = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_text(size = 12),
                 axis.text.y = ggplot2::element_text(size = 12),
                 legend.position = "none") +
  ggplot2::labs(fill = '', x = "", y = " ") +
  ggplot2::scale_y_discrete(limits = paste0('S', 1:dim(signatures)[1])) +
  ggplot2::ggtitle("")

g <- g + ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                        axis.ticks.y = ggplot2::element_blank(),
                        axis.text.x = ggplot2::element_text(size = 12),
                        axis.title.x = ggplot2::element_text(size = 12),
                        axis.text.y = ggplot2::element_blank(),
                        legend.position = "none") +
  ggplot2::coord_flip()

sample <- read.csv("~/Local_Documents/CN_Signatures/data/sample.csv", sep="")

cancer_type <- read.csv("~/Local_Documents/CN_Signatures/outputs/cancer_type_rescored_histotypes_187.csv")
gene_tracks <- read.csv("~/Local_Documents/CN_Signatures/outputs/gene_CN_status/X/gene_tracks_187.csv")
gene_tracks_top_6 <- read.csv("~/Local_Documents/CN_Signatures/outputs/gene_CN_status/X/gene_tracks_top_6_187.csv")
gene_tracks_targeted_panel_seq <- read.csv("~/Local_Documents/CN_Signatures/outputs/targeted_panel_seq_status/gene_tracks_187_targeted_panel_seq_top_5.csv")
survival_tracks <- read.csv("~/Local_Documents/CN_Signatures/outputs/survival_track_Three_Year_PFS_187_rescored.csv")
ploidy <- read.csv("~/Local_Documents/CN_Signatures/outputs/ploidy_classification_187.csv")
p53abn_MUT <- read.csv("~/Local_Documents/CN_Signatures/outputs/p53abn_mut_188.csv")
HER2_Track <- readxl::read_xlsx("~/Downloads/HER2_Track_simplified.xlsx")
tumour_stage <- read.csv("~/Local_Documents/CN_Signatures/outputs/tumour_stage_187.csv")

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
  scale_fill_manual(values = c("black", "chartreuse3", "lightblue", "darkgoldenrod2", "white")) +
  theme(plot.title = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.position = "none") +
  labs(fill = '', x = "", y = "", title = "")  + coord_flip()

r <- ggplot(p53abn_MUT) +
  geom_tile(aes(x = sample, y = 30, fill = p53_mut_type), width = 1, height = 0.2) +
  scale_fill_manual(values = c("red", "chartreuse3", "blue", "darkgoldenrod2", "grey", "white","lightblue", "purple", "darkblue")) +
  theme(plot.title = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.position = "none") +
  labs(fill = '', x = "", y = "", title = "")  + coord_flip()

t <- ggplot(ploidy) +
  geom_tile(aes(x = sample, y = 30, fill = classification), width = 1, height = 0.2) +
  scale_fill_manual(
    values = c(
      "1" = "red", "2" = "grey", "3" = "lightblue", "4" = "blue", "5+" = "black"
    )
  ) +
  theme(plot.title = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.position = "none") +
  labs(fill = '', x = "", y = "", title = "") + coord_flip()

w <- ggplot(cancer_type) +
  geom_tile(aes(x = sample, y = 25, fill = str_wrap(cancer_type, 20)), width = 1, height = 0.2) +
  scale_fill_manual(values = c("black", "chartreuse3", "darkgoldenrod2", "red", "lightblue")) +
  theme(plot.title = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size = 12),
      axis.text.y = element_blank(),
      legend.position = "none") +
  labs(fill = '', x = "", y = "", title = "")  + coord_flip()

y <- ggplot(tumour_stage) +
  geom_tile(aes(x = sample, y = 30, fill = stage_full), width = 1, height = 0.2) +
  scale_fill_manual(values = c("black", "green", "purple", "lightblue", "white")) +
  theme(plot.title = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.position = "none") +
  labs(fill = '', x = "", y = "", title = "")  + coord_flip()

gene_plot <- ggplot2::ggplot(gene_tracks, ggplot2::aes(sample, gene, fill=classification)) +
  ggplot2::geom_tile() +
  scale_fill_manual(
    values = c(
      "Gain" = "lightblue", "Amplification" = "blue", "No Data" = "white", "Normal" = "grey", "Loss" = "red", "High Amplification" = "black"
    ), breaks = c('High Amplification', 'Amplification', 'Gain', 'Normal', 'Loss', 'No Data')
  ) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_text(size = 10),
                 axis.text.y = ggplot2::element_text(size = 8),
                 legend.position = "none") +
  ggplot2::labs(fill = 'Status', x = "", y = " ") +
  ggplot2::ggtitle("")

gene_plot <- gene_plot + ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                                        axis.ticks.y = ggplot2::element_blank(),
                                        axis.text.x = ggplot2::element_text(size = 8),
                                        axis.title.x = ggplot2::element_text(size = 10),
                                        axis.text.y = ggplot2::element_blank(),
                                        legend.position = "none") +
  ggplot2::coord_flip()

gene_plot <- gene_plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

gene_plot_top_6 <- ggplot2::ggplot(gene_tracks_top_6, ggplot2::aes(sample, gene, fill=classification)) +
  ggplot2::geom_tile() +
  scale_fill_manual(
    values = c(
      "Gain" = "lightblue", "Amplification" = "blue", "No Data" = "white", "Normal" = "grey", "Loss" = "red", "High Amplification" = "black"
    ), breaks = c('High Amplification', 'Amplification', 'Gain', 'Normal', 'Loss', 'No Data')
  ) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_text(size = 10),
                 axis.text.y = ggplot2::element_text(size = 8),
                 legend.position = "none") +
  ggplot2::labs(fill = 'Status', x = "", y = " ") +
  ggplot2::ggtitle("")

gene_plot_top_6 <- gene_plot_top_6 + ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                                                    axis.ticks.y = ggplot2::element_blank(),
                                                    axis.text.x = ggplot2::element_text(size = 8),
                                                    axis.title.x = ggplot2::element_text(size = 10),
                                                    axis.text.y = ggplot2::element_blank(),
                                                    legend.position = "none") +
  ggplot2::coord_flip()

gene_plot_top_6 <- gene_plot_top_6 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

gene_targeted_panel_seq_plot <- ggplot2::ggplot(gene_tracks_targeted_panel_seq, ggplot2::aes(sample, gene, fill=classification)) +
  ggplot2::geom_tile() +
  scale_fill_manual(
    values = c(
      "Absent" = "white", "Present" = "black", "No Data" = 'lightgrey'
    )
  ) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_text(size = 10),
                 axis.text.y = ggplot2::element_text(size = 8),
                 legend.position = 'none') +
  ggplot2::labs(fill = '', x = "", y = " ") +
  ggplot2::ggtitle("")

gene_targeted_panel_seq_plot <- gene_targeted_panel_seq_plot + ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                                                                              axis.ticks.y = ggplot2::element_blank(),
                                                                              axis.text.x = ggplot2::element_text(size = 8),
                                                                              axis.title.x = ggplot2::element_text(size = 10),
                                                                              axis.text.y = ggplot2::element_blank(),
                                                                              legend.position = 'none') +
  ggplot2::coord_flip()

gene_targeted_panel_seq_plot <- gene_targeted_panel_seq_plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

survival_track <- ggplot2::ggplot(survival_tracks, ggplot2::aes(sample, year, fill=classification)) +
  ggplot2::geom_tile() +
  scale_fill_manual(
    values = c(
      "No Event" = "aquamarine3", "No Data" = "white", "Any Event" = "brown1"
    )
  ) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_text(size = 10),
                 axis.text.y = ggplot2::element_text(size = 6),
                 legend.position = "none") +
  ggplot2::labs(fill = 'Survival Status', x = "", y = " ") +
  ggplot2::ggtitle("")

survival_track <- survival_track + ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                                                  axis.ticks.y = ggplot2::element_blank(),
                                                  axis.text.x = ggplot2::element_text(size = 8),
                                                  axis.title.x = ggplot2::element_text(size = 10),
                                                  axis.text.y = ggplot2::element_blank(),
                                                  legend.position = "none") +
  ggplot2::coord_flip()

survival_track <- survival_track + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

merged_plot <- r + gene_targeted_panel_seq_plot + gene_plot_top_6 + gene_plot + t + g + survival_track + w + y + q + guide_area() +
  plot_layout(nrow = 1, ncol = 11,
              widths = c(0.5, 3.5, 4, 3, 0.5, 6, 0.5, 0.5, 0.5, 0.5),
              heights = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4),
              guides = 'collect')
merged_plot <- merged_plot + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
merged_plot

ggplot2::ggsave("~/Local_Documents/CN_Signatures/plots/sigs_final_Vancouver.png", plot = merged_plot, width = 16, height = 10)
