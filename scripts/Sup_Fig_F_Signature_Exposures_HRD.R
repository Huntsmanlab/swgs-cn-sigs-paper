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
  ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_text(size = 12),
                 axis.text.y = ggplot2::element_text(size = 12),
                 legend.title = ggplot2::element_text(size = 12, face = 'bold'),
                 legend.key.size = unit(0.4, "cm")) +
  ggplot2::labs(fill = 'Signature \nExposure', x = "", y = " ") +
  ggplot2::scale_y_discrete(limits = paste0('S', 1:dim(signatures)[1])) +
  ggplot2::ggtitle("Signature Exposures")

g <- g + ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                        axis.ticks.y = ggplot2::element_blank(),
                        axis.text.x = ggplot2::element_text(size = 12),
                        axis.title.x = ggplot2::element_text(size = 12),
                        axis.text.y = ggplot2::element_blank(),
                        legend.title = ggplot2::element_text(size = 12),
                        legend.key.size = unit(0.4, "cm")) +
  ggplot2::coord_flip()

sample <- read.csv("~/Local_Documents/CN_Signatures/data/sample.csv", sep="")

gene_tracks <- read.csv("~/Downloads/gene_Mutation_187_HRD.csv")
gene_tracks_targeted_panel_seq <- read.csv("~/Downloads/gene_tracks_targeted_panel_seq_HRD_187.csv")
ploidy <- read.csv("~/Local_Documents/CN_Signatures/outputs/ploidy_classification_187.csv")

colnames(gene_tracks)[1] <- 'sample'
colnames(gene_tracks_targeted_panel_seq)[1] <- 'sample'
colnames(ploidy)[1] <- 'sample'

gene_tracks$sample <- gsub("-", "\\.", gene_tracks$sample)
gene_tracks_targeted_panel_seq$sample <- gsub("-", "\\.", gene_tracks_targeted_panel_seq$sample)
ploidy$sample <- gsub("-", "\\.", ploidy$sample)

gene_tracks$sample <- factor(gene_tracks$sample, levels(g[["data"]][["X"]]))
gene_tracks_targeted_panel_seq$sample <- factor(gene_tracks_targeted_panel_seq$sample, levels(g[["data"]][["X"]]))
ploidy$sample <- factor(ploidy$sample, levels(g[["data"]][["X"]]))

t <- ggplot(ploidy) +
  geom_tile(aes(x = sample, y = 30, fill = classification), width = 1, height = 0.2) +
  scale_fill_manual(
    values = c(
      "1" = "red", "2" = "grey", "3" = "lightblue", "4" = "blue", "5+" = "black"
    )
  ) +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 12, face = 'bold'),
        legend.direction = "vertical",
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(0.01, 0.01, 0.01, 0.01, "cm")) +
  labs(fill = 'Ploidy', x = "", y = "", title = "Ploidy") + coord_flip()

gene_plot <- ggplot2::ggplot(gene_tracks, ggplot2::aes(sample, gene, fill=classification)) +
  ggplot2::geom_tile() +
  scale_fill_manual(
    values = c(
      "Gain" = "lightblue", "Amplification" = "blue", "No Data" = "white", "Normal" = "grey", "Loss" = "red", "High Amplification" = "black"
    ), breaks = c('High Amplification', 'Amplification', 'Gain', 'Normal', 'Loss', 'No Data')
  ) +
  scale_y_discrete(labels = expression(italic(ATM), italic(ATRX), italic(BAP1), italic(BLM), italic(BRCA1),
                                       italic(BRCA2), italic(BRIP1), italic(CHEK1), italic(CHEK2),
                                       italic(EXO1), italic(FANCC), italic(FANCD2), italic(FANCE),
                                       italic(FANCF), italic(FANCG), italic(MRE11A), italic(PALB2),
                                       italic(RAD50), italic(RAD51))) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_text(size = 12),
                 axis.text.y = ggplot2::element_text(size = 9),
                 legend.title = ggplot2::element_text(size = 12, face = 'bold'),
                 legend.key.size = unit(0.4, "cm")) +
  ggplot2::labs(fill = 'Copy Number Status', x = "", y = " ") +
  ggplot2::ggtitle("Copy Number Status: \nHRD genes")

gene_plot <- gene_plot + ggplot2::theme(plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
                                        axis.ticks.y = ggplot2::element_blank(),
                                        axis.text.x = ggplot2::element_text(size = 9),
                                        axis.title.x = ggplot2::element_text(size = 12),
                                        axis.text.y = ggplot2::element_blank(),
                                        legend.title = ggplot2::element_text(size = 12),
                                        legend.key.size = unit(0.4, "cm")) +
  ggplot2::coord_flip()

gene_plot <- gene_plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

gene_targeted_panel_seq_plot <- ggplot2::ggplot(gene_tracks_targeted_panel_seq, ggplot2::aes(sample, gene, fill=classification)) +
  ggplot2::geom_tile() +
  scale_fill_manual(
    values = c(
      "Absent" = "white", "Present" = "black", "No Data" = 'lightgrey'
    )
  ) +
  scale_y_discrete(labels = expression(italic(ATM), italic(ATRX), italic(BAP1), italic(BLM), italic(BRCA1),
                                       italic(BRCA2), italic(BRIP1),
                                       italic(EXO1), italic(FANCA), italic(FANCC), italic(FANCD2), italic(FANCE),
                                       italic(FANCG), italic(PALB2),
                                       italic(RAD50), italic(RAD51))) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_text(size = 12),
                 axis.text.y = ggplot2::element_text(size = 9),
                 legend.title = ggplot2::element_text(size = 12, face = 'bold'),
                 legend.key.size = unit(0.4, "cm"),
                 legend.direction = 'vertical') +
  ggplot2::labs(fill = 'Mutation Status', x = "", y = " ") +
  ggplot2::ggtitle("Targeted Panel mutation: \nHRD genes") + guides(color = guide_legend(override.aes = list(linetype = c(1, 0, 0))))

gene_targeted_panel_seq_plot <- gene_targeted_panel_seq_plot + ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                                                                              axis.ticks.y = ggplot2::element_blank(),
                                                                              axis.text.x = ggplot2::element_text(size = 9),
                                                                              axis.title.x = ggplot2::element_text(size = 12),
                                                                              axis.text.y = ggplot2::element_blank(),
                                                                              legend.title = ggplot2::element_text(size = 12),
                                                                              legend.key.size = unit(0.4, "cm"),
                                                                              legend.direction = 'vertical') +
  guides(color = guide_legend(override.aes = list(shape = c(15, 15, 15), size = 2))) +
  ggplot2::coord_flip()

gene_targeted_panel_seq_plot <- gene_targeted_panel_seq_plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

merged_plot <- g + t + gene_plot + gene_targeted_panel_seq_plot +
  plot_layout(nrow = 1, ncol = 4,
              widths = c(6, 0.5, 4, 3.5),
              heights = c(4, 4, 4, 4))
merged_plot <- merged_plot + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
merged_plot

ggplot2::ggsave("~/Local_Documents/CN_Signatures/plots/sigs_final_Vancouver_HRD.png", plot = merged_plot, width = 16, height = 10)
