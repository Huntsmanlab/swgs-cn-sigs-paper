# Create an oncoplot
# Use maftools

suppressPackageStartupMessages({
  library(maftools)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(readr)
  library(magrittr)
  library(tidyr)
  library(patchwork)
})

# Create MAF file from custom csv for p53abn cn-sigs project
maff <- readxl::read_xls('./mutation_calls/JSB_reclass_finalmutationcalls_p53AbNcohort.xls')

TPanel_list_samples <- readxl::read_xlsx("./fig_1/TPanel_list_samples.xlsx")
TPanel_list_samples <- TPanel_list_samples[-which(TPanel_list_samples$sample_id == 'VOA888A'), ]
samples_of_interest <- data.frame(TPanel_list_samples$sample_id)
colnames(samples_of_interest) <- 'Tumor_Sample_Barcode'

colnames(maff)[1] <- 'sample_id'
colnames(maff)[2] <- 'chromosome'
colnames(maff)[3] <- 'start'
maff$end <- maff$start + nchar(maff$REF) - 1
colnames(maff)[9] <- 'ref'
colnames(maff)[10] <- 'alt'
maff <- maff[,-c(13)]
maff <- maff[,-c(20)]
colnames(maff)[20] <- 'consequence'
colnames(maff)[22] <- 'gene_symbol'

maff <- maff[,c(1:3, 37, 4:36)]

maff <- maff %>% dplyr::rename(Hugo_Symbol = gene_symbol,
                               Entrez_Gene_Id = Gene_ID,
                               Chromosome = chromosome,
                               Start_Position = start,
                               End_Position = end,
                               Variant_Classification = consequence,
                               Variant_Type = Oncoplot_Variant_Type,
                               Reference_Allele = ref,
                               Tumor_Seq_Allele2 = alt,
                               Tumor_Sample_Barcode = sample_id,
                               n_depth = DP)
maff$Center <- 'GSC'
maff$NCBI_Build <- 37
maff$Strand <- '+'

maff <- samples_of_interest %>% left_join(maff)
maff$Tumor_Sample_Barcode[which(maff$Tumor_Sample_Barcode == 'VOA923A')] <- 'YW-EC002'
maff$Variant_Classification[which(maff$Variant_Classification == 'Unknown')] <- 'InDel'

maff = read.maf(maf = maff, vc_nonSyn = c("5'UTR", "Frame_Shift_Del",
                                          "Frame_Shift_Ins", "In_Frame_Del",
                                          "In_Frame_Ins", "Missense_Mutation",
                                          "Nonsense_Mutation", "Nonstop_Mutation",
                                          "Splice_Site", "Translation_Start_Site",
                                          "InDel"))

CCNE1 <- read.csv("./cn_calls/CCNE1.csv")
HER2 <- read.csv("./cn_calls/HER2.csv")
cancer_type <- read.csv("./clinical_data/cancer_type_rescored_histotypes_187.csv")
HER2_Track <- readxl::read_xlsx("./clinical_data/HER2_Track_simplified.xlsx")
ploidy <- read.csv("./cn_calls/ploidy_classification_187.csv")
survival_data <- read.csv("./clinical_data/survival_track_Three_Year_PFS_187_rescored.csv")
tumour_stage <- read.csv("./clinical_data/tumour_stage_187.csv")

colnames(ploidy)[1] <- 'sample'
colnames(cancer_type)[1] <- 'sample'
colnames(HER2_Track)[1] <- 'sample'
colnames(survival_data)[1] <- 'sample'
colnames(tumour_stage)[1] <- 'sample'
gene_tracks <- rbind(HER2, CCNE1)

pdf(file = "~/Local_Documents/CN_Signatures/plots/Oncoplot.pdf")
pp <- oncoplot(maf = maff, draw_titv = FALSE, top = 25, drawColBar = FALSE)
dev.off()

oncoplot_samples <- data.frame(pp)
colnames(oncoplot_samples)[1] <- 'sample'

cancer_type <- oncoplot_samples %>% left_join(cancer_type)
cancer_type$cancer_type[which(cancer_type$sample == 'CC-RJH0177')] <- 'serous'
cancer_type$cancer_type[which(cancer_type$sample == 'YW-EC042')] <- 'clear cell'
cancer_type$cancer_type[which(cancer_type$sample == 'CC-LAV-0649')] <- 'serous'
cancer_type$cancer_type[which(cancer_type$sample == 'CC-NSH-0337')] <- 'endometrioid'

HER2_Track$sample[which(HER2_Track$sample == 'YW-EC002(VOA923A)')] <- 'YW-EC002'
HER2_Track <- oncoplot_samples %>% left_join(HER2_Track)
HER2_Track$HER2_IHC_score[which(HER2_Track$HER2_IHC_score == 'NA')] <- 'No Data'

gene_tracks$sample <- gsub("\\.", "-", gene_tracks$sample)
gene_tracks <- oncoplot_samples %>% left_join(gene_tracks)

gene_tracks[21, ] <- c("CC-RJH0177", "No Data", "HER2")
gene_tracks[74, ] <- c("YW-EC042", "No Data", "HER2")
gene_tracks[249, ] <- c("CC-LAV-0649", "No Data", "HER2")
gene_tracks[270, ] <- c("CC-NSH-0337", "No Data", "HER2")

gene_tracks[nrow(gene_tracks) + 1,] <- c("CC-RJH0177", "No Data", "CCNE1")
gene_tracks[nrow(gene_tracks) + 1,] <- c("YW-EC042", "No Data", "CCNE1")
gene_tracks[nrow(gene_tracks) + 1,] <- c("CC-LAV-0649", "No Data", "CCNE1")
gene_tracks[nrow(gene_tracks) + 1,] <- c("CC-NSH-0337", "No Data", "CCNE1")

ploidy <- oncoplot_samples %>% left_join(ploidy)
ploidy$classification[which(is.na(ploidy$classification))] <- 'No Data'

survival_data$sample <- gsub("\\.", "-", survival_data$sample)
survival_data <- oncoplot_samples %>% left_join(survival_data)

survival_data[11, ] <- c("CC-RJH0177", "No Data", "Three Year PFS")
survival_data[38, ] <- c("YW-EC042", "No Data", "Three Year PFS")
survival_data[126, ] <- c("CC-LAV-0649", "No Data", "Three Year PFS")
survival_data[137, ] <- c("CC-NSH-0337", "No Data", "Three Year PFS")

tumour_stage$sample <- gsub("\\.", "-", tumour_stage$sample)
tumour_stage <- oncoplot_samples %>% left_join(tumour_stage)
tumour_stage$stage_full[11] <- "No Data"
tumour_stage$stage_full[38] <- "No Data"
tumour_stage$stage_full[126] <- "No Data"
tumour_stage$stage_full[137] <- "No Data"

cancer_type$sample <- factor(cancer_type$sample, pp)
HER2_Track$sample <- factor(HER2_Track$sample, pp)
gene_tracks$sample <- factor(gene_tracks$sample, pp)
ploidy$sample <- factor(ploidy$sample, pp)
survival_data$sample <- factor(survival_data$sample, pp)
tumour_stage$sample <- factor(tumour_stage$sample, pp)

t <- ggplot(ploidy) +
  geom_tile(aes(x = sample, y = 30, fill = classification), width = 1, height = 0.2) +
  scale_fill_manual(
    values = c(
      "1" = "red", "2" = "grey", "3" = "lightblue", "4" = "blue", "5+" = "black", "No Data" = "white"
    )
  ) +
  theme(plot.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 18, angle = 0, vjust = 0.5, hjust = 1, margin = margin(r = -150)),
        legend.position = "none") +
  labs(fill = '', x = "", y = "")

q <- ggplot(HER2_Track) +
  geom_tile(aes(x = sample, y = 30, fill = HER2_IHC_score), width = 1, height = 0.2) +
  scale_fill_manual(values = c("black", "chartreuse3", "lightblue", "darkgoldenrod2", "white")) +
  theme(plot.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 18, angle = 0, vjust = 0.5, hjust = 1, margin = margin(r = -130)),
        legend.position = "none") +
  labs(fill = '', x = "", y = "")

w <- ggplot(cancer_type) +
  geom_tile(aes(x = sample, y = 25, fill = str_wrap(cancer_type, 20)), width = 1, height = 0.2) +
  scale_fill_manual(values = c(
    "carcinosarcoma" = "black",
    "clear cell" = "chartreuse3",
    "endometrioid" = "darkgoldenrod2",
    "other" = "red",
    "serous" = "lightblue")) +
  theme(plot.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 18, angle = 0, vjust = 0.5, hjust = 1, margin = margin(r = -105)),
        legend.position = "none") +
  labs(fill = '', x = "", y = "")

y <- ggplot(tumour_stage) +
  geom_tile(aes(x = sample, y = 30, fill = stage_full), width = 1, height = 0.2) +
  scale_fill_manual(values = c("black", "green", "purple", "lightblue", "white")) +
  theme(plot.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 18, angle = 0, vjust = 0.5, hjust = 1, margin = margin(r = -90)),
        legend.position = "none") +
  labs(fill = '', x = "", y = "")


gene_plot <- ggplot2::ggplot(gene_tracks, ggplot2::aes(sample, gene, fill=classification)) +
  ggplot2::geom_tile() +
  scale_fill_manual(
    values = c(
      "Gain" = "lightblue", "Amplification" = "blue", "No Data" = "white", "Normal" = "grey", "Loss" = "red", "High Amplification" = "black"
    ), breaks = c('High Amplification', 'Amplification', 'Gain', 'Normal', 'Loss', 'No Data')
  ) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 18),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_text(size = 18),
                 axis.text.y = ggplot2::element_text(size = 15),
                 legend.position = "none") +
  ggplot2::labs(fill = 'Status', x = "", y = " ") +
  ggplot2::ggtitle("")

survival_track <- ggplot2::ggplot(survival_data, ggplot2::aes(sample, year, fill=classification)) +
  ggplot2::geom_tile() +
  scale_fill_manual(
    values = c(
      "No Event" = "aquamarine3", "No Data" = "white", "Any Event" = "brown1"
    )
  ) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 18),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_text(size = 18),
                 axis.text.y = ggplot2::element_text(size = 15),
                 legend.position = 'none') +
  ggplot2::labs(fill = 'Status', x = "", y = " ") +
  ggplot2::ggtitle("")


merged_plot <- w + y + survival_track + q + gene_plot + t +
  plot_layout(nrow = 6, ncol = 1,
              widths = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
              heights = c(0.5, 0.5, 0.5, 0.5, 1, 0.5),
              guides = 'auto')
merged_plot <- merged_plot + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
merged_plot

ggplot2::ggsave("~/Local_Documents/CN_Signatures/plots/oncoplot_tracks.png", plot = merged_plot, width = 16, height = 10)
