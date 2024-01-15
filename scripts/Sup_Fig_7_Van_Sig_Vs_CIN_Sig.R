# Script to compute the signature exposures for the p53abn samples to the "pan-cancer compendium of chromosomal instability" (CIN) signatures and compare to our signatures
# Relevant paper: https://doi.org/10.1038/s41586-022-04789-9

# Required packages
library(CINSignatureQuantification)
library(dplyr)
library(easyalluvial)
library(ggplot2)

# Read in the copy-number data and existing signature exposures
acn_data <- readRDS(file = "~/Downloads/30kb_aCN_comCNVfilt_187_filter_false.rds")
van_exp <- read.delim(file = "~/Downloads/agglomerated_exposures_table.csv", sep = ",", header = TRUE)

data("Drews2022_TCGA_Signatures")
cin_sig_by_comp <- Drews2022_TCGA_Signatures

# Format the tables neatly and as required by the feature extraction and signature exposure steps
acn_data_join <- acn_data %>%
  bind_rows(.id = "sample") %>%
  relocate(sample, .after = last_col())

van_exp <- van_exp %>%
  rename("sample" = "X")

# Quantify signature exposures; two samples will be dropped due to low segment counts (CC.JGH.0469, CC.JGH.0487)
cin_exp <- quantifyCNSignatures(object = acn_data_join, method = "drews", build = "hg19")
cin_act <- data.frame(getActivities(object = cin_exp, type = "threshold")) %>%
  mutate(sample = rownames(.))

# Add our new signature exposures to the existing exposure table
all_exp <- full_join(x = van_exp, y = cin_act, by = join_by(sample)) %>%
  mutate(max_cin_sig = names(x = .[, 16:32])[max.col(m = .[, 16:32], ties.method = "first")])

# Save table
write.table(x = all_exp, file = "../Output/Reviewers/all_sig_exposures.csv", sep = ",", row.names = FALSE, col.names = TRUE)

# Generate an alluvial plot, displaying the correspondence between the Vancouver and CIN signatures
plot_data <- all_exp %>%
  filter(!is.na(max_cin_sig)) %>%
  select(max_van_sig, max_cin_sig)

van_vs_cin_plot <- alluvial_wide(data = plot_data, fill_by = "first_variable", stratum_label_size = 3.5) +
  ggplot2::labs(y = "Vancouver Cohort Samples") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        title = element_blank())

# An activities plot, showing per sample exposure distributions
plotActivities(object = cin_exp, type = "threshold")

# Save plot(s)
ggsave(filename = "alluvial_van_vs_cin_sig_exp.png", plot = van_vs_cin_plot, path = "../Output/Reviewers")

