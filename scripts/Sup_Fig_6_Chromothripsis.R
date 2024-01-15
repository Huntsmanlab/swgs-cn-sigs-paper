### Script to compute the association between signatures and predicted chromothripsis
### Chromothripsis determined via CNV_chromothripsis models

### Required packages
library(dplyr)
library(ggpubr)

#### Read in the data
# Signature exposures
sig_data <- read.delim(file = "./cn_calls/agglomerated_exposures_table.csv", sep = ",", row.names = 1, header = TRUE) %>%
  mutate(sample_id = rownames(.)) %>%
  relocate(sample_id)

cnv_chromo_data_all <- read.delim("./cn_calls/30kb_aCN_CNV_chromothripsis_predictions_all_models.csv", sep = ",", row.names = 1) %>%
  rename(sample_id = sample)

### CNV_chromothripsis model
# Joining the signature exposures with the CNV_chromothripsis predictions
merged_data_cnv_all <- left_join(x = cnv_chromo_data_all, y = sig_data, by = join_by(sample_id))

# Calculate the correlation (Pearson) between each signature and the predicted probability of chromothripsis
cors_all <- sapply(merged_data_cnv_all[, c("VS1", "VS2", "VS3", "VS4", "VS5", "BS1", "BS2", "BS3", "BS4", "BS5", "BS6", "BS7")], function(x) cor.test(x = x, y = merged_data_cnv_all$pred_prob, method = "pearson")$estimate)

p_vals_all <- sapply(merged_data_cnv_all[, c("VS1", "VS2", "VS3", "VS4", "VS5", "BS1", "BS2", "BS3", "BS4", "BS5", "BS6", "BS7")], function(x) cor.test(x = x, y = merged_data_cnv_all$pred_prob, method = "pearson")$p.value)

p_vals_adj_all <- p.adjust(p = p_vals_all, method = "holm")
p_vals_adj_all <- signif(x = p_vals_adj_all, digits = 2)

# Examine the correlation between each Vancouver signature and the predicted probability of chromothripsis (all models); plot with adjusted p-value
vs1_sp_all <- ggscatter(data = merged_data_cnv_all, x = "VS1", y = "pred_prob", add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, xlim = c(0, 1), ylim = c(0, 1)) +
  stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p_vals_adj_all["VS1"])), method = "pearson", label.x = 0.5, label.y = 1)
vs2_sp_all <- ggscatter(data = merged_data_cnv_all, x = "VS2", y = "pred_prob", add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, xlim = c(0, 1), ylim = c(0, 1)) +
  stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p_vals_adj_all["VS2"])), method = "pearson", label.x = 0.5, label.y = 1)
vs3_sp_all <- ggscatter(data = merged_data_cnv_all, x = "VS3", y = "pred_prob", add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, xlim = c(0, 1), ylim = c(0, 1)) +
  stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p_vals_adj_all["VS3"])), method = "pearson", label.x = 0.5, label.y = 1)
vs4_sp_all <- ggscatter(data = merged_data_cnv_all, x = "VS4", y = "pred_prob", add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, xlim = c(0, 1), ylim = c(0, 1)) +
  stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p_vals_adj_all["VS4"])), method = "pearson", label.x = 0.5, label.y = 1)
vs5_sp_all <- ggscatter(data = merged_data_cnv_all, x = "VS5", y = "pred_prob", add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, xlim = c(0, 1), ylim = c(0, 1)) +
  stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p_vals_adj_all["VS5"])), method = "pearson", label.x = 0.5, label.y = 1)

vs_sig_chromo_plot_all <- cowplot::plot_grid(vs1_sp_all, vs2_sp_all, vs3_sp_all, vs4_sp_all, vs5_sp_all, ncol = 2, labels = "A")
ggsave(plot = vs_sig_chromo_plot_all, filename = "van_sigs_correlation_chromothripsis_all_models.png", path = "../Output", device = "png", height = 14, width = 9)

# Examine the correlation between each Brenton signature and the predicted probability of chromothripsis (all models)
bs1_sp_all <- ggscatter(data = merged_data_cnv_all, x = "BS1", y = "pred_prob", add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, xlim = c(0, 1), ylim = c(0, 1)) +
  stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p_vals_adj_all["BS1"])), method = "pearson", label.x = 0.5, label.y = 1)
bs2_sp_all <- ggscatter(data = merged_data_cnv_all, x = "BS2", y = "pred_prob", add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, xlim = c(0, 1), ylim = c(0, 1)) +
  stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p_vals_adj_all["BS2"])), method = "pearson", label.x = 0.5, label.y = 1)
bs3_sp_all <- ggscatter(data = merged_data_cnv_all, x = "BS3", y = "pred_prob", add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, xlim = c(0, 1), ylim = c(0, 1)) +
  stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p_vals_adj_all["BS3"])), method = "pearson", label.x = 0.5, label.y = 1)
bs4_sp_all <- ggscatter(data = merged_data_cnv_all, x = "BS4", y = "pred_prob", add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, xlim = c(0, 1), ylim = c(0, 1)) +
  stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p_vals_adj_all["BS4"])), method = "pearson", label.x = 0.5, label.y = 1)
bs5_sp_all <- ggscatter(data = merged_data_cnv_all, x = "BS5", y = "pred_prob", add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, xlim = c(0, 1), ylim = c(0, 1)) +
  stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p_vals_adj_all["BS5"])), method = "pearson", label.x = 0.5, label.y = 1)
bs6_sp_all <- ggscatter(data = merged_data_cnv_all, x = "BS6", y = "pred_prob", add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, xlim = c(0, 1), ylim = c(0, 1)) +
  stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p_vals_adj_all["BS6"])), method = "pearson", label.x = 0.5, label.y = 1)
bs7_sp_all <- ggscatter(data = merged_data_cnv_all, x = "BS7", y = "pred_prob", add = "reg.line", add.params = list(color = "blue"), conf.int = TRUE, xlim = c(0, 1), ylim = c(0, 1)) +
  stat_cor(aes(label = paste0(..r.label.., "~`,`~`p =`~", p_vals_adj_all["BS7"])), method = "pearson", label.x = 0.5, label.y = 1)

bs_sig_chromo_plot_all <- cowplot::plot_grid(bs1_sp_all, bs2_sp_all, bs3_sp_all, bs4_sp_all, bs5_sp_all, bs6_sp_all, bs7_sp_all, ncol = 2, labels = "B")
ggsave(plot = bs_sig_chromo_plot_all, filename = "brenton_sigs_correlation_chromothripsis_all_models.png", path = "../Output", device = "png", height = 14, width = 9)
