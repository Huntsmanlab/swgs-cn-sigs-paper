# survival.R
#
# Survival analysis of sWGS CN-Sigs samples grouped by predicted signature groups.
#

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(data.table)
  library(purrr)
  library(readr)
  library(survival)
  library(ranger)
  library(ggplot2)
  library(ggfortify)
  library(tibble)
})


###########################
### Functions
###########################

survival_data <- fread(file = '~/Downloads/survival_data_3year.csv', sep = ',')
exposures <- fread(file = '~/Documents/projects/cn_sigs_swgs/paper_data_and_figures/signatureexposures/signature_exposures_component_models_britroc_aCNs_component_by_signature_britroc_aCNs_2023-06-28.csv', sep = ',')
cexps <- fread(file = '~/Documents/projects/cn_sigs_swgs/paper_data_and_figures/signatureexposures/signature_exposures_component_models_britroc_aCNs_component_by_signature_britmodelsvansigs5_aCNs_2023-06-26.csv', sep = ',')
metadata <- data.table::fread(file = '~/Documents/projects/cn_sigs_swgs/metadata/metadata.tsv', sep = '\t', header = TRUE)
gene_info <- data.table::fread(file = '~/Documents/projects/cn_sigs_swgs/paper_data_and_figures/plot_annotation_tracks/long_gene_tracks_188.csv', sep = ',')
tgene_info <- data.table::fread(file = '~/Documents/projects/cn_sigs_swgs/paper_data_and_figures/plot_annotation_tracks/long_gene_tracks_188_targeted_panel_seq.csv', sep = ',')
p53_info <- data.table::fread(file = '~/Documents/projects/cn_sigs_swgs/paper_data_and_figures/plot_annotation_tracks/p53abn_track_simplified.csv', sep = ',')
shrd <- data.table::fread(file = '~/Documents/projects/cn_sigs_swgs/paper_data_and_figures/plot_annotation_tracks/30kb_combined_Xchr_shallowHRD.tsv',
                          sep = '\t', header = TRUE)

# Progression-free-survival data corrections/calculations

# Most of the pfs_yrs values should be from the last follow up date (these can be follow-ups or deaths)
survival_data$pfs_yrs <- as.numeric((strptime(survival_data$date_lfu, format = "%y-%m-%d")-
                          strptime(survival_data$date_dx, format = "%y-%m-%d"))/31557600)
# However, sometimes there are progression events, and these super-cede any LFU date (death or otherwise)
survival_data <- survival_data %>% dplyr::mutate(pfs_yrs = if_else(date_prog == '-',
                                                                   pfs_yrs,
                                                                   as.numeric((strptime(date_prog, format = "%y-%m-%d")-
                                                                   strptime(date_dx, format = "%y-%m-%d"))/31557600)),
                                                 pfs_sts = if_else(survival_status == 'dead/disease',
                                                                   'pfs.event',
                                                                   if_else(date_prog == '-',
                                                                           'pfs.censored',
                                                                           'pfs.event')))


# survival_data$pfs_yrs <- as.numeric((strptime(survival_data$date_prog, format = "%y-%m-%d")-
#                           strptime(survival_data$date_dx, format = "%y-%m-%d"))/31557600)


# Basic data-prep for survival analysis
gene_info <- gene_info %>% tidyr::spread(gene, classification)
tgene_info <- tgene_info %>% tidyr::spread(gene, classification)
colnames(gene_info)[2:dim(gene_info)[2]] <- paste0(colnames(gene_info)[2:dim(gene_info)[2]], '_cn')
gene_info$sample <- str_replace_all(gene_info$sample, "-", ".")
tgene_info$sample_id <- str_replace_all(tgene_info$sample_id, "-", ".")
shrd$sample_id <- str_replace_all(shrd$sample_id, "-", ".")
p53_info$sample_id <- str_replace_all(p53_info$sample_id, "-", ".")
p53_info$sample_id[192] <- 'YW.EC002'

exposures[,1] <- NULL
exposures <- as.data.frame(t(exposures))
colnames(exposures) <- paste0('S', 1:dim(exposures)[2])
exposures$sample_id <- rownames(exposures)
# PlotSignatureExposures(exposures)
cexps[,1] <- NULL
cexps <- as.data.frame(t(cexps))
colnames(cexps) <- paste0('cS', 1:dim(cexps)[2])
cexps$sample_id <- rownames(cexps)

metadata$sample_id <- str_replace_all(metadata$sample_id, "-", ".")
survival_data$sample_id <- str_replace_all(survival_data$sample_id, "-", ".")
survival_data$agebins <- '80+'
survival_data$agebins[(survival_data$age_surg >= 70) & (survival_data$age_surg < 80)] <- '70-80'
survival_data$agebins[(survival_data$age_surg >= 60) & (survival_data$age_surg < 70)] <- '60-70'
survival_data$agebins[survival_data$age_surg < 60] <- 'under 60'
survival_data$binaryAge <- ifelse(survival_data$age_surg >= 65, "65+", "under 65")
survival_data$os_bin <- ifelse(survival_data$os_sts == "os.event", 1, 0)
survival_data$dss_bin <- ifelse(survival_data$dss_sts == "dss.event", 1, 0)


# survival analysis on the p53_abn subset
survival_data <- survival_data[survival_data$sample_id %in% exposures$sample_id, ]
survival_data_p53 <- survival_data %>% dplyr::left_join(exposures, by = c('sample_id'))
survival_data_p53 <- survival_data_p53 %>% dplyr::left_join(cexps, by = c('sample_id'))
survival_data_p53 <- survival_data_p53 %>% dplyr::left_join(gene_info, by = c('sample_id' = 'sample'))
survival_data_p53 <- survival_data_p53 %>% dplyr::left_join(tgene_info, by = c('sample_id'))
survival_data_p53 <- survival_data_p53 %>% dplyr::left_join(p53_info, by = c('sample_id'))
survival_data_p53 <- survival_data_p53 %>% dplyr::left_join(shrd, by = c('sample_id'))

# Make some new columns
survival_data_p53$max_bsig <- apply(exposures[,c(1:7)], 1, function(x) which.max(x))
survival_data_p53$max_vsig <- apply(cexps[,c(1:5)], 1, function(x) which.max(x))
survival_data_p53$os_yrs <- as.numeric(survival_data_p53$os_yrs)
survival_data_p53$dss_yrs <- as.numeric(survival_data_p53$dss_yrs)
survival_data_p53$age_dx <- as.numeric(survival_data_p53$age_dx)
survival_data_p53$hrd <- apply(survival_data_p53[,c('ATM', 'ATRX', 'BRCA1', 'BRCA2', 'BAP1',
                                                    'BLM', 'BRIP1', 'EXO1', 'FANCA', 'FANCC',
                                                    'FANCD2', 'FANCE', 'FANCG',
                                                    'PALB2', 'RAD50', 'RAD51')], 1,
                               function(x) any(x == 'Present'))
table(survival_data_p53$hrd)



# Make 3 yr and 4 yr pfs survival columns
survival_data_p53$pfs_3yrs <- survival_data_p53$pfs_yrs
survival_data_p53$pfs_4yrs <- survival_data_p53$pfs_yrs
survival_data_p53$pfs_3sts <- survival_data_p53$pfs_sts
survival_data_p53$pfs_4sts <- survival_data_p53$pfs_sts

survival_data_p53$pfs_3yrs[survival_data_p53$pfs_yrs > 3] <- 3
survival_data_p53$pfs_3sts[survival_data_p53$pfs_yrs > 3] <- 'pfs.censored'
survival_data_p53$pfs_4yrs[survival_data_p53$pfs_yrs > 4] <- 4
survival_data_p53$pfs_4sts[survival_data_p53$pfs_yrs > 4] <- 'pfs.censored'

write.table(survival_data_p53, file = '~/Downloads/x.tsv', quote = FALSE,
            sep = '\t', row.names = FALSE, col.names = TRUE)

# survival_data_p53 <- survival_data %>% dplyr::filter(status == 'p53_abn')
# survival_data_p53 <- survival_data_p53 %>% dplyr::filter(survival_data_p53$sample_id %in% colnames(exposures))
# sig_groups <- data.frame(sample_id = colnames(exposures), max_sigs = as.character(apply(exposures, 2, function(x) which.max(x))) )
# survival_data_p53 <- survival_data_p53 %>% dplyr::left_join(sig_groups, by = c('sample_id'))
# survival_data_p53 <- survival_data_p53 %>%
#                       dplyr::left_join(tibble::rownames_to_column(as.data.frame(t(exposures)), "sample_id"), by = c('sample_id'))

### KM Plots
### Removing the signature classes that don't have enough observations for the km plot to make sense
survival_data_in <- survival_data_p53 %>% dplyr::filter(max_bsig %in% c('1','3','4','5','7'))

km_trt_fit <- survfit(Surv(dss_yrs, dss_bin) ~ max_bsig, data=survival_data_in)
p <- autoplot(km_trt_fit, conf.int = F, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by max. custom sig. exp.") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
p
pp <- autoplot(km_trt_fit, conf.int = T, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by max. custom sig. exp.") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
pp
ggsave(filename = '~/Documents/projects/cn_signatures_shallowWGS/plotting/survival_plots/p53_abn/kapMeier_QDNA_30kb_autosomesOnly_p53_abn_MADabsbritrocSignature_exposures_3.png',
       plot = p, width = 8, height = 8, units = 'in')

### Tests of Statistical Significance - Accelerated Failure Time Models
# https://www.drizopoulos.com/courses/emc/survival%20analysis%20in%20r%20companion
logrank <- survdiff(Surv(os_yrs, os_bin) ~ max_sigs, data = survival_data_in, subset = (survival_data_in$max_sigs %in% c(1,4)))
logrank <- survdiff(Surv(os_yrs, os_bin) ~ cS2 + cS3 + cS4 + cS5 + binaryAge, data = survival_data_p53)

# Cox Proportional hazards
cox <- coxph(Surv(os_yrs, os_bin) ~ max_sigs + hist + binaryAge + grade, data = survival_data_in)
cox <- coxph(Surv(os_yrs, os_bin) ~ S1 + S3 + S4 + S5 + S7 + cS1 + cS2 + cS3 + cS4 + cS5 + binaryAge, data = survival_data_p53)
cox <- coxph(Surv(dss_yrs, dss_bin) ~ S1 + S3 + S4 + S5 + S7 + age_dx + nLGAs, data = survival_data_p53)
cox <- coxph(Surv(dss_yrs, dss_bin) ~ cS2 + cS3 + cS4 + cS5 + age_dx, data = survival_data_p53_3yr)
cox <- coxph(Surv(dss_yrs, dss_bin) ~ BRCA1_cn + age_dx, data = survival_data_p53)
summary(cox)
cox_fit <- survfit(cox)
#plot(cox_fit, main = "cph model", xlab="Days")
autoplot(cox_fit)






# metadata <- metadata %>% dplyr::filter(str_sub(sample_id,-1,-1) != 'N')
# metadata$sample_id[str_sub(metadata$sample_id,-2,-1) == '-T'] <- str_sub(metadata$sample_id,1,11)[str_sub(metadata$sample_id,-2,-1) == '-T']
# metadata$sample_id[str_sub(metadata$sample_id,-2,-1) == '-T'] <- str_sub(metadata$sample_id,1,8)[str_sub(metadata$sample_id,-2,-1) == '-T']
# metadata$sample_id[str_sub(metadata$sample_id,-1,-1) == 'T'] <- str_sub(metadata$sample_id,1,11)[str_sub(metadata$sample_id,-1,-1) == 'T']
# metadata$sample_id[25] <- 'VOA2050'
#
# long_data <- gather(exposures)
# long_data$max_sig <- rep(apply(exposures, 2, function(x) which.max(x)), times = 1, each = nsigs)
# long_data$sigs <- rep(1:nsigs,dim(exposures)[2])
# colnames(long_data) <- c('X', 'Z', 'max_sig', 'Y')
# long_data <- long_data %>% arrange(max_sig)
# long_data$X <- factor(long_data$X, levels = unique(long_data$X))


#######################
### CCNE1 Investigation
#######################
logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ CCNE1_cn, data = survival_data_p53_3yr)
logrank
survD_ccne1 <- survival_data_p53_3yr %>% dplyr::filter(!(CCNE1_cn %in% c('Loss', 'No Data')))
logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ CCNE1_cn, data = survD_ccne1)
logrank

km_trt_fit <- survfit(Surv(dss_yrs, dss_bin) ~ CCNE1_cn, data=survD_ccne1)
p <- autoplot(km_trt_fit, conf.int = F, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "Kaplan Meier survival curves for patients stratified by CCNE1 Copy-Number") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
p
pp <- autoplot(km_trt_fit, conf.int = TRUE, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "Kaplan Meier survival curves for patients stratified by CCNE1 Copy-Number") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
pp

#######################
### HER2 Investigation
#######################
logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ HER2_cn, data = survival_data_p53_3yr)
logrank
survD_her2 <- survival_data_p53_3yr %>% dplyr::filter(!(HER2_cn %in% c('No Data')))
logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ HER2_cn, data = survD_her2)
logrank

km_trt_fit <- survfit(Surv(dss_yrs, dss_bin) ~ HER2_cn, data=survD_her2)
p <- autoplot(km_trt_fit, conf.int = F, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "Kaplan Meier survival curves for patients stratified by HER2 Copy-Number") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
p
pp <- autoplot(km_trt_fit, conf.int = TRUE, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "Kaplan Meier survival curves for patients stratified by HER2 Copy-Number") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
pp


#######################
### FBXW7 Investigation
#######################
logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ CCNE1_cn, data = survival_data_p53_3yr)
logrank
survD_ccne1 <- survival_data_p53_3yr %>% dplyr::filter(!(CCNE1_cn %in% c('Loss', 'No Data')))
logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ CCNE1_cn, data = survD_ccne1)
logrank

km_trt_fit <- survfit(Surv(dss_yrs, dss_bin) ~ CCNE1_cn, data=survD_ccne1)
p <- autoplot(km_trt_fit, conf.int = F, xlab = 'time (years)') +
  labs(title = "Kaplan Meier survival curves for patients stratified by CCNE1 Copy-Number") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
p
pp <- autoplot(km_trt_fit, conf.int = TRUE, xlab = 'time (years)') +
  labs(title = "Kaplan Meier survival curves for patients stratified by CCNE1 Copy-Number") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
pp


#######################
### BRCA1 Investigation
#######################
logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ BRCA1_cn, data = survival_data_p53_3yr)
logrank
survD_brca1 <- survival_data_p53_3yr %>% dplyr::filter(!(BRCA1_cn %in% c('No Data')))
logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ BRCA1_cn, data = survD_brca1)
logrank

km_trt_fit <- survfit(Surv(dss_yrs, dss_bin) ~ BRCA1_cn, data=survD_brca1)
p <- autoplot(km_trt_fit, conf.int = F, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "Kaplan Meier survival curves for patients stratified by BRCA1 Copy-Number") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
p
pp <- autoplot(km_trt_fit, conf.int = TRUE, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "Kaplan Meier survival curves for patients stratified by BRCA1 Copy-Number") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
pp


#######################
### BRCA2 Investigation
#######################
logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ BRCA2_cn, data = survival_data_p53_3yr)
logrank
survD_brca2 <- survival_data_p53_3yr %>% dplyr::filter(!(BRCA2_cn %in% c('Amplification')))
logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ BRCA2_cn, data = survD_brca2)
logrank

km_trt_fit <- survfit(Surv(dss_yrs, dss_bin) ~ BRCA2_cn, data=survD_brca2)
p <- autoplot(km_trt_fit, conf.int = F, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "Kaplan Meier survival curves for patients stratified by BRCA2 Copy-Number") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
p
pp <- autoplot(km_trt_fit, conf.int = TRUE, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "Kaplan Meier survival curves for patients stratified by BRCA2 Copy-Number") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
pp


#######################
### HRD Investigation
#######################
logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ hrd, data = survival_data_p53_3yr)
logrank

logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ hrd_status, data = survival_data_p53_3yr)
logrank

survival_data_in <- survival_data_p53_3yr
survival_data_in$hrd_status[survival_data_p53_3yr$hrd_status == 'Borderline [15;19]'] <- 'Yes (>= 20)'
logrank <- survdiff(Surv(dss_yrs, dss_bin) ~ hrd_status, data = survival_data_in)
logrank

km_trt_fit <- survfit(Surv(dss_yrs, dss_bin) ~ BRCA2_cn, data=survD_brca2)
p <- autoplot(km_trt_fit, conf.int = F, xlab = 'time (years)') +
  labs(title = "Kaplan Meier survival curves for patients stratified by BRCA2 Copy-Number") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
p
pp <- autoplot(km_trt_fit, conf.int = TRUE, xlab = 'time (years)') +
  labs(title = "Kaplan Meier survival curves for patients stratified by BRCA2 Copy-Number") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
pp


#######################
### PCA of Signatures
#######################
library(factoextra)

# Our signatures
pc1 <- prcomp(t(survival_data_p53_3yr[,40:44]),
              center = TRUE,
              scale. = F)
coord <- t(t(pc1$rotation) * pc1$sdev)
contrib <- t(t(coord ^ 2) / colSums(coord ^ 2)) * 100
var <- get_pca_var(pc1)
fviz_contrib(pc1, choice = "ind", axes = 1,
             ggtheme = theme_classic(), fill = 'black', sort.val = 'none')
summary(pc1)
temp <- pc1$rotation[,c(1:4)]
colnames(temp) <- c('VPC1', 'VPC2', 'VPC3', 'VPC4')
survival_data_p53_3yr <- cbind(survival_data_p53_3yr, temp)

# Brenton signatures
pc2 <- prcomp(t(survival_data_p53_3yr[,33:39]),
              center = TRUE,
              scale. = F)
coord <- t(t(pc2$rotation) * pc2$sdev)
contrib <- t(t(coord ^ 2) / colSums(coord ^ 2)) * 100
var <- get_pca_var(pc2)
fviz_contrib(pc2, choice = "ind", axes = 1,
             ggtheme = theme_classic(), fill = 'black', sort.val = 'none')
summary(pc2)
temp <- pc2$rotation[,c(1:5)]
colnames(temp) <- c('BPC1', 'BPC2', 'BPC3', 'BPC4', 'BPC5')
survival_data_p53_3yr <- cbind(survival_data_p53_3yr, temp)


cox <- coxph(Surv(dss_yrs, dss_bin) ~ VPC1 + VPC2 + VPC3 + VPC4 + age_dx, data = survival_data_p53_3yr)
summary(cox)

cox <- coxph(Surv(dss_yrs, dss_bin) ~ BPC1 + BPC2 + BPC3 + BPC4 + BPC5 + BRCA2_cn + age_dx, data = survival_data_p53_3yr)
summary(cox)

ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity")

