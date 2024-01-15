# survival_analysis.R
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
  library(factoextra)
  library(survminer)
})


###########################
### Prep Data
###########################

#### Read in data
survival_data <- fread(file = '~/Downloads/survival_data.csv', sep = ',', header = T,
                       colClasses=c(date_dx='character',
                                    date_lfu='character',
                                    date_surg='character',
                                    date_prog='character'))
exposures <- fread(file = '~/Downloads/agglomerated_exposures_table.csv', sep = ',')
histotypes <- fread(file = '~/Downloads/cancer_type_rescored_histotypes_187.csv', sep = ',')
ploidy <- fread(file = '~/Downloads/ploidy_classification_187.csv', sep = ',')


#### Calculate survival times from dates and call survival event status

# Progression-free-survival data corrections/calculations
# Most of the pfs_yrs values should be from the last follow up date (these can be follow-ups or deaths)
survival_data$pfs_yrs <- as.numeric((strptime(survival_data$date_lfu, format = "%Y-%m-%d")-
                                       strptime(survival_data$date_dx, format = "%Y-%m-%d"))/31557600)
# However, sometimes there are progression events, and these super-cede any LFU date (death or otherwise)
# Note I am using time in days as the divisor here because when NA values are present strptime
# returns time in days and not seconds
survival_data <- survival_data %>% dplyr::mutate(pfs_yrs = if_else(date_prog == '-',
                                                                   pfs_yrs,
                                                                   as.numeric((strptime(date_prog, format = "%Y-%m-%d")-
                                                                                 strptime(date_dx, format = "%Y-%m-%d"))/365.2422)),
                                                 pfs_sts = if_else(status == 'dead/disease',
                                                                   'pfs.event',
                                                                   if_else(date_prog == '-',
                                                                           'pfs.censored',
                                                                           'pfs.event')))
# Disease specific-survival data corrections/calculations
# All of the dss_yrs values should be from the last follow up date
# These can be follow-ups or deaths, in many cases the value will just be censored
survival_data$dss_yrs <- as.numeric((strptime(survival_data$date_lfu, format = "%Y-%m-%d")-
                                       strptime(survival_data$date_dx, format = "%Y-%m-%d"))/31557600)
survival_data <- survival_data %>% dplyr::mutate(dss_sts = if_else(status == 'dead/disease',
                                                                   'dss.event',
                                                                   'dss.censored'))


#### Fix-up and format data for joining
colnames(survival_data) <- c('sample_id', colnames(survival_data)[2:length(colnames(survival_data))])
survival_data$sample_id[191] <- 'YW.EC002'
survival_data$sample_id <- str_replace_all(survival_data$sample_id, "-", ".")
colnames(exposures) <- c('sample_id', colnames(exposures)[2:length(colnames(exposures))])
histotypes$sample_id <- str_replace_all(histotypes$sample_id, "-", ".")
#exposures <- exposures %>% dplyr::filter(cohort == 'VanP53ABNendo')
ploidy <- ploidy %>%
  mutate_at(.vars = c("sample"), .funs = gsub, pattern = "-", replacement = ".") %>%
  rename(sample_id = sample, ploidy = classification)

# Brief aside to check for what the overall survival for our cohort looks like
mean(survival_data$pfs_yrs[!(survival_data$pfs_yrs > 3)])
mean(survival_data$dss_yrs[!(survival_data$dss_yrs > 3)])

# joins
survival_data <- survival_data[survival_data$sample_id %in% histotypes$sample_id,]
survival_data_p53 <- survival_data %>% dplyr::left_join(histotypes, by = c('sample_id'))
survival_data_p53 <- survival_data_p53 %>% dplyr::left_join(exposures, by = c('sample_id'))
survival_data_p53 <- survival_data_p53 %>% dplyr::left_join(ploidy, by = "sample_id")


#### Make 3 yr and 4 yr survival columns
survival_data_p53$os_3yrs <- survival_data_p53$os_yrs
survival_data_p53$dss_3yrs <- survival_data_p53$dss_yrs
survival_data_p53$pfs_3yrs <- survival_data_p53$pfs_yrs
survival_data_p53$pfs_4yrs <- survival_data_p53$pfs_yrs
survival_data_p53$os_3sts <- survival_data_p53$os_sts
survival_data_p53$dss_3sts <- survival_data_p53$dss_sts
survival_data_p53$pfs_3sts <- survival_data_p53$pfs_sts
survival_data_p53$pfs_4sts <- survival_data_p53$pfs_sts

survival_data_p53$os_yrs <- as.numeric(survival_data_p53$os_yrs)
survival_data_p53$os_3yrs <- as.numeric(survival_data_p53$os_3yrs)
survival_data_p53$dss_yrs <- as.numeric(survival_data_p53$dss_yrs)
survival_data_p53$dss_3yrs <- as.numeric(survival_data_p53$dss_3yrs)
survival_data_p53$pfs_yrs <- as.numeric(survival_data_p53$pfs_yrs)
survival_data_p53$pfs_3yrs <- as.numeric(survival_data_p53$pfs_3yrs)
survival_data_p53$pfs_4yrs <- as.numeric(survival_data_p53$pfs_4yrs)
survival_data_p53$age_dx <- as.numeric(survival_data_p53$age_dx)

survival_data_p53$os_3yrs[survival_data_p53$os_yrs > 3] <- 3
survival_data_p53$os_3sts[survival_data_p53$os_yrs > 3] <- 'os.censored'
survival_data_p53$dss_3yrs[survival_data_p53$dss_yrs > 3] <- 3
survival_data_p53$dss_3sts[survival_data_p53$dss_yrs > 3] <- 'dss.censored'
survival_data_p53$pfs_3yrs[survival_data_p53$pfs_yrs > 3] <- 3
survival_data_p53$pfs_3sts[survival_data_p53$pfs_yrs > 3] <- 'pfs.censored'
survival_data_p53$pfs_4yrs[survival_data_p53$pfs_yrs > 4] <- 4
survival_data_p53$pfs_4sts[survival_data_p53$pfs_yrs > 4] <- 'pfs.censored'

survival_data_p53$os_bin <- ifelse(survival_data_p53$os_sts == "os.event", 1, 0)
survival_data_p53$os_3bin <- ifelse(survival_data_p53$os_3sts == "os.event", 1, 0)
survival_data_p53$dss_bin <- ifelse(survival_data_p53$dss_sts == "dss.event", 1, 0)
survival_data_p53$dss_3bin <- ifelse(survival_data_p53$dss_3sts == "dss.event", 1, 0)
survival_data_p53$pfs_bin <- ifelse(survival_data_p53$pfs_sts == "pfs.event", 1, 0)
survival_data_p53$pfs_3bin <- ifelse(survival_data_p53$pfs_3sts == "pfs.event", 1, 0)
survival_data_p53$pfs_4bin <- ifelse(survival_data_p53$pfs_4sts == "pfs.event", 1, 0)


####################################
### Kaplan Meier Plots and Logrank Tests
####################################

#### By signature

# Use the first line to remove the signature classes that don't have enough observations, if you want to plot the Brenton signatures; otherwise keep all with line two
# survival_data_in <- survival_data_p53 %>% dplyr::filter(max_brenton_sig %in% c('BS1','BS3','BS4','BS5','BS7'))
survival_data_in <- survival_data_p53

# OS
logrank <- survdiff(Surv(os_3yrs, os_3bin) ~ max_van_sig, data = survival_data_in)
# Significant overall; let's do a pairwise test to see which signatures differ
logrank_pairs <- pairwise_survdiff(Surv(os_3yrs, os_3bin) ~ max_van_sig, data = survival_data_in, p.adjust.method = "BH")

km_trt_fit <- survfit(Surv(os_3yrs, os_3bin) ~ max_van_sig, data=survival_data_in)
os_p <- autoplot(km_trt_fit, conf.int = F, xlab = 'Time (years)', ylab = 'Proportion surviving', ylim = c(0.2, 1), surv.size = 1) +
  labs(title = "KM survival curves for patients grouped by max. sig. exp.", subtitle = paste0("Overall survival, p = ", round(logrank$pvalue, 3))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12),
        plot.subtitle = element_text(size = 14),
        axis.text = element_text(size = 12))
os_p
os_ci_p <- autoplot(km_trt_fit, conf.int = T, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by max. custom sig. exp.", subtitle = paste0("Overall survival, p = ", round(logrank$pvalue, 3))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
os_ci_p

# DSS
logrank <- survdiff(Surv(dss_3yrs, dss_3bin) ~ max_van_sig, data = survival_data_in)

km_trt_fit <- survfit(Surv(dss_3yrs, dss_3bin) ~ max_van_sig, data=survival_data_in)
dss_p <- autoplot(km_trt_fit, conf.int = F, xlab = 'Time (years)', ylab = 'Proportion surviving', ylim = c(0.2, 1), surv.size = 1) +
  labs(title = "KM survival curves for patients grouped by max. custom sig. exp.", subtitle = paste0("Disease-specific survival, p = ", round(logrank$pvalue, 3))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12),
        plot.subtitle = element_text(size = 14),
        axis.text = element_text(size = 12))
dss_p
dss_ci_p <- autoplot(km_trt_fit, conf.int = T, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by max. custom sig. exp.", subtitle = paste0("Disease-specific survival, p = ", round(logrank$pvalue, 3))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
dss_ci_p

# PFS
logrank <- survdiff(Surv(pfs_3yrs, pfs_3bin) ~ max_van_sig, data = survival_data_in)

km_trt_fit <- survfit(Surv(pfs_3yrs, pfs_3bin) ~ max_van_sig, data=survival_data_in)
pfs_p <- autoplot(km_trt_fit, conf.int = F, xlab = 'Time (years)', ylab = 'Proportion surviving', ylim = c(0.2, 1), surv.size = 1) +
  labs(title = "KM survival curves for patients grouped by max. custom sig. exp.", subtitle = paste0("Progression-free survival, p = ", round(logrank$pvalue, 3))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12),
        plot.subtitle = element_text(size = 14),
        axis.text = element_text(size = 12))
pfs_p
pfs_ci_p <- autoplot(km_trt_fit, conf.int = T, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by max. custom sig. exp.", subtitle = paste0("Progression-free survival, p = ", round(logrank$pvalue, 3))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
pfs_ci_p

# Tweak plots for combined grid plot
os_p <- os_p +
  xlab("") +
  ylab("") +
  theme(legend.position = "none")
dss_p <- dss_p +
  labs(title = "") +
  xlab("") +
  theme(legend.position = "none")
pfs_p <- pfs_p +
  labs(title = "") +
  ylab("") +
  theme(legend.position = "bottom", legend.title = element_blank())

survival_fig <- cowplot::plot_grid(os_p, dss_p, pfs_p, ncol = 1)

#### By histotype

#OS
km_trt_fit <- survfit(Surv(os_3yrs, os_3bin) ~ cancer_type, data=survival_data_in)
p <- autoplot(km_trt_fit, conf.int = F, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by histotype", subtitle = "Overall survival") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
p
pp <- autoplot(km_trt_fit, conf.int = T, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by histotype", subtitle = "Overall survival") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
pp

# DSS
km_trt_fit <- survfit(Surv(dss_3yrs, dss_3bin) ~ cancer_type, data=survival_data_in)
p <- autoplot(km_trt_fit, conf.int = F, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by histotype", subtitle = "Disease-specific survival") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
p
pp <- autoplot(km_trt_fit, conf.int = T, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by histotype", subtitle = "Disease-specific survival") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
pp

# PFS
km_trt_fit <- survfit(Surv(pfs_3yrs, pfs_3bin) ~ cancer_type, data=survival_data_in)
p <- autoplot(km_trt_fit, conf.int = F, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by histotype", subtitle = "Progression-free survival") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
p
pp <- autoplot(km_trt_fit, conf.int = T, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by histotype", subtitle = "Progression-free survival") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
pp

#### By ploidy
# OS
logrank <- survdiff(Surv(os_3yrs, os_3bin) ~ ploidy, data = survival_data_in)

km_trt_fit <- survfit(Surv(os_3yrs, os_3bin) ~ ploidy, data=survival_data_in)
os_ploidy <- autoplot(km_trt_fit, conf.int = F, xlab = 'Time (years)', ylab = 'Proportion surviving', ylim = c(0.2, 1), surv.size = 1) +
  labs(title = "KM survival curves for patients grouped by ploidy", subtitle = paste0("Overall survival, p = ", round(logrank$pvalue, 3))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12),
        plot.subtitle = element_text(size = 14),
        axis.text = element_text(size = 12))
os_ploidy
os_ci_ploidy <- autoplot(km_trt_fit, conf.int = T, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by ploidy", subtitle = paste0("Overall survival, p = ", round(logrank$pvalue, 3))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
os_ci_ploidy

# DSS
logrank <- survdiff(Surv(dss_3yrs, dss_3bin) ~ ploidy, data = survival_data_in)

km_trt_fit <- survfit(Surv(dss_3yrs, dss_3bin) ~ ploidy, data=survival_data_in)
dss_ploidy <- autoplot(km_trt_fit, conf.int = F, xlab = 'Time (years)', ylab = 'Proportion surviving', ylim = c(0.2, 1), surv.size = 1) +
  labs(title = "KM survival curves for patients grouped by ploidy", subtitle = paste0("Disease-specific survival, p = ", round(logrank$pvalue, 3))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12),
        plot.subtitle = element_text(size = 14),
        axis.text = element_text(size = 12))
dss_ploidy
dss_ci_ploidy <- autoplot(km_trt_fit, conf.int = T, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by ploidy", subtitle = paste0("Disease-specific survival, p = ", round(logrank$pvalue, 3))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
dss_ci_ploidy

# PFS
logrank <- survdiff(Surv(pfs_3yrs, pfs_3bin) ~ ploidy, data = survival_data_in)

km_trt_fit <- survfit(Surv(pfs_3yrs, pfs_3bin) ~ ploidy, data=survival_data_in)
pfs_ploidy <- autoplot(km_trt_fit, conf.int = F, xlab = 'Time (years)', ylab = 'Proportion surviving', ylim = c(0.2, 1), surv.size = 1) +
  labs(title = "KM survival curves for patients grouped by ploidy", subtitle = paste0("Progression-free survival, p = ", round(logrank$pvalue, 3))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12),
        plot.subtitle = element_text(size = 14),
        axis.text = element_text(size = 12))
pfs_ploidy
pfs_ci_ploidy <- autoplot(km_trt_fit, conf.int = T, xlab = 'time (years)', ylim = c(0.2, 1)) +
  labs(title = "KM survival curves for patients grouped by ploidy", subtitle = paste0("Progression-free survival, p = ", round(logrank$pvalue, 3))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
pfs_ci_ploidy

# Can adjust to test signatures or cancer types
logrank <- survdiff(Surv(os_3yrs, os_3bin) ~ cancer_type, data = survival_data_in)
logrank
logrank <- survdiff(Surv(dss_3yrs, dss_3bin) ~ cancer_type, data = survival_data_in)
logrank
logrank <- survdiff(Surv(pfs_3yrs, pfs_3bin) ~ cancer_type, data = survival_data_in)
logrank

survival_data_in_filtered <- survival_data_in %>% dplyr::filter(!(cancer_type %in% c('other')))
logrank <- survdiff(Surv(os_3yrs, os_3bin) ~ cancer_type, data = survival_data_in_filtered)
logrank
logrank <- survdiff(Surv(dss_3yrs, dss_3bin) ~ cancer_type, data = survival_data_in_filtered)
logrank
logrank <- survdiff(Surv(pfs_3yrs, pfs_3bin) ~ cancer_type, data = survival_data_in_filtered)
logrank



#######################
### PCA of Signatures and Cox Proportional Hazards Model
#######################

# Our signatures
pc1 <- prcomp(t(survival_data_p53[,21:25]),
              center = TRUE,
              scale. = F)
coord <- t(t(pc1$rotation) * pc1$sdev)
contrib <- t(t(coord ^ 2) / colSums(coord ^ 2)) * 100
var <- get_pca_var(pc1)
# Can iterate through this call to get the contribution plot for each signature to each component:
fviz_contrib(pc1, choice = "ind", axes = 1,
             ggtheme = theme_classic(), fill = 'black', sort.val = 'none')
summary(pc1)
temp <- pc1$rotation[,c(1:4)]
colnames(temp) <- c('VPC1', 'VPC2', 'VPC3', 'VPC4')
survival_data_p53 <- cbind(survival_data_p53, temp)

# Brenton signatures
pc2 <- prcomp(t(survival_data_p53[,27:33]),
              center = TRUE,
              scale. = F)
coord <- t(t(pc2$rotation) * pc2$sdev)
contrib <- t(t(coord ^ 2) / colSums(coord ^ 2)) * 100
var <- get_pca_var(pc2)
# Can iterate through this call to get the contribution plot for each signature to each component:
fviz_contrib(pc2, choice = "ind", axes = 1,
             ggtheme = theme_classic(), fill = 'black', sort.val = 'none')
summary(pc2)
temp <- pc2$rotation[,c(1:5)]
colnames(temp) <- c('BPC1', 'BPC2', 'BPC3', 'BPC4', 'BPC5')
survival_data_p53 <- cbind(survival_data_p53, temp)


cox <- coxph(Surv(dss_yrs, dss_bin) ~ VPC1 + VPC2 + VPC3 + VPC4 + age_dx, data = survival_data_p53)
summary(cox)

