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

survival_data <- fread(file = '~/Documents/projects/cn_sigs_swgs/metadata/survival_data.csv', sep = ',')
exposures <- fread(file = '~/repos/cnsignatures/output/batch1-13_output/signature_exposures_QDNA_30kb_autosomesOnly_p53_abn_VAFabsbritrocSignatures.csv', sep = ',')
exposures <- fread(file = '~/repositories/cnsignatures/data/batch1-13_output/signature_exposures_QDNA_30kb_autosomesOnly_p53_abn_absbritrocSignatures.csv', sep = ',')
metadata <- data.table::fread(file = '~/Documents/projects/cn_sigs_swgs/metadata/metadata.tsv', sep = '\t', header = TRUE)
exposures[,1] <- NULL

# Basic data-prep for survival analysis
metadata$sample_id <- str_replace_all(metadata$sample_id, "-", ".")
survival_data$sample_id <- str_replace_all(survival_data$sample_id, "-", ".")
survival_data$agebins <- '80+'
survival_data$agebins[(survival_data$age_surg >= 70) & (survival_data$age_surg < 80)] <- '70-80'
survival_data$agebins[(survival_data$age_surg >= 60) & (survival_data$age_surg < 70)] <- '60-70'
survival_data$agebins[survival_data$age_surg < 60] <- 'under 60'
survival_data$binaryAge <- ifelse(survival_data$age_surg >= 65, "65+", "under 65")
survival_data$os_bin <- ifelse(survival_data$os_sts == "os.event", 1, 0)
survival_data$dss_bin <- ifelse(survival_data$dss_sts == "os.event", 1, 0)

# survival analysis on the p53_abn subset
survival_data_p53 <- survival_data %>% dplyr::filter(status == 'p53_abn')
survival_data_p53 <- survival_data_p53 %>% dplyr::filter(survival_data_p53$sample_id %in% colnames(exposures))
sig_groups <- data.frame(sample_id = colnames(exposures), max_sigs = as.character(apply(exposures, 2, function(x) which.max(x))) )
survival_data_p53 <- survival_data_p53 %>% dplyr::left_join(sig_groups, by = c('sample_id'))
survival_data_p53 <- survival_data_p53 %>% 
                      dplyr::left_join(tibble::rownames_to_column(as.data.frame(t(exposures)), "sample_id"), by = c('sample_id'))

### KM Plots
### Removing the signature classes that don't have enough observations for the km plot to make sense
survival_data_in <- survival_data_p53 %>% dplyr::filter(max_sigs %in% c('1','3','4','5','7'))

km_trt_fit <- survfit(Surv(os_yrs, os_bin) ~ max_sigs, data=survival_data_in)
p <- autoplot(km_trt_fit, conf.int = F, xlab = 'time (years)') +
  labs(title = "Kaplan Meier survival curves for patients grouped by age group") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12))
p
p <- autoplot(km_trt_fit, conf.int = TRUE, xlab = 'time (years)') +
  labs(title = "Kaplan Meier survival curves for patients grouped by age") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12))
p
ggsave(filename = '~/Documents/projects/cn_signatures_shallowWGS/plotting/survival_plots/p53_abn/kapMeier_QDNA_30kb_autosomesOnly_p53_abn_MADabsbritrocSignature_exposures_3.png', 
       plot = p, width = 8, height = 8, units = 'in')

### Tests of Statistical Significance - Accelerated Failure Time Models
# https://www.drizopoulos.com/courses/emc/survival%20analysis%20in%20r%20companion
logrank <- survdiff(Surv(os_yrs, os_bin) ~ max_sigs, data = survival_data_in, subset = (survival_data_in$max_sigs %in% c(1,4)))


# Cox Proportional hazards
cox <- coxph(Surv(os_yrs, os_bin) ~ max_sigs + hist + binaryAge + grade, data = survival_data_in)
summary(cox)
cox_fit <- survfit(cox)
#plot(cox_fit, main = "cph model", xlab="Days")
autoplot(cox_fit)






metadata <- metadata %>% dplyr::filter(str_sub(sample_id,-1,-1) != 'N')
metadata$sample_id[str_sub(metadata$sample_id,-2,-1) == '-T'] <- str_sub(metadata$sample_id,1,11)[str_sub(metadata$sample_id,-2,-1) == '-T']
metadata$sample_id[str_sub(metadata$sample_id,-2,-1) == '-T'] <- str_sub(metadata$sample_id,1,8)[str_sub(metadata$sample_id,-2,-1) == '-T']
metadata$sample_id[str_sub(metadata$sample_id,-1,-1) == 'T'] <- str_sub(metadata$sample_id,1,11)[str_sub(metadata$sample_id,-1,-1) == 'T']
metadata$sample_id[25] <- 'VOA2050'

long_data <- gather(exposures)
long_data$max_sig <- rep(apply(exposures, 2, function(x) which.max(x)), times = 1, each = nsigs)
long_data$sigs <- rep(1:nsigs,dim(exposures)[2])
colnames(long_data) <- c('X', 'Z', 'max_sig', 'Y')
long_data <- long_data %>% arrange(max_sig)
long_data$X <- factor(long_data$X, levels = unique(long_data$X))

