# Script to assess if the distribution of CN signatures differ between serous + carcinosarcoma and non-serous carcinoma

# 1) Does the distribution of CN signatures differ between serous + carcinosarcoma and non-serous carcinoma

# Required packages
library(dplyr)
library(tidyr)

# Read in the data
van_exp <- read.delim(file = "~/Downloads/agglomerated_exposures_table.csv", sep = ",", header = TRUE)
histotypes <- read.delim(file = "~/Downloads/Data/cancer_type_rescored_histotypes_187.csv", sep = ",", header = TRUE)

# Clean up the data and join
van_exp_clean <- van_exp %>%
  rename(sample = X) %>%
  mutate_at(.vars = c("sample"), .funs = gsub, pattern = "\\.", replacement = "-")

van_exp_hist <- inner_join(x = van_exp_clean, y = histotypes, by = join_by(sample == sample_id))

# Generate contingency tables and run tests
van_exp_hist_summary <- van_exp_hist %>%
  mutate(cancer_type = if_else(condition = cancer_type == "serous" | cancer_type == "carcinosarcoma", true = "serous/carcinosarcoma", false = "non-serous")) %>%
  summarise(.by = c("cancer_type", "max_van_sig"), count = n()) %>%
  arrange(cancer_type, max_van_sig)

van_exp_hist_matrix <- matrix(data = van_exp_hist_summary$count, nrow = 2, ncol = 5, dimnames = list(c("Non-serous", "Serous/carcinosarcoma"), paste0("VS", 1:5)), byrow = TRUE)
van_exp_hist_chi <- chisq.test(x = van_exp_hist_matrix)
van_exp_hist_fisher <- fisher.test(x = van_exp_hist_matrix)

vs1_vs2 <- fisher.test(x = van_exp_hist_matrix[, 1:2], conf.level = 0.99)
vs2_vs3 <- fisher.test(x = van_exp_hist_matrix[, 2:3], conf.level = 0.99)
vs2_vs4 <- fisher.test(x = van_exp_hist_matrix[, c(2, 4)], conf.level = 0.99)
vs2_vs5 <- fisher.test(x = van_exp_hist_matrix[, c(2, 5)], conf.level = 0.99)

# Same thing with logistic regression
van_exp_hist_recode <- van_exp_hist_summary %>%
  pivot_wider(names_from = cancer_type, values_from = count) %>%
  rename("non_serous" = "non-serous", "serous" = "serous/carcinosarcoma")

mod <- glm(cbind(non_serous, serous) ~ max_van_sig, data = van_exp_hist_recode, family = binomial)
binom_log_results <- emmeans(mod, pairwise ~ max_van_sig, type = "response")
