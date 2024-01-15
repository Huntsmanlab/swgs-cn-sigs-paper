# Code to generate an upset plot, showing the intersection of all possible combinations of mutations present in the data

# Required packages
library(dplyr)
library(ComplexUpset)
library(readxl)

# Read in the data and generate the combination matrix
mut_table <- read_xlsx(path = "./mutation_calls/mutations_consolidated_recode_vs5.xlsx")
sig_exposures <- read.csv(file = "./sig_exposures/agglomerated_exposures_table.csv", row.names = 1) %>%
  mutate(sample_id = rownames(.))

mut_table_recode <- mut_table %>%
  mutate_at(vars(!(matches("sample_id"))), ~ case_match(., c("Negative", "Absent") ~ "FALSE", c("Positive", "Present") ~ "TRUE"), "No Data" ~ NA) %>%
  mutate_at(vars(!(matches("sample_id"))), as.logical) %>%
  left_join(y = dplyr::select(.data = sig_exposures, sample_id, max_van_sig, max_brenton_sig), by = join_by(sample_id)) %>%
  data.frame()

rownames(mut_table_recode) <- mut_table_recode[, 1]
names(mut_table_recode)[5] <- "FBXW7/PPP2R1A"
mut_table_recode <- mut_table_recode[, -1]
mut_table_recode <-mut_table_recode[complete.cases(mut_table_recode), ] # Only keep cases with mutation information for all categories



# Generate the plot, with an exclusive intersection. I.e, in sets A and B but not in C, or in A but not in any other set
upset(mut_table_recode, names(mut_table_recode)[1:4], mode = "distinct", name = "Group combination",
      set_sizes = (upset_set_size() + theme(text = element_text(size = 14))),
      base_annotations = list("Intersection size" = (intersection_size() + theme(text = element_text(size = 14))))) +
  theme(text = element_text(size = 14))
upset(mut_table_hrd_recode, names(mut_table_hrd_recode)[-1], mode = "distinct", name = "Group combination")
