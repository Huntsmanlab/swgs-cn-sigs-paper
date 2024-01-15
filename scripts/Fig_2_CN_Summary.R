# Script to generate various CN summary information.

# Required packages.
library(QDNAseq)
library(utanos)

# Read in the CN data, as stored in a QDNAseq object.
rcn_data <- readRDS(file = "./cn_objects/30kb_rCN_comCNVfilt_196_cutoff.rds") # CNs called with cutoffs: deletions < -0.1, amplifications > 0.1

# Generate the summary plot; we don't mask the plot because we have already selected a CN cutoff and want to visualize everything that remains.
SummaryCNPlot(rcn_data, main = "Relative Copy-number Frequency Plot", maskprob = 0, maskaberr = 0)
