suppressPackageStartupMessages({
  library(tidyr)
  library(ggplot2)
  library(magrittr)
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

# Sup_Fig_A
# Read in data
segments <- readRDS(file = '~/Documents/projects/cn_sigs_swgs/copy_number_objects/combined_30kb_Xchr_included/30kb_aCN_comCNVfilt_187_filter_false.rds')
components <- '~/repos/utanosmodellingdata/component_models/30kb_ovarian/component_models_britroc_aCNs.rds'
sigsB_path <- '~/repos/utanosmodellingdata/signatures/30kb_ovarian/component_by_signature_britroc_aCNs.rds'
segments_noX <- segments
for (i in 1:length(segments)) {
  segments_noX[[i]] <- segments_noX[[i]] %>% dplyr::filter(!(chromosome %in% c('X', 'Y')))
}

copy_numbers_noX <- SegmentsToCopyNumber(segments_noX, 1000000,
                                     genome = 'hg19', Xincluded = FALSE)
VanSpBrSigs <- CallSignatureExposures(segments_noX,
                                      component_models = components,
                                      signatures = sigsB_path,
                                      refgenome = 'hg19')
reord_britroc <- as.integer(c(2,6,5,4,7,3,1))
VanSpBrSigs <- VanSpBrSigs[reord_britroc,]

plot1 <- PlotSignatureExposures(VanSpBrSigs, transpose = TRUE)


pp <- SwgsCnvHeatmaps(reads = copy_numbers_noX, order = plot1$ordering, ret_order = T)
pp$plot

library(cowplot)


p1 <- plot_grid(pp$plot, plot1$plot, nrow = 1, rel_widths = c(3,2) )

png(filename = '~/Desktop/attempt1.png', width=40, height=20, units = 'in', res = 400, type = 'cairo-png')
p1
dev.off()


# Sup_Fig_B
sigsV <- readRDS(file = '~/repos/utanosmodellingdata/signatures/30kb_endometrial/component_by_signature_britmodelsvansigs5_aCNs.rds')
sig_pat_mat <- NMF::scoef(sigsV)
sig_pat_mat_hq <- NormaliseMatrix(sig_pat_mat)

copy_numbers <- SegmentsToCopyNumber(segments, 1000000,
                                     genome = 'hg19', Xincluded = FALSE)

segs_lq <- segments[!(names(segments) %in% colnames(sig_pat_mat_hq))]
sig_pat_mat_lq <- CallSignatureExposures(segs_lq,
                                         component_models = components,
                                         signatures = sigsV_path,
                                         refgenome = 'hg19')
VanSpVanSigs <- cbind(sig_pat_mat_hq, sig_pat_mat_lq)

plot2 <- PlotSignatureExposures(VanSpVanSigs, transpose = TRUE)
pp2 <- SwgsCnvHeatmaps(reads = copy_numbers, order = plot2$ordering, ret_order = T)
p2 <- plot_grid(pp2$plot, plot2$plot, nrow = 1, rel_widths = c(3,2) )

png(filename = '~/Desktop/attempt2.png', width=40, height=20, units = 'in', res = 400, type = 'cairo-png')
p2
dev.off()
