suppressPackageStartupMessages({
  library(tidyr)
  library(ggplot2)
  library(magrittr)
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(easyalluvial)
})

######
######
###### Get/calculate Signature exposures for each cohort and each set of signatures
###### 4x comparisons:
###### Britroc HGSOC samples -> Britroc HGSOC Sigs
###### Britroc HGSOC samples -> Vancouver p53abn Endometrial Sigs
###### Vancouver p53abn Endometrial samples -> Britroc HGSOC Sigs
###### Vancouver p53abn Endometrial samples -> Vancouver p53abn Endometrial Sigs

# Part 1
# ------
# Grab Britroc HGSOC signatures so we can get the exposures for the 91 samples used for modelling/creation
sigsB <- readRDS(file = './utanosmodellingdata/signatures/30kb_ovarian/component_by_signature_britroc_aCNs.rds')
sig_pat_mat_hq <- NMF::scoef(sigsB)
reord_britroc <- as.integer(c(2,6,5,4,7,3,1))
sig_pat_mat_hq <- sig_pat_mat_hq[reord_britroc,]
sig_pat_mat_hq <- NormaliseMatrix(sig_pat_mat_hq)

# Determine exposures to the lower-quality (2-star) samples
components <- './utanosmodellingdata/component_models/30kb_ovarian/component_models_britroc_aCNs.rds'
sigsB_path <- './utanosmodellingdata/signatures/30kb_ovarian/component_by_signature_britroc_aCNs.rds'

# Get annotations for Britroc samples From Macintyre et al. 2018
# Simplified/copy-pasted code in britrocSampleProcessing.R script
# Run said code so that we have the 'result' dataframe in memory
ids <- result %>% filter(star_rating==2)
ids <- ids$IM.JBLAB_ID
lq_CN <- all_CN[, colnames(all_CN) %in% ids]

sig_pat_mat_lq <- CallSignatureExposures(lq_CN,
                                         component_models = components,
                                         signatures = sigsB_path,
                                         refgenome = 'hg19')
sig_pat_mat_lq <- sig_pat_mat_lq[reord_britroc,]

BrSpBrSigs <- cbind(sig_pat_mat_hq, sig_pat_mat_lq)

# Part 2
# ------
# Determine exposures for all Britroc samples together to p53abn endo signatures
sigsV_path <- './utanosmodellingdata/signatures/30kb_endometrial/component_by_signature_britmodelsvansigs5_aCNs.rds'
ids <- result %>% filter(star_rating %in% c(2, 3))
ids <- ids$IM.JBLAB_ID
hqlq_CN <- all_CN[, colnames(all_CN) %in% ids]

sig_pat_mat <- CallSignatureExposures(hqlq_CN,
                                         component_models = components,
                                         signatures = sigsV_path,
                                         refgenome = 'hg19')
BrSpVanSigs <- sig_pat_mat


######
###### Get/calculate Vancouver p53abn Endometrial samples signature exposures
###### for both p53abn endometrial and HGSOC signatures

# Part 3
# ------

sigsV <- readRDS(file = './utanosmodellingdata/signatures/30kb_endometrial/component_by_signature_britmodelsvansigs5_aCNs.rds')
sig_pat_mat <- NMF::scoef(sigsV)
sig_pat_mat_hq <- NormaliseMatrix(sig_pat_mat)

segs <- readRDS(file = './cn_objects/combined_30kb_Xchr_included/30kb_aCN_comCNVfilt_187_filter_false.rds')
segs_lq <- segs[!(names(segs) %in% colnames(sig_pat_mat_hq))]
sig_pat_mat_lq <- CallSignatureExposures(segs_lq,
                                         component_models = components,
                                         signatures = sigsV_path,
                                         refgenome = 'hg19')
VanSpVanSigs <- cbind(sig_pat_mat_hq, sig_pat_mat_lq)

# Part 4
# ------
for (i in 1:length(segs)) {
  segs[[i]] <- segs[[i]] %>% dplyr::filter(!(chromosome %in% c('X', 'Y')))
}
sig_pat_mat <- CallSignatureExposures(segs,
                                      component_models = components,
                                      signatures = sigsB_path,
                                      refgenome = 'hg19')
VanSpBrSigs <- sig_pat_mat[reord_britroc,]

# Part 5
# ------
# Combine signature exposures into a single dataframe for convenience/ease-of-use

VanExposures1 <- as.data.frame(t(VanSpVanSigs))
VanExposures2 <- as.data.frame(t(BrSpVanSigs))
VanExposures <- rbind(VanExposures1, VanExposures2)
colnames(VanExposures) <- paste0('VS', 1:dim(VanExposures)[2])
VanExposures$max_van_sig <- apply(VanExposures, 1, function(x) which.max(x))
VanExposures$max_van_sig <- paste0('VS', VanExposures$max_van_sig)
VanExposures$max_van_sig <- factor(VanExposures$max_van_sig, levels = c('VS2', 'VS4',
                                                                        'VS3', 'VS1',
                                                                        'VS5'))
VanExposures$sample_ids <- rownames(VanExposures)

BrExposures1 <- as.data.frame(t(VanSpBrSigs))
BrExposures2 <- as.data.frame(t(BrSpBrSigs))
BrExposures <- rbind(BrExposures1, BrExposures2)
colnames(BrExposures) <- paste0('BS', 1:dim(BrExposures)[2])
BrExposures$max_brenton_sig <- apply(BrExposures, 1, function(x) which.max(x))
BrExposures$max_brenton_sig <- paste0('BS', BrExposures$max_brenton_sig)
BrExposures$sample_ids <- rownames(BrExposures)

plotting_data <- VanExposures %>% dplyr::left_join(BrExposures, by = c('sample_ids'))
rownames(plotting_data) <- plotting_data$sample_ids
plotting_data$sample_ids <- NULL

write.csv(plotting_data, file = './sig_exposures/agglomerated_exposures_table.csv')



plotting_data$cohort <- c(rep('VanP53ABNendo',187), rep('BritrocHGSOC',117))



p9 <- alluvial_wide( select(plotting_data[1:187,], max_brenton_sig, max_van_sig),
                     fill_by = 'first_variable',
                     stratum_label_size = 3.5) +
  ggplot2::labs(y = "Vancouver Cohort Samples", caption = '') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        title = element_blank())
p10 <- alluvial_wide( select(plotting_data[188:304,], max_brenton_sig, max_van_sig),
                      fill_by = 'first_variable',
                      stratum_label_size = 3.5) +
  ggplot2::labs(y = "Britroc Cohort Samples", caption = '') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        title = element_blank())

# Bar plots Britroc Signatures
sigs <- readRDS(file = './utanosmodellingdata/signatures/30kb_ovarian/component_by_signature_britroc_aCNs.rds')
reord_britroc <- as.integer(c(2,6,5,4,7,3,1))
feat_sig_mat <- NMF::basis(sigs)
feat_sig_mat <- feat_sig_mat[,reord_britroc]
colnames(feat_sig_mat) <- paste0("s",1:7)
sig_feat_mat <- t(feat_sig_mat)
temp <- as.data.frame(sig_feat_mat)
norm_const <- apply(temp,2,sum)
temp <- data.frame(t(apply(temp,1,function(x){x/norm_const})))
temp$sig <- rownames(sig_feat_mat)
pdat <- reshape2::melt(temp,id.vars="sig")
vars <- gsub("^\\d+|\\d+$", "",as.character(pdat$variable))
pdat <- cbind(pdat,vars)
colnames(pdat) <- c("sig","Feature","value","Distribution")
pdat$sig <- factor(pdat$sig,levels=paste0("s",1:7))
pdat$Distribution <- plyr::revalue(pdat$Distribution,
                                 c(bp10MB="Breakpoint number",
                                   copynumber="Copy-number",
                                   bpchrarm="Breakpoints per chr arm",
                                   changepoint="CN changepoint",
                                   segsize="Segment size",
                                   osCN="Oscilating CN length"))
pdat$Distribution <- factor(pdat$Distribution,levels=c("Breakpoint number","Copy-number","CN changepoint","Breakpoints per chr arm","Oscilating CN length","Segment size"))
pdat$sig <- plyr::revalue(pdat$sig,
                        c(s1=1, s2=2, s3=3, s4=4, s5=5, s6=6, s7=7))

bp1 <- ggplot(pdat[pdat$sig == 1,], aes(x = interaction(Feature, Distribution),
                               y = value, fill = Distribution,
                               group = Distribution))+
  geom_col(position="dodge") +
  scale_x_discrete(labels=c(1:3,1:8,1:7,1:5,1:3,1:10)) +
  theme(legend.position="none", axis.text = element_text(size = 6),
        axis.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0,1)) +
  ggplot2::labs(title = "B. Signature 1") +
  scale_fill_viridis_d(option="turbo")
bp2 <- ggplot(pdat[pdat$sig == 2,], aes(x = interaction(Feature, Distribution),
                                      y = value, fill = Distribution,
                                      group = Distribution))+
  geom_col(position="dodge") +
  scale_x_discrete(labels=c(1:3,1:8,1:7,1:5,1:3,1:10)) +
  theme(legend.position="none", axis.text = element_text(size = 6),
        axis.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0,1)) +
  ggplot2::labs(title = "B. Signature 2") +
  scale_fill_viridis_d(option="turbo")
bp3 <- ggplot(pdat[pdat$sig == 3,], aes(x = interaction(Feature, Distribution),
                                        y = value, fill = Distribution,
                                        group = Distribution))+
  geom_col(position="dodge") +
  scale_x_discrete(labels=c(1:3,1:8,1:7,1:5,1:3,1:10)) +
  theme(legend.position="none", axis.text = element_text(size = 6),
        axis.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0,1)) +
  ggplot2::labs(title = "B. Signature 3") +
  scale_fill_viridis_d(option="turbo")
bp4 <- ggplot(pdat[pdat$sig == 4,], aes(x = interaction(Feature, Distribution),
                                        y = value, fill = Distribution,
                                        group = Distribution))+
  geom_col(position="dodge") +
  scale_x_discrete(labels=c(1:3,1:8,1:7,1:5,1:3,1:10)) +
  theme(legend.position="none", axis.text = element_text(size = 6),
        axis.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0,1)) +
  ggplot2::labs(title = "B. Signature 4") +
  scale_fill_viridis_d(option="turbo")
bp5 <- ggplot(pdat[pdat$sig == 5,], aes(x = interaction(Feature, Distribution),
                                        y = value, fill = Distribution,
                                        group = Distribution))+
  geom_col(position="dodge") +
  scale_x_discrete(labels=c(1:3,1:8,1:7,1:5,1:3,1:10)) +
  theme(legend.position="none", axis.text = element_text(size = 6),
        axis.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0,1)) +
  ggplot2::labs(title = "B. Signature 5") +
  scale_fill_viridis_d(option="turbo")
bp6 <- ggplot(pdat[pdat$sig == 6,], aes(x = interaction(Feature, Distribution),
                                        y = value, fill = Distribution,
                                        group = Distribution))+
  geom_col(position="dodge") +
  scale_x_discrete(labels=c(1:3,1:8,1:7,1:5,1:3,1:10)) +
  theme(legend.position="none", axis.text = element_text(size = 6),
        axis.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0,1)) +
  ggplot2::labs(title = "Signature 6") +
  scale_fill_viridis_d(option="turbo")
bp7 <- ggplot(pdat[pdat$sig == 7,], aes(x = interaction(Feature, Distribution),
                                        y = value, fill = Distribution,
                                        group = Distribution))+
  geom_col(position="dodge") +
  scale_x_discrete(labels=c(1:3,1:8,1:7,1:5,1:3,1:10)) +
  theme(legend.position="none", axis.text = element_text(size = 6),
        axis.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0,1)) +
  ggplot2::labs(title = "B. Signature 7") +
  scale_fill_viridis_d(option="turbo")


# Bar plots Vancouver Signatures
signatures <- readRDS(file = './utanosmodellingdata/signatures/30kb_endometrial/component_by_signature_britmodelsvansigs5_aCNs.rds')
feat_sig_mat <- NMF::basis(signatures)
colnames(feat_sig_mat) <- paste0("s",1:5)
sig_feat_mat <- t(feat_sig_mat)
temp <- as.data.frame(sig_feat_mat)
norm_const <- apply(temp,2,sum)
temp <- data.frame(t(apply(temp,1,function(x){x/norm_const})))
temp$sig <- rownames(sig_feat_mat)
pdat <- reshape2::melt(temp,id.vars="sig")
vars <- gsub("^\\d+|\\d+$", "",as.character(pdat$variable))
pdat <- cbind(pdat,vars)
colnames(pdat) <- c("sig","Feature","value","Distribution")
pdat$sig <- factor(pdat$sig,levels=paste0("s",1:5))
pdat$Distribution <- plyr::revalue(pdat$Distribution,
                                   c(bp10MB="Breakpoint number",
                                     copynumber="Copy-number",
                                     bpchrarm="Breakpoints per chr arm",
                                     changepoint="CN changepoint",
                                     segsize="Segment size",
                                     osCN="Oscilating CN length"))
pdat$Distribution <- factor(pdat$Distribution,levels=c("Breakpoint number","Copy-number","CN changepoint","Breakpoints per chr arm","Oscilating CN length","Segment size"))
pdat$sig <- plyr::revalue(pdat$sig,
                          c(s1=1, s2=2, s3=3, s4=4, s5=5))

bp8 <- ggplot(pdat[pdat$sig == 1,], aes(x = interaction(Feature, Distribution),
                                        y = value, fill = Distribution,
                                        group = Distribution))+
  geom_col(position="dodge") +
  scale_x_discrete(labels=c(1:3,1:8,1:7,1:5,1:3,1:10)) +
  theme(legend.position="none", axis.text = element_text(size = 6),
        axis.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0,1)) +
  ggplot2::labs(title = "V. Signature 1") +
  scale_fill_viridis_d(option="turbo")
bp9 <- ggplot(pdat[pdat$sig == 2,], aes(x = interaction(Feature, Distribution),
                                        y = value, fill = Distribution,
                                        group = Distribution))+
  geom_col(position="dodge") +
  scale_x_discrete(labels=c(1:3,1:8,1:7,1:5,1:3,1:10)) +
  theme(legend.position="none", axis.text = element_text(size = 6),
        axis.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0,1)) +
  ggplot2::labs(title = "V. Signature 2") +
  scale_fill_viridis_d(option="turbo")
bp10 <- ggplot(pdat[pdat$sig == 3,], aes(x = interaction(Feature, Distribution),
                                        y = value, fill = Distribution,
                                        group = Distribution))+
  geom_col(position="dodge") +
  scale_x_discrete(labels=c(1:3,1:8,1:7,1:5,1:3,1:10)) +
  theme(legend.position="none", axis.text = element_text(size = 6),
        axis.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0,1)) +
  ggplot2::labs(title = "V. Signature 3") +
  scale_fill_viridis_d(option="turbo")
bp11 <- ggplot(pdat[pdat$sig == 4,], aes(x = interaction(Feature, Distribution),
                                        y = value, fill = Distribution,
                                        group = Distribution))+
  geom_col(position="dodge") +
  scale_x_discrete(labels=c(1:3,1:8,1:7,1:5,1:3,1:10)) +
  theme(legend.position="none", axis.text = element_text(size = 6),
        axis.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0,1)) +
  ggplot2::labs(title = "V. Signature 4") +
  scale_fill_viridis_d(option="turbo")
bp12 <- ggplot(pdat[pdat$sig == 5,], aes(x = interaction(Feature, Distribution),
                                        y = value, fill = Distribution,
                                        group = Distribution))+
  geom_col(position="dodge") +
  scale_x_discrete(labels=c(1:3,1:8,1:7,1:5,1:3,1:10)) +
  theme(legend.position="bottom", axis.text = element_text(size = 6),
        axis.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0,1)) +
  ggplot2::labs(title = "V. Signature 5") +
  scale_fill_viridis_d(option="turbo")

library(patchwork)

A <- (bp7 / bp6 / bp5 / bp4 / bp3 / bp2 / bp1)
B <- ((p9 / p10) | (bp12 / bp11 / bp10 / bp9 / bp8 )) + plot_layout(widths = c(3, 2))
B <- B / guide_area() + plot_layout(guides = 'collect', heights = c(20, 1))
alluvial_fig <- (A | B) + plot_layout(widths = c(2, 5))

png(filename = '~/Desktop/attempt3.png', width=9, height=11, units = 'in', res = 400, type = 'cairo-png')
alluvial_fig
dev.off()






