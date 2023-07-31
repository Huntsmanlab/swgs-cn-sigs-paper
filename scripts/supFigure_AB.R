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

############################
############################
############################ Plot both signature exposures from Brenton and van together, but order by BrentonSig3, then by BrentonSig7
PlotSpecialSignatureExposures <- function (signatures, save_path = FALSE,
                                           obj_name = 'sig_exposures_obj',
                                           order = FALSE, transpose = FALSE) {
  nsigs <- nrow(signatures)
  long_data <- tidyr::gather(signatures)
  long_data$max_sig <- rep(apply(signatures, 2, function(x) which.max(x)),
                           times = 1,
                           each = nsigs)
  long_data$sigs <- rep(1:nsigs,dim(signatures)[2])
  colnames(long_data) <- c('X', 'Z', 'max_sig', 'Y')
  long_data <- long_data %>% dplyr::arrange(max_sig)
  long_data$X <- factor(long_data$X, levels = unique(long_data$X))

  if (!isFALSE(order)) {
    long_data$X <- factor(long_data$X, levels = order)
  }

  g <- ggplot2::ggplot(long_data, ggplot2::aes(X, Y, fill=Z)) +
    ggplot2::geom_tile() +
    viridis::scale_fill_viridis(discrete=FALSE) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                   axis.ticks.x = ggplot2::element_blank(),
                   # axis.text.x = element_text(size = 15, angle = 75, vjust = 0.5, hjust=0.5),
                   axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size = 18),
                   axis.text.y = ggplot2::element_text(size = 18),
                   legend.title = ggplot2::element_text(size = 18)) +
    ggplot2::labs(fill = 'Signature \nExposure', x = "Samples", y = " ") +
    ggplot2::scale_y_discrete(limits = c(paste0('BS', 1:7), paste0('VS', 1:5))) +
    ggplot2::ggtitle("Signature exposures called per sample")
  g

  if (transpose != FALSE) {
    g <- g + ggplot2::theme(plot.title = ggplot2::element_text(size = 28),
                            axis.ticks.y = ggplot2::element_blank(),
                            axis.text.x = ggplot2::element_text(size = 18),
                            axis.title.x = ggplot2::element_text(size = 18),
                            axis.text.y = ggplot2::element_blank(),
                            legend.title = ggplot2::element_text(size = 18)) +
      ggplot2::coord_flip()
  }
  if (save_path != FALSE) {
    if (transpose != FALSE) {
      ggplot2::ggsave(paste0(save_path, "/signatures_heatmap_", obj_name,".png"), plot = g, width = 10, height = 15)
    } else {
      ggplot2::ggsave(paste0(save_path, "/signatures_heatmap_", obj_name,".png"), plot = g, width = 15, height = 10)
    }
  }

  output <- list(plot = g, ordering = levels(long_data$X))

  return(output)
}

# Arrange samples as intended
temp <- as.data.frame(t(BrSpBrSigs))
colnames(temp) <- paste0('BS', 1:7)
temp <- temp %>% dplyr::arrange(BS7)

temp2 <- rbind(BrSpBrSigs, BrSpVanSigs)

p7 <- PlotSpecialSignatureExposures(temp2, order = rownames(temp),
                                    transpose = TRUE)
p7

png(filename = '~/Desktop/attempt2.png', width=20, height=20, units = 'in', res = 400, type = 'cairo-png')
p7$plot
dev.off()



temp <- as.data.frame(t(VanSpBrSigs))
colnames(temp) <- paste0('BS', 1:7)
temp <- temp %>% dplyr::arrange(BS3)

temp2 <- rbind(VanSpBrSigs, VanSpVanSigs)

p8 <- PlotSpecialSignatureExposures(temp2, order = rownames(temp),
                                    transpose = TRUE)
p8

png(filename = '~/Desktop/attempt1.png', width=20, height=20, units = 'in', res = 400, type = 'cairo-png')
p8$plot
dev.off()


############################
############################
############################ Plot component-by-signature heatmaps for Brenton and Van Juxtaposed
reord_vancouver <- as.integer(c(1,2,5,3,4))
NMF::basismap(signatures, Rowv=NA, Colv=reord_vancouver, main="VancouverP53ABN", tracks=NA)
NMF::basismap(sigs, Rowv=NA, Colv=reord_britroc, main="BritROC", tracks=NA)


############################
############################
############################
############################
############################
############################
############################
############################
############################
############################
############################
############################ Make and layout the figure






feat_sig_mat <- basis(sigs)



############################
############################
############################ Lets try visualizing using CN-Diversity plots first

# First need to pull segmented aCNs out of Brenton QDNAseq object

# Continuing From the britrocSampleProcessing script
brit_CN <- all_CN[,colnames(all_CN) %in% colnames(sig_pat_mat_britroc)]

swide <- as.data.frame(brit_CN@assayData[["segmented"]])
swide$indices <- rownames(swide)
slong <- gather(swide, sample, segmented, IM_100:`JBLAB-4282`)
tbins <- as.data.frame(stringr::str_split(slong$indices, pattern = ':|-', simplify = TRUE))
colnames(tbins) <- c('chromosome', 'start', 'end')
slong <- cbind(tbins, slong)
slong$indices <- NULL
slong <- slong %>% dplyr::filter(!(chromosome %in% c('Y', 'X')))                # Brenton group doesn't have X chromosome data
slong$start <- as.integer(slong$start)
slong$end <- as.integer(slong$end)
cns <- CopyNumberSegments(slong)                                                # Collapse to segments
segs <- cns[,c(2:5,8)]

# Now convert segmented aCNs to segs list and then convert back to long
slist <- list()
for (i in unique(segs$sample)) {
  slist[[i]] <- segs %>% dplyr::filter(sample == i) %>%
    dplyr::select(chromosome, start, end, copy_number) %>%
    dplyr::rename(segVal = copy_number)
}

segs2 <- SegmentsToCopyNumber(slist, 1000000,
                              genome = 'hg19', Xincluded = FALSE)
pp_britov <- SwgsCnvHeatmaps(reads = segs2, ret_order = T)

p <- plot_grid(pp_britov$plot, pp_vanendo$plot, nrow = 1)
png(filename = '~/Desktop/attempt2.png', width=40, height=20, units = 'in', res = 400, type = 'cairo-png')
p
dev.off()

#####
#####
##### CN-Diversity + Signature Exposures Plot

p2 <- SwgsCnvHeatmaps(reads = segs2, order = p1$ordering)
p3 <- plot_grid(p1$plot, p2, nrow = 1)
png(filename = '~/Desktop/attempt5.png', width=40, height=20, units = 'in', res = 400, type = 'cairo-png')
p3
dev.off()
