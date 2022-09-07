library(dplyr)
library(stringr)
library(tidyr)
library(magrittr)
library(data.table)
library(Biobase)
library(QDNAseq)

swgs_dirs <- list.dirs('~/Downloads/npzs_redo')
swgs_dirs <- swgs_dirs[2:length(swgs_dirs)]
samplenames <- sub(swgs_dirs, pattern = '/Users/mdouglas/Downloads/npzs_redo/', replacement = '')
samplenames <- sub(samplenames, pattern = '.plots', replacement = '')

for (i in swgs_dirs) {
  samplename <- sub(i, pattern = '/Users/mdouglas/Downloads/npzs_redo/', replacement = '')
  samplename <- sub(samplename, pattern = '.plots', replacement = '')
  system(command = paste('mv', paste0(i, '/genome_wide.png'), paste0(i, '/', samplename, '_plot.png')))
  print(paste('mv', paste0(i, '/genome_wide.png'), paste0(i, '/', samplename, '_plot.png')))
} 

qdnaobj <- `30kb_rCN_comCNVfilt`
plot(qdnaobj[,1])



s1_s <- data.table::fread(input = '~/Documents/projects/cn_sigs_swgs/copy_number_objects/wisecondorX_Xchr_included/CC-CHM-1341_segments.bed')
s1_b <- data.table::fread(input = '~/Documents/projects/cn_sigs_swgs/copy_number_objects/wisecondorX_Xchr_included/CC-CHM-1341_bins.bed')
s2_s <- data.table::fread(input = '~/Documents/projects/cn_sigs_swgs/copy_number_objects/wisecondorX_Xchr_included/CC-CHM-1347_segments.bed')
s2_b <- data.table::fread(input = '~/Documents/projects/cn_sigs_swgs/copy_number_objects/wisecondorX_Xchr_included/CC-CHM-1347_bins.bed')

colnames(s1_b) <- c('chromosome', 'start', 'end', 'id', 'ratio', 'zscore')
bins <- matrix(NA_integer_, nrow=nrow(s1_b), ncol=4,
                 dimnames=list(s1_b$id, c('chromosome', 'start', 'end', 'zscore')))
bins[,1:4] <- as.matrix(s1_b[,c(1,2,3,6)])
bins <- AnnotatedDataFrame(as.data.frame(bins))
bins@data$start <- as.integer(bins@data$start)
bins@data$end <- as.integer(bins@data$end)
bins@data$zscore <- as.numeric(bins@data$zscore)

counts <- matrix(NA_integer_, nrow=nrow(s1_b), ncol=2,
                 dimnames=list(s1_b$id, c('CC-CHM-1341', 'CC-CHM-1347')))
counts[, 1] <- s1_b$ratio
counts[, 2] <- s2_b$ratio

segments <- matrix(NA_integer_, nrow=nrow(s1_b), ncol=2,
                   dimnames=list(s1_b$id, c('CC-CHM-1341', 'CC-CHM-1347')))
#
# Expand segments using cn-sigs-utils function
#

phenodata <- data.frame(name=c('CC-CHM-1341', 'CC-CHM-1347'), row.names=c('CC-CHM-1341', 'CC-CHM-1347'),
                        stringsAsFactors=FALSE)
phenodata$total.reads <- colSums(counts, na.rm = T)
# phenodata$used.reads <- colSums2(counts, rows=condition, useNames=FALSE)

object <- new('QDNAseqCopyNumbers', bins=bins, copynumber=counts,
              phenodata=phenodata)
assayDataElement(object, "segmented") <- copynumber2

plot(object, logTransform=FALSE, ylim=c(-3, 5))
# object$expected.variance <- expectedVariance(object)

# Plot wisecondorX CN profiles using the same styling and technique as QDNAseq to ease comparability
# Parameters:
# input_path - Path to the wisecondorX output. Expects a directory.

plotWisecondorProfiles <- function (input_path) {
  
  # Read in files
  bin_files <- list.files(input_path, pattern = '*_bins.bed', full.names = TRUE)                            # grab just 1
  samples <- list.files(input_path, pattern = '*_bins.bed')
  samples <- sub('_bins.bed', '', samples)
  seg_files <- list.files(input_path, pattern = '*_segments.bed', full.names = TRUE)
  stats_files <- list.files(input_path, pattern = '*_statistics.txt', full.names = TRUE)
  cns <- plyr::ldply(seq_along(bin_files), function(x) {
    df <- data.table::fread(bin_files[[x]], sep = '\t', col.names = c('chromosome', 'start', 'end', 'id', 'ratio', 'zscore'), nThread = 4)
    df <- df[,c(1,2,3,5)]
    df$sample_id <- samples[x]
    colnames(df) = c('chromosome', 'start', 'end', 'state', 'sample_id')
    df})
  segs <- lapply(seg_files, function(x) {
    df <- data.table::fread(x, sep = '\t', col.names = c('chromosome', 'start', 'end', 'segVal', 'zscore'))
    df <- df[,c(1,2,3,4)]
    df})
  names(segs) <- samples
  stats <- plyr::ldply(seq_along(stats_files), function(x) {
    # df <- data.table::fread(stats_files[[x]], sep = '\t', )
    read_count <- as.integer(sub("Number of reads: ", "", readLines(stats_files[[x]])[26]))
    df <- data.frame(name = samples[x], 
                     total_reads = read_count, 
                     stringsAsFactors=FALSE)
    df})
  rownames(stats) <- samples
  

  # Create the copynumber object
  cns_wide <- spread(cns, sample_id, state)
  dim_names <- paste0(cns_wide$chromosome, ':', cns_wide$start, '-', cns_wide$end)       # create bin dimnames
  copynumbers <- matrix(NA_real_, nrow=nrow(cns_wide), ncol=length(bin_files),
                   dimnames=list(dim_names, samples))
  copynumbers[,1:length(bin_files)] <- as.matrix(cns_wide[,4:(length(bin_files)+3)])
  # Create the segment object
  segs_long <- segments_to_copy_number(segs, bin_size = 30000, Xincluded = TRUE)
  segs_wide <- spread(segs_long, sample_id, state)
  segments <- matrix(NA_real_, nrow=nrow(segs_wide), ncol=length(seg_files),
                        dimnames=list(dim_names, samples))
  segments[,1:length(seg_files)] <- as.matrix(segs_wide[,4:(length(seg_files)+3)])
  # Create the bins
  bins <- matrix(NA_integer_, nrow=nrow(cns_wide), ncol=4,
                 dimnames=list(dim_names, c('chromosome', 'start', 'end', 'use')))
  bins[,1:3] <- as.matrix(cns_wide[,1:3])
  bins[,4] <- complete.cases(segments)
  bins <- AnnotatedDataFrame(as.data.frame(bins))
  bins@data$chromosome <- factor(bins@data$chromosome, levels = c(as.character(c(1:22)),'X'))
  bins@data$start <- as.integer(bins@data$start)
  bins@data$end <- as.integer(bins@data$end)
  bins@data$use <- as.logical(bins@data$use)

  # Assemble QDNAseq object
  wx_qdnaobj <- new('QDNAseqCopyNumbers', bins=bins, copynumber=copynumbers,
                phenodata=stats)
  assayDataElement(wx_qdnaobj, "segmented") <- segments
  
  browser()
  plot(wx_qdnaobj[,1], logTransform=FALSE, ylim=c(-3, 5))
  plot(wx_qdnaobj[,2], logTransform=FALSE, ylim=c(-3, 5))
  
}

.makeSegmentsMod <- function(data,chrdata) {
  previous    <- 2000
  chrpr       <- -100
  values      <- c()
  start       <- c()
  end         <- c()
  el          <- length(data)
  data <- c(data,-10000) #add value to allow data[i+1]
  for (i in 1:el) {
    browser()
    if ((data[i] != previous & previous != data[i+1]) | chrdata[i] != chrpr) { #bug repaired 12/06/09
      start   <- c(start, i)
      last    <- i - 1
      if (last > 0) end <- c(end, last)
      values  <- c(values, data[i])
    }
    previous    <- data[i]
    chrpr <- chrdata[i]
  }
  end     <- c(end, el)
  result  <- cbind(values, start, end)
  result
}

nsmpgroup <- c('EC053', 'EC070', 'EC083', 'VOA658', 'VOA707', 'VOA844', 'EC019', 
               'EC059', 'EC146', 'VOA1111', 'VS03-21607', 'VS05-32611', 'VS06-19921')
input_path <- '~/Documents/projects/cn_sigs_swgs/wisecondorX_QDNA_output_comparison/raw_wisecondorX_output/'
plotWisecondorProfiles(input_path)
