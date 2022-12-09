library(tidyr)
library(ggplot2)
library(magrittr)
library(dplyr)
library(stringr)
library(readr)
library(rascal)

options(dplyr.summarise.inform = FALSE)

variants <- '~/Documents/projects/cn_sigs_swgs/targeted_panel_seq/gatk_and_annotation_output/all.variants.custom_params_GATK.clinvar.cosmic.exons.csv'
variants <- data.table::fread(variants)
xx <- variants %>% dplyr::group_by(samplename) %>% filter(VAF1 == max(VAF1))
table(xx$genecode_gene_name)


#' Calculate Absolute Copy Numbers
#'
#' Use the rascal package in R to do this transformation.
#' Based on instructions in the vignette: \cr
#' https://github.com/crukci-bioinformatics/rascal/blob/master/vignettes/rascal.Rmd \cr \cr
#'
#' CalculateACNs() calculates the absolute copy numbers (ACNs) from the relative copy numbers of one or more samples.
#' There are several included options by which to do this.
#' Note: If not providing a table of of variant allele frequencies (VAFs) then 'mad' is the only method available.
#'
#' @param relative_cns A tsv of the relative copy-numbers.
#' @param rascal_sols A tsv of the calculated rascal solutions.
#' @param variants A dataframe of the variants including variant allele frequencies per gene and per sample.
#' If using variant allele frequencies from targeted panel sequencing or some other technology: \cr
#' - The variants must must be in a datatable/dataframe. \cr
#' - Required columns: sample_id, chromosome, start, end, gene_name, ref, alt, vaf.
#' - Each row of said table must correspond to a unique variant. \cr
#' - Each variant must have an associated variant allele frequency. \cr
#' - Each row must also be associated with a specific sample. \cr
#' @param acnmethod The method by which to calculate ACNs. Can be one of: \cr
#' "maxvaf" - Calculate ACNs assuming the maximum discovered VAF for the sample is an appropriate representation for the tumour fraction. \cr
#' char vector - Same as above but rather than using the max vaf provide a character vector of the genes from which to pull VAFs.
#' The genes are assumed to be in order of decreasing precedence. ex. c('TP53', 'KRAS', 'PTEN') \cr
#' "mad" - Calculate ACNs using the mean absolute difference (MAD) column from the solutions table. \cr
#' @param acn_save_path (optional) The output path (absolute path recommended) where to save the result.
#' @returns A list of dataframes. One DF for each sample where an Absolute Copy Number profile was successfully found.
#' @examples
#' solutions <- "~/Documents/.../rascal_solutions.csv"
#' rcn_segs <- "~/Documents/.../rCN_segs.tsv"
#' variants <- "~/Documents/.../allvariants.clinvar.cosmic.exons.csv"
#' save_path <- "~/Documents/.../rascal_ACN_segments.rds"
#' variants <- data.table::fread(file = variants, sep = ',')
#' variants <- variants %>% dplyr::rename(chromosome = chr,
#'                                        gene_name = genecode_gene_name)
#' variants$sample_id <- stringr::str_replace_all(variants$sample_id, "-", ".")
#' output <- CalculateACNs(relative_cns = rcn_segs,
#'                         rascal_sols = solutions,
#'                         variants = variants,
#'                         acnmethod = 'maxvaf')
#' output <- CalculateACNs(relative_cns = rcn_segs,
#'                         rascal_sols = solutions,
#'                         variants = variants,
#'                         acnmethod = c('TP53', 'KRAS', 'BRCA1',
#'                                       'BRCA2', 'PIK3CA', 'PTEN'),
#'                         acn_save_path = save_path)
#' @export
CalculateACNs <- function (relative_cns, rascal_sols, variants = FALSE, acnmethod, acn_save_path = FALSE) {

  relative_cns <- read.table(file = relative_cns, header = TRUE)
  relative_segments <- relative_cns %>% gather(sample, segmented, 5:dim(.)[2], factor_key=TRUE)    # Convert to long
  segments <- CopyNumberSegments(relative_segments)                            # Collapse to continuous segments

  rascal_batch_solutions <- read.table(file = rascal_sols, sep = ',', header = TRUE)
  rascal_batch_solutions$sample <- str_replace_all(rascal_batch_solutions$sample, "-", ".")

  if (variants != FALSE) {
    if(class(variants)[1] == 'character') {
      variants <- data.table::fread(file = variants, header = TRUE, sep = ",", fill = TRUE)
    } else if ('data.frame' %in% class(variants)) { # do nothing
    } else {
      stop("Invalid value passed to the 'variants' parameter. Please supply either a path or dataframe.")
    }

    variants <- variants %>% dplyr::select(chromosome, start, end, sample_id, gene_name,
                                           starts_with('ref.'),
                                           starts_with('alt.'),
                                           starts_with('vaf')) %>%
      dplyr::mutate(vaf = dplyr::select(., starts_with('vaf')) %>% rowSums(na.rm = TRUE)) %>%
      dplyr::group_by(sample_id, gene_name) %>%
      summarise(sample_id = dplyr::first(sample_id),
                gene_name = dplyr::first(gene_name),
                vaf = max(vaf))
    variants$sample_id <- str_replace_all(variants$sample_id, "-", ".")
  }

  if ((length(acnmethod) == 1) && (acnmethod == 'maxvaf')) {
    variants <- variants %>% group_by(sample_id) %>% dplyr::filter(vaf == max(vaf))
    variants <- variants[!duplicated(variants$sample_id),]
    output <- GenVafAcns(segments, rascal_batch_solutions, variants)
  } else if (length(acnmethod) > 1) {
    variants <- variants %>% group_by(sample_id) %>%
                             dplyr::filter(gene_name %in% acnmethod | vaf == max(vaf, na.rm = TRUE)) %>%
                             dplyr::filter(gene_name == c(intersect(acnmethod, gene_name), setdiff(gene_name, acnmethod))[1])
    output <- GenVafAcns(segments, rascal_batch_solutions, variants)
  } else if (acnmethod == 'mad') {

  } else {
    stop("Invalid value passed to the 'acnmethod' parameter. \n
         Please supply an accepted option.")
  }

  if (acn_save_path != FALSE) {
    saveRDS(output, file = acn_save_path)
  }
  return(output)
}

GenVafAcns <- function (segments, rascal_batch_solutions, variants) {

  chosenSegmentTablesList <- list()
  j <- 1
  for (i in unique(rascal_batch_solutions$sample)) {
    sample_segments <- dplyr::filter(segments, sample == i)
    solutions <- rascal_batch_solutions %>% dplyr::filter(sample == i)
    vafs <- variants %>% dplyr::filter(sample_id == i)
    if (is.na(vafs$sample_id[1])) next                                          # If there isn't a vaf for the sample skip finding an ACN
    rcn_obj <- sample_segments
    suppressWarnings({
      vaf_gene_rcn <- GetSegmentedRcnForGene(rcn_obj, vafs$gene_name)
    })
    if (vaf_gene_rcn == FALSE) { message(paste('No segmented CN call found for',
                                               rcn_obj$sample[1], 'VAF gene.',
                                               sep = ' ')); next}               # If there isn't a CN for the VAF skip finding an ACN
    solution_set <- solutions %>%                                               # Get from running find best fit solution
      dplyr::select(ploidy, cellularity) %>%
      dplyr::mutate(absolute_copy_number = relative_to_absolute_copy_number(vaf_gene_rcn, ploidy, cellularity)) %>%
      dplyr::mutate(tumour_fraction = tumour_fraction(absolute_copy_number, cellularity))

    sol_idx <- which.min(abs(as.double(vafs$vaf) - (solution_set$tumour_fraction*100)))
    solution_set <- solution_set[,]
    absolute_segments <- mutate(sample_segments,
                                copy_number = relative_to_absolute_copy_number(copy_number,
                                                                               solution_set[sol_idx,]$ploidy,
                                                                               solution_set[sol_idx,]$cellularity))
    absolute_segments <- absolute_segments %>% dplyr::select(chromosome=chromosome, start=start, end=end, segVal=copy_number)
    chosenSegmentTablesList[[i]] <- as.data.frame(absolute_segments)
    j <- j+1
  }
  return(chosenSegmentTablesList)
}


####################
# Refactoring - create a new QDNAseq object out of multiple other QDNAseq objects
####################
QD_X <- readRDS(file = '~/Documents/projects/cn_sigs_swgs/copy_number_objects/qdnaseq_Xchr_included/30kb_rCN_comCNVfilt.rds')
WX_X <- readRDS(file = '~/Documents/projects/cn_sigs_swgs/copy_number_objects/wisecondorX_Xchr_included/30kb_rCN_comCNVfilt.rds')

sample_selections <- '~/Documents/projects/cn_sigs_swgs/quality/sample_quality_decision_caller.csv'
sample_selections <- data.table::fread(sample_selections, sep = ',', header = TRUE)
sels <- sample_selections[,c(1,8)]
colnames(sels) <- c('sample_id', 'selection')

BestRelProfile <- function (selection, robjects) {

  if (length(robjects) < 2) {
    stop("Please pass a minimum of two QDNAseq objects to this function.")
  }

  feature_lengths <- sapply(robjects, nrow)
  idx <- which.max(feature_lengths)
  n <- dim(selection)[1]
  dim_names <- rownames(robjects[[idx]]@assayData[["copynumber"]])
  template <- matrix(NA_real_,
                     nrow=max(feature_lengths),
                     ncol=dim(selection)[1],
                     dimnames=list(dim_names, selection$sample_id))
  tbins <- as.data.frame(stringr::str_split(dim_names, pattern = ':|-', simplify = TRUE))
  colnames(bins) <- c('chromosome', 'start', 'end')
  tbins <- tbins %>% dplyr::filter(chromosome != 'Y')
  cns <- tbins
  segs <- tbins
  objstats <- selection
  for (i in names(robjects)) {

    # Pull what needed out of each object in turn to fill the new slots for:
    # copy-number
    mask <- robjects[[i]]@phenoData@data[["name"]] %in% selection$sample_id
    tempmat <- robjects[[i]]@assayData[["copynumber"]][,mask]
    tempmat <- tempmat[,selection$sample_id]
    tempmat <- tempmat[,selection$selection == i]
    tempmat <- cbind(as.data.frame(tempmat),
                     as.data.frame(stringr::str_split(dim_names,
                                                      pattern = ':|-',
                                                      simplify = TRUE)))
    tempmat <- tempmat %>% dplyr::rename(chromosome = V1, start = V2, end = V3)
    cns <- cns %>% dplyr::left_join(tempmat, by = c('chromosome' = 'V1',
                                                    'start' = 'V2',
                                                    'end' = 'V3'))
    # segments
    tempmat <- robjects[[i]]@assayData[["segmented"]][,mask]
    tempmat <- tempmat[,selection$sample_id]
    tempmat <- tempmat[,selection$selection == i]
    tempmat <- cbind(as.data.frame(tempmat),
                     as.data.frame(stringr::str_split(dim_names,
                                                      pattern = ':|-',
                                                      simplify = TRUE)))
    tempmat <- tempmat %>% dplyr::rename(chromosome = V1, start = V2, end = V3)
    segs <- segs %>% dplyr::left_join(tempmat, by = c('chromosome' = 'V1',
                                                      'start' = 'V2',
                                                      'end' = 'V3'))
    # stats
    objstats <- rbind(objstats,
                      robjects[[i]]@phenoData@data[mask, c('total.reads',
                                                           'expected.variance')])

  }

  # Create copy-numbers array
  copynumbers <- matrix(NA_real_, nrow=nrow(cns), ncol=n,
                        dimnames=list(dim_names, selection$sample_id))
  copynumbers[,1:n] <- as.matrix(cns[,4:(dim(cns)[2]+3)])
  # Create segments array
  segments <- matrix(NA_real_, nrow=nrow(segs), ncol=n,
                     dimnames=list(dim_names, selection$sample_id))
  segments[,1:n] <- as.matrix(segs[,4:dim(segs)[2]+3])
  # bins
  bins <- matrix(NA_integer_, nrow=nrow(cns_wide), ncol=4,
                 dimnames=list(dim_names, c('chromosome', 'start', 'end', 'use')))
  bins[,1:3] <- as.matrix(tbins[,1:3])
  bins[,4] <- complete.cases(segments)
  bins <- Biobase::AnnotatedDataFrame(as.data.frame(bins))
  bins@data$chromosome <- factor(bins@data$chromosome, levels = c(as.character(c(1:22)),'X'))
  bins@data$start <- as.integer(bins@data$start)
  bins@data$end <- as.integer(bins@data$end)
  bins@data$use <- as.logical(bins@data$use)



  # Assemble QDNAseq object
  wx_qdnaobj <- new('QDNAseqCopyNumbers', bins=bins, copynumber=copynumbers, phenodata=stats)
  Biobase::assayDataElement(wx_qdnaobj, "segmented") <- segments








  # Create bins annotated dataframe
  bins <- matrix(NA_integer_, nrow=nrow(cns_wide), ncol=4,
                 dimnames=list(dim_names, c('chromosome', 'start', 'end', 'use')))
  bins[,1:3] <- as.matrix(cns_wide[,1:3])
  bins[,4] <- complete.cases(segments)
  bins <- Biobase::AnnotatedDataFrame(as.data.frame(bins))
  bins@data$chromosome <- factor(bins@data$chromosome, levels = c(as.character(c(1:22)),'X'))
  # Assemble QDNAseq object
  wx_qdnaobj <- new('QDNAseqCopyNumbers', bins=bins, copynumber=copynumbers, phenodata=stats)
  assayDataElement(wx_qdnaobj, "segmented") <- segments
  # Change to total.reads
  colnames(wx_qdnaobj@phenoData@data) <- c("name", "total.reads")
  # Calculate expected variance
  expected.variance <- rep(NA_real_, 262)
  for (i in seq_len(262)) {
    expected.variance[i] <- (sum(QDNAseq:::binsToUse(wx_qdnaobj[,i])) / wx_qdnaobj[,i]$total.reads)
  }
  metadata <-  matrix(NA_integer_, nrow= 3, ncol= 1,
                      dimnames=list(c('name', 'total.reads', 'expected.variance')))
  colnames(metadata) <- "labelDescription"
  wx_qdnaobj@phenoData@varMetadata <- as.data.frame(metadata)
  return(wx_qdnaobj)


}

BestRelProfile(sels, list(QDNA = QD_X, WX = WX_X))
robjects <- list(QDNA = QD_X, WX = WX_X)




