library(dplyr)
library(stringr)
library(tidyr)
library(magrittr)
library(data.table)
library(Biobase)
library(QDNAseq)


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
  # NaNs to NAs
  cns$state[is.nan(cns$state)] <- NA
  # Exponeniate 
  cns$state <- exp(cns$state)
  segs <- lapply(seg_files, function(x) {
    df <- data.table::fread(x, sep = '\t', col.names = c('chromosome', 'start', 'end', 'segVal', 'zscore'))
    df <- df[,c(1,2,3,4)]
    df})
  names(segs) <- samples
  segs <- segs %>% 
    purrr::map(~mutate_at(.x, vars("segVal"), function(x) {exp(x)}))
  stats <- plyr::ldply(seq_along(stats_files), function(x) {
    df <- data.table::fread(stats_files[[x]], sep = '\t', )
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
  segs_long <- segments_to_copy_number(segs, bin_size = 30000, genome = 'Xchr_included')
  segs_wide <- spread(segs_long, sample_id, state)
  browser()
  segments <- matrix(NA_real_, nrow=nrow(segs_wide), ncol=length(seg_files),
                        dimnames=list(dim_names, samples))
  segments[,1:length(seg_files)] <- as.matrix(segs_wide[,4:(length(seg_files)+3)])
  # Create the bins
  bins <- matrix(NA_integer_, nrow=nrow(cns_wide), ncol=4,
                 dimnames=list(dim_names, c('chromosome', 'start', 'end', 'use')))
  bins[,1:3] <- as.matrix(cns_wide[,1:3])
  bins[,4] <- complete.cases(segments)
  bins <- Biobase::AnnotatedDataFrame(as.data.frame(bins))
  bins@data$chromosome <- factor(bins@data$chromosome, levels = c(as.character(c(1:22)),'X'))
  bins@data$start <- as.integer(bins@data$start)
  bins@data$end <- as.integer(bins@data$end)
  bins@data$use <- as.logical(bins@data$use)

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

# New plotting function
setMethod("plot", signature(x="QDNAseqSignals", y="missing"),
          function (x, y, main=NULL, includeReadCounts=TRUE,
                    logTransform=TRUE, scale=TRUE, sdFUN="sdDiffTrim",
                    delcol=getOption("QDNAseq::delcol", "darkred"),
                    losscol=getOption("QDNAseq::losscol", "red"),
                    gaincol=getOption("QDNAseq::gaincol", "blue"),
                    ampcol=getOption("QDNAseq::ampcol", "darkblue"),
                    pointcol=getOption("QDNAseq::pointcol", "black"),
                    segcol=getOption("QDNAseq::segcol", "chocolate"),
                    misscol=getOption("QDNAseq::misscol", NA),
                    pointpch=getOption("QDNAseq::pointpch", 1L),
                    pointcex=getOption("QDNAseq::pointcex", 0.1),
                    xlab=NULL, ylab=NULL, ylim=NULL, xaxt="s", yaxp=NULL,
                    showDataPoints=TRUE, showSD=TRUE, doSegments=TRUE, doCalls=TRUE, ...,
                    verbose=getOption("QDNAseq::verbose", TRUE)) {
            
            oopts <- options("QDNAseq::verbose"=verbose)
            on.exit(options(oopts))
            
            ## Import private functions
            ns <- asNamespace("CGHbase")
            .getChromosomeLengths <- get(".getChromosomeLengths", envir=ns, mode="function")
            .makeSegments <- get(".makeSegments", envir=ns, mode="function")
            
            if (inherits(x, c("QDNAseqCopyNumbers", "QDNAseqReadCounts"))) {
              condition <- QDNAseq:::binsToUse(x)
            } else {
              condition <- rep(TRUE, times=nrow(x))
            }
            baseLine <- NA_real_
            doCalls <- "calls" %in% assayDataElementNames(x) & doCalls
            doSegments <- "segmented" %in% assayDataElementNames(x) & doSegments
            if (doCalls) {
              if (is.null(ylim))
                if (logTransform) {
                  ylim <- c(-5, 5)
                } else {
                  ylim <- c(-2, 4)
                }
            }
            if ("copynumber" %in% assayDataElementNames(x)) {
              copynumber <- assayDataElement(x, "copynumber")[condition, , drop=FALSE]
              if (is.null(ylab))
                ylab <- ifelse(logTransform, expression(log[2]~ratio), "ratio")
              if (is.null(ylim))
                if (logTransform) {
                  ylim <- c(-3, 5)
                } else {
                  ylim <- c(0, 4)
                }
              if (is.null(yaxp))
                yaxp <- c(ylim[1], ylim[2], ylim[2]-ylim[1])
              baseLine <- ifelse(logTransform, 0, 1)
            } else {
              copynumber <- assayDataElement(x, "counts")[condition, , drop=FALSE]
              if (is.null(ylab))
                ylab <- ifelse(logTransform, expression(log[2]~read~count),
                               "read count")
              if (is.null(ylim))
                if (logTransform) {
                  ylim <- c(0, max(QDNAseq:::log2adhoc(copynumber)))
                } else {
                  ylim <- range(copynumber)
                }
            }
            if (is.null(main))
              main <- sampleNames(x)
            if (includeReadCounts && "total.reads" %in% names(pData(x)))
              main <- paste(main, " (",
                            format(x$total.reads, trim=TRUE, big.mark=","), " reads)", sep="")
            if (length(ylab) == 1)
              ylab <- rep(ylab, times=ncol(x))
            all.chrom <- QDNAseq:::chromosomes(x)
            if (is.integer(all.chrom)) # when x is a cghRaw, cghSeg, or cghCall object
              all.chrom <- as.character(all.chrom)
            chrom <- all.chrom[condition]
            uni.chrom <- unique(chrom)
            uni.chrom <- as.character(str_sort(uni.chrom, numeric = TRUE))
            chrom.num <- as.integer(factor(chrom, levels= uni.chrom, ordered=TRUE))
            uni.chrom.num <- unique(chrom.num)
            uni.chrom.num <- sort( uni.chrom.num )
            if (!scale) {
              pos <- pos2 <- 1:sum(condition)
              chrom.ends <- aggregate(pos,
                                      by=list(chromosome=chrom), FUN=max)$x
            } else {
              if (inherits(x, c("cghRaw", "cghSeg", "cghCall"))) {
                chrom.lengths <- .getChromosomeLengths("GRCh37")
              } else {
                all.chrom.lengths <- aggregate(bpend(x),
                                               by=list(chromosome=all.chrom), FUN=max)
                chrom.lengths <- all.chrom.lengths$x
                names(chrom.lengths) <- all.chrom.lengths$chromosome
              }
              pos <- as.numeric(bpstart(x)[condition])
              pos2 <- as.numeric(bpend(x)[condition])
              chrom.lengths <- chrom.lengths[uni.chrom]
              chrom.ends <- integer()
              cumul <- 0
              for (i in seq_along(uni.chrom)) {
                pos[chrom.num > uni.chrom.num[i]] <-
                  pos[chrom.num > uni.chrom.num[i]] +
                  chrom.lengths[uni.chrom[i]]
                pos2[chrom.num > uni.chrom.num[i]] <-
                  pos2[chrom.num > uni.chrom.num[i]] +
                  chrom.lengths[uni.chrom[i]]
                cumul <- cumul + chrom.lengths[uni.chrom[i]]
                chrom.ends <- c(chrom.ends, cumul)
              }
              names(chrom.ends) <- names(chrom.lengths)
            }
            if (length(uni.chrom) == 1) {
              xax <- pretty(c(0, chrom.lengths[uni.chrom]))
              xaxlab <- xax / 1e6L
              if (is.null(xlab))
                xlab <- paste0("chromosome ", uni.chrom, ", Mbp")
            } else {
              xax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
              xaxlab <- uni.chrom
              if (is.null(xlab))
                xlab <- "chromosome"
            }
            if (inherits(x, c("cghRaw", "cghSeg", "cghCall")))
              copynumber <- log2adhoc(copynumber, inv=TRUE)
            if (is.character(sdFUN) && sdFUN == "sdDiffTrim") {
              symbol <- quote(hat(sigma)[Delta^"*"])
            } else if (is.character(sdFUN) && length(grep("Diff", sdFUN)) == 1) {
              symbol <- quote(hat(sigma)[Delta])
            } else {
              symbol <- quote(hat(sigma))
            }
            sdFUN <- QDNAseq:::sdDiffTrim
            noise <- apply(copynumber, MARGIN=2L, FUN=sdFUN, na.rm=TRUE)
            if (logTransform)
              copynumber <- QDNAseq:::log2adhoc(copynumber)
            for (i in seq_len(ncol(x))) {
              QDNAseq:::vmsg("Plotting sample ", main[i], " (", i, " of ", ncol(x), ") ...",
                             appendLF=FALSE)
              cn <- copynumber[, i]
              if (doSegments) {
                segmented <- assayDataElement(x, "segmented")[condition, i]
                if (inherits(x, c("cghRaw", "cghSeg", "cghCall")))
                  segmented <- QDNAseq:::log2adhoc(segmented, inv=TRUE)
                if (logTransform)
                  segmented <- QDNAseq:::log2adhoc(segmented)
                segment <- .makeSegments(segmented, chrom)
              }
              if (doCalls) {
                losses <- probloss(x)[condition, i]
                gains <- probgain(x)[condition, i]
                if (!is.null(probdloss(x)))
                  losses <- losses + probdloss(x)[condition, i]
                if (!is.null(probamp(x)))
                  gains <- gains + probamp(x)[condition, i]
                par(mar=c(5, 4, 4, 4) + 0.2)
                plot(NA, main=main[i], xlab=NA, ylab=NA, las=1,
                     xlim=c(0, max(chrom.ends)), ylim=ylim, xaxs="i", xaxt="n",
                     yaxp=c(ylim[1], ylim[2], ylim[2]-ylim[1]), yaxs="i")
                lim <- par("usr")
                lim[3:4] <- c(0, 1)
                par(usr=lim)
                dticks <- seq(0, 1, by=0.2)
                axis(4, at=dticks, labels=NA, srt=270, las=1, cex.axis=1,
                     cex.lab=1, tck=-0.015)
                axis(4, at=dticks, labels=dticks, srt=270, las=1, cex.axis=1,
                     cex.lab=1, line=-0.4, lwd=0)
                mtext("probability", side=4, line=2, cex=par("cex"))
                if (!is.na(misscol)) {
                  rect(0, -1, max(pos2), 1, col=misscol, border=NA)
                  rect(pos, -1, pos2, 1, col="white", border=NA)
                }
                rect(pos[segment[,2]], 0, pos2[segment[,3]], losses[segment[,2]],
                     col=losscol, border=losscol)
                if (!is.null(probdloss(x)))
                  rect(pos[segment[,2]], 0, pos2[segment[,3]],
                       probdloss(x)[condition, i][segment[,2]],
                       col=delcol, border=delcol)
                rect(pos[segment[,2]], 1, pos2[segment[,3]], 1-gains[segment[,2]],
                     col=gaincol, border=gaincol)
                if (!is.null(probamp(x)))
                  rect(pos[segment[,2]], 1, pos2[segment[,3]],
                       1-probamp(x)[condition, i][segment[,2]],
                       col=ampcol, border=ampcol)
                axis(3, at=pos[which(probamp(x)[condition,i] >= 0.5)],
                     labels=FALSE, col=ampcol, col.axis="black", srt=270, las=1,
                     cex.axis=1, cex.lab=1)
                axis(1, at=pos[which(probdloss(x)[condition,i] >= 0.5)],
                     labels=FALSE, col=delcol, col.axis="black", srt=270, las=1,
                     cex.axis=1, cex.lab=1)
                box()
                lim[3:4] <- ylim
                par(usr=lim)
                points(pos, cn, cex=pointcex, col=pointcol, pch=pointpch)
              } else {
                plot(pos, cn, cex=pointcex, col=pointcol, main=main[i],
                     xlab=NA, ylab=NA, ylim=ylim, xaxt="n", xaxs="i", yaxs="i",
                     yaxp=yaxp, tck=-0.015, las=1, pch=pointpch)
              }
              mtext(text=xlab, side=1, line=2, cex=par("cex"))
              mtext(text=ylab[i], side=2, line=2, cex=par("cex"))
              abline(h=baseLine)
              abline(v=chrom.ends[-length(chrom.ends)], lty="dashed")
              if (!is.na(xaxt) && xaxt != "n") {
                axis(side=1, at=xax, labels=NA, cex=.2, lwd=.5, las=1,
                     cex.axis=1, cex.lab=1, tck=-0.015)
                axis(side=1, at=xax, labels=xaxlab, cex=.2, lwd=0, las=1,
                     cex.axis=1, cex.lab=1, tck=-0.015, line=-0.4)
              }
              if (doSegments) {
                for (jjj in seq_len(nrow(segment))) {
                  segments(pos[segment[jjj,2]], segment[jjj,1],
                           pos[segment[jjj,3]], segment[jjj,1], col=segcol, lwd=3)
                }
              }
              par(xpd=TRUE)
              amps <- cn
              amps[amps <= ylim[2]] <- NA_real_
              amps[!is.na(amps)] <- ylim[2] + 0.01 * (ylim[2]-ylim[1])
              dels <- cn
              dels[dels >= ylim[1]] <- NA_real_
              dels[!is.na(dels)] <- ylim[1] - 0.01 * (ylim[2]-ylim[1])
              points(pos, amps, pch=24, col=pointcol, bg=pointcol, cex=0.5)
              points(pos, dels, pch=25, col=pointcol, bg=pointcol, cex=0.5)
              if (doSegments) {
                amps <- assayDataElement(x, "segmented")[condition, i]
                if (logTransform)
                  amps <- QDNAseq:::log2adhoc(amps)
                amps[amps <= ylim[2]] <- NA_real_
                amps[!is.na(amps)] <- ylim[2] + 0.01 * (ylim[2]-ylim[1])
                dels <- assayDataElement(x, "segmented")[condition, i]
                if (logTransform)
                  dels <- QDNAseq:::log2adhoc(dels)
                dels[dels >= ylim[1]] <- NA_real_
                dels[!is.na(dels)] <- ylim[1] - 0.01 * (ylim[2]-ylim[1])
                points(pos, amps, pch=24, col=segcol, bg=segcol, cex=0.5)
                points(pos, dels, pch=25, col=segcol, bg=segcol, cex=0.5)
              }
              par(xpd=FALSE)
              ### estimate for standard deviation
              if (showSD) {
                if (!is.na(x$expected.variance[i])) {
                  sdexp <- substitute(paste(E~sigma==e, ", ", symbol==sd),
                                      list(e=sprintf("%.3g", sqrt(x$expected.variance[i])),
                                           symbol=symbol, sd=sprintf("%.3g", noise[i])))
                } else {
                  sdexp <- substitute(symbol==sd,
                                      list(symbol=symbol, sd=sprintf("%.3g", noise[i])))
                }
                mtext(sdexp, side=3, line=0, adj=1, cex=par("cex"))
              }
              ### number of data points
              if (showDataPoints) {
                str <- paste(round(sum(condition) / 1000), "k x ", sep="")
                probe <- median(bpend(x)-bpstart(x)+1)
                if (probe < 1000) {
                  str <- paste(str, probe, " bp", sep="")
                } else {
                  str <- paste(str, round(probe / 1000), " kbp", sep="")
                }
                if (doSegments)
                  str <- paste(str, ", ", nrow(segment), " segments", sep="")
                mtext(str, side=3, line=0, adj=0, cex=par("cex"))
              }
              QDNAseq:::vmsg()
            }
            options("QDNAseq::plotLogTransform"=logTransform)
            options("QDNAseq::plotScale"=scale)
          })

# For loop I used to generate the plots
# Can optimize to remove this as well as other for loop
for (i in seq_len(262)) {
  png(filename = paste0( colnames(wx_qdnaobj)[i] ,"_segs_plot.png"), 
      height = 25 , width = 30, res = 300, units = "cm")
  plot(condor_file[,i], ylim=c(-3, 5), doSegment = TRUE, showSD = TRUE)
  dev.off()
}



