suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(magrittr)
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input_path"), type="character", default='input_bam_header.sam',
              help="Path to the bam header to be re-structured. [default %default]"),
  make_option(c("-o", "--output_path"), type="character", default='output_bam_header.sam',
              help="Path where you would like the new header written. [default %default]")
)
args <- parse_args(OptionParser(usage="%prog [options] file", option_list=option_list))


processFile = function(filepath, output_path) {
  con <- file(filepath, "r")
  fileConn <- file(output_path, open = 'a')
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    if (strsplit(line, '\t')[[1]][2] == 'SN:chrM') {
      line <- sub('chrM', 'MT', line)
    }
    
    line <- sub('chr', '', line)
    
    temp <- strsplit(line, '\t')[[1]]
    if (grepl("_", temp[2], fixed = TRUE)) {
      next
    }
    
    writeLines(line, con = fileConn)
  }
  close(fileConn)
}

processFile(args$input_path, args$output_path)
