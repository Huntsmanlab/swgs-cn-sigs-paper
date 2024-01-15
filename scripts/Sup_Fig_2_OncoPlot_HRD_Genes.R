# Create an oncoplot
# Use maftools

suppressPackageStartupMessages({
  library(maftools)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(readr)
  library(magrittr)
  library(tidyr)
  library(patchwork)
})

# Create MAF file from custom xls for p53abn cn-sigs project
maff <- readxl::read_xls('./mutation_calls/JSB_reclass_finalmutationcalls_p53AbNcohort.xls')

TPanel_list_samples <- readxl::read_xlsx("./clinical_data/TPanel_list_samples.xlsx")
TPanel_list_samples <- TPanel_list_samples[-which(TPanel_list_samples$sample_id == 'VOA888A'), ]
samples_of_interest <- data.frame(TPanel_list_samples$sample_id)
colnames(samples_of_interest) <- 'Tumor_Sample_Barcode'

Hugo_Symbol <- c("ATM","ATRX","BAP1","BLM","BRCA1","BRCA2","BRIP1",
                 "CHEK1","CHEK2","EXO1","FANCA","FANCC","FANCD2","FANCE",
                 "FANCF","FANCG","MRE11A","PALB2","RAD50","RAD51")
Hugo_Symbol <- data.frame(Hugo_Symbol)

colnames(maff)[1] <- 'sample_id'
colnames(maff)[2] <- 'chromosome'
colnames(maff)[3] <- 'start'
maff$end <- maff$start + nchar(maff$REF) - 1
colnames(maff)[9] <- 'ref'
colnames(maff)[10] <- 'alt'
maff <- maff[,-c(13)]
maff <- maff[,-c(20)]
colnames(maff)[20] <- 'consequence'
colnames(maff)[22] <- 'gene_symbol'

maff <- maff[,c(1:3, 37, 4:36)]

maff <- maff %>% dplyr::rename(Hugo_Symbol = gene_symbol,
                               Entrez_Gene_Id = Gene_ID,
                               Chromosome = chromosome,
                               Start_Position = start,
                               End_Position = end,
                               Variant_Classification = consequence,
                               Variant_Type = Oncoplot_Variant_Type,
                               Reference_Allele = ref,
                               Tumor_Seq_Allele2 = alt,
                               Tumor_Sample_Barcode = sample_id,
                               n_depth = DP)
maff$Center <- 'GSC'
maff$NCBI_Build <- 37
maff$Strand <- '+'

maff <- Hugo_Symbol %>% left_join(maff)
maff <- samples_of_interest %>% left_join(maff)
maff$Variant_Classification[which(maff$Variant_Classification == 'Unknown')] <- 'InDel'

maff = read.maf(maf = maff, vc_nonSyn = c("5'UTR", "Frame_Shift_Del",
                                          "Frame_Shift_Ins", "In_Frame_Del",
                                          "In_Frame_Ins", "Missense_Mutation",
                                          "Nonsense_Mutation", "Nonstop_Mutation",
                                          "Splice_Site", "Translation_Start_Site",
                                          "InDel"))

pdf(file = "~/Downloads/Oncoplot_HRD.pdf")

pp <- oncoplot(maf = maff, draw_titv = FALSE, top = 20, drawColBar = FALSE)

dev.off()
