suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(magrittr)
  library(data.table)
  library(vcfR)
  library(VariantAnnotation)
  library(optparse)
})


option_list <- list(
  make_option(c("-d", "--depth_limit"), type="integer", default=250,
              help="Filter variants with a read count below this number [default %default]"),
  make_option(c("-v", "--vaf_min"), type="double", default=0.1,
              help="The minimum variant allele frequency. Filter variants that occur less frequently than this. [default %default]"),
  make_option(c("-o", "--output_path"), type="character", default='variants_after_depth_filtering.tsv',
              help="The path where you would like th filtered variants written. [default %default]")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)

args <- parse_args(parser, positional_arguments = 1)
opt <- args$options
vfile <- args$args

vcff <- read.vcfR(vfile, verbose = FALSE )
vcf2 <- readVcf(vfile, "hg19")

# Move funcotator output to table
funcotator_fields <- 'Gencode_34_hugoSymbol|Gencode_34_ncbiBuild|Gencode_34_chromosome|Gencode_34_start|Gencode_34_end|Gencode_34_variantClassification|Gencode_34_secondaryVariantClassification|Gencode_34_variantType|Gencode_34_refAllele|Gencode_34_tumorSeqAllele1|Gencode_34_tumorSeqAllele2|Gencode_34_genomeChange|Gencode_34_annotationTranscript|Gencode_34_transcriptStrand|Gencode_34_transcriptExon|Gencode_34_transcriptPos|Gencode_34_cDnaChange|Gencode_34_codonChange|Gencode_34_proteinChange|Gencode_34_gcContent|Gencode_34_referenceContext|Gencode_34_otherTranscripts|Achilles_Top_Genes|CGC_Name|CGC_GeneID|CGC_Chr|CGC_Chr_Band|CGC_Cancer_Somatic_Mut|CGC_Cancer_Germline_Mut|CGC_Tumour_Types__(Somatic_Mutations)|CGC_Tumour_Types_(Germline_Mutations)|CGC_Cancer_Syndrome|CGC_Tissue_Type|CGC_Cancer_Molecular_Genetics|CGC_Mutation_Type|CGC_Translocation_Partner|CGC_Other_Germline_Mut|CGC_Other_Syndrome/Disease|ClinVar_HGMD_ID|ClinVar_SYM|ClinVar_TYPE|ClinVar_ASSEMBLY|ClinVar_rs|ClinVar_VCF_AF_ESP|ClinVar_VCF_AF_EXAC|ClinVar_VCF_AF_TGP|ClinVar_VCF_ALLELEID|ClinVar_VCF_CLNDISDB|ClinVar_VCF_CLNDISDBINCL|ClinVar_VCF_CLNDN|ClinVar_VCF_CLNDNINCL|ClinVar_VCF_CLNHGVS|ClinVar_VCF_CLNREVSTAT|ClinVar_VCF_CLNSIG|ClinVar_VCF_CLNSIGCONF|ClinVar_VCF_CLNSIGINCL|ClinVar_VCF_CLNVC|ClinVar_VCF_CLNVCSO|ClinVar_VCF_CLNVI|ClinVar_VCF_DBVARID|ClinVar_VCF_GENEINFO|ClinVar_VCF_MC|ClinVar_VCF_ORIGIN|ClinVar_VCF_RS|ClinVar_VCF_SSR|ClinVar_VCF_ID|ClinVar_VCF_FILTER|Cosmic_overlapping_mutations|CosmicFusion_fusion_genes|CosmicFusion_fusion_id|CosmicTissue_total_alterations_in_gene|CosmicTissue_tissue_types_affected|DNARepairGenes_Activity_linked_to_OMIM|DNARepairGenes_Chromosome_location_linked_to_NCBI_MapView|DNARepairGenes_Accession_number_linked_to_NCBI_Entrez|Familial_Cancer_Genes_Syndrome|Familial_Cancer_Genes_Synonym|Familial_Cancer_Genes_Reference|Gencode_XHGNC_hgnc_id|Gencode_XRefSeq_mRNA_id|Gencode_XRefSeq_prot_acc|HGNC_HGNC_ID|HGNC_Approved_Name|HGNC_Status|HGNC_Locus_Type|HGNC_Locus_Group|HGNC_Previous_Symbols|HGNC_Previous_Name|HGNC_Synonyms|HGNC_Name_Synonyms|HGNC_Chromosome|HGNC_Date_Modified|HGNC_Date_Symbol_Changed|HGNC_Date_Name_Changed|HGNC_Accession_Numbers|HGNC_Enzyme_IDs|HGNC_Entrez_Gene_ID|HGNC_Ensembl_Gene_ID|HGNC_Pubmed_IDs|HGNC_RefSeq_IDs|HGNC_Gene_Family_ID|HGNC_Gene_Family_Name|HGNC_CCDS_IDs|HGNC_Vega_ID|HGNC_Entrez_Gene_ID(supplied_by_NCBI)|HGNC_OMIM_ID(supplied_by_OMIM)|HGNC_RefSeq(supplied_by_NCBI)|HGNC_UniProt_ID(supplied_by_UniProt)|HGNC_Ensembl_ID(supplied_by_Ensembl)|HGNC_UCSC_ID(supplied_by_UCSC)|Oreganno_Build|Oreganno_ID|Oreganno_Values|Simple_Uniprot_uniprot_entry_name|Simple_Uniprot_DrugBank|Simple_Uniprot_alt_uniprot_accessions|Simple_Uniprot_uniprot_accession|Simple_Uniprot_GO_Biological_Process|Simple_Uniprot_GO_Cellular_Component|Simple_Uniprot_GO_Molecular_Function|dbSNP_ASP|dbSNP_ASS|dbSNP_CAF|dbSNP_CDA|dbSNP_CFL|dbSNP_COMMON|dbSNP_DSS|dbSNP_G5|dbSNP_G5A|dbSNP_GENEINFO|dbSNP_GNO|dbSNP_HD|dbSNP_INT|dbSNP_KGPhase1|dbSNP_KGPhase3|dbSNP_LSD|dbSNP_MTP|dbSNP_MUT|dbSNP_NOC|dbSNP_NOV|dbSNP_NSF|dbSNP_NSM|dbSNP_NSN|dbSNP_OM|dbSNP_OTH|dbSNP_PM|dbSNP_PMC|dbSNP_R3|dbSNP_R5|dbSNP_REF|dbSNP_RS|dbSNP_RSPOS|dbSNP_RV|dbSNP_S3D|dbSNP_SAO|dbSNP_SLO|dbSNP_SSR|dbSNP_SYN|dbSNP_TOPMED|dbSNP_TPA|dbSNP_U3|dbSNP_U5|dbSNP_VC|dbSNP_VLD|dbSNP_VP|dbSNP_WGT|dbSNP_WTD|dbSNP_dbSNPBuildID|dbSNP_ID|dbSNP_FILTER'
func <- str_split(funcotator_fields, pattern = '\\|')
funcotated <- matrix(NA_integer_, nrow=length(vcf2@info@listData[["FUNCOTATION"]]), ncol=length(func[[1]]),
                     dimnames=list(1:length(vcf2@info@listData[["FUNCOTATION"]]), func[[1]]))
for (i in 1:length(vcf2@info@listData[["FUNCOTATION"]])) {
  funcotated[i,] <- str_split(info(vcf2)$FUNCOTATION[[i]], pattern = '\\|')[[1]]
}


vars <- as.data.frame(vcff@fix)
vars <- cbind(vars, funcotated)
vars <- vars[vars$FILTER == 'PASS',]
vars <- vars %>% separate(INFO, c("AS_FilterStatus", "AS_SB_TABLE", "DP", "ECNT",
                                  "FUNCOTATION", "GERMQ", "MBQ", "MFRL", "MMQ", 
                                  "MPOS", "POPAF", "ROQ", "RPA", "RU", "STR", 
                                  "STRQ", "TLOD" ), ";")
vars$DP <- sub('DP=', '', vars$DP)
vars$DP <- as.integer(vars$DP)
vars <- vars %>% dplyr::filter(DP > opt$depth_limit)
vars$Gencode_34_hugoSymbol <- sub('\\[', '', vars$Gencode_34_hugoSymbol)
temp <- sub('AS_SB_TABLE=', '', vars$AS_SB_TABLE)
vars$vaf <- sub('\\|', ',', temp) 
vars <- vars %>%  separate(vaf, c("ref_pos_count", "ref_neg_count", 
                                  "alt_pos_count", "alt_neg_count"), ",", convert = TRUE)
vars$vaf <- (vars$alt_pos_count + vars$alt_neg_count)/(vars$ref_pos_count + vars$ref_neg_count)
output <- vars %>% dplyr::filter(vaf > opt$vaf_min)
write.table(output, file = opt$output_path, sep = '\t', row.names = FALSE)


# vcf <- read.vcfR('~/Downloads/CC-CHM-1347_mutectRun/CC-CHM-1347.qiaseq.filtered.vcf', verbose = FALSE )
# variants <- as.data.frame(vcf@fix)
# variants <- variants[variants$FILTER == 'PASS',]
# variants <- variants %>% separate(INFO, c("AS_FilterStatus", "AS_SB_TABLE", "DP", "ECNT", "GERMQ", 
#                    "MBQ", "MFRL", "MMQ", "MPOS", "POPAF", "ROQ", 
#                    "RPA", "RU", "STR", "STRQ", "TLOD" ), ";")
# variants$DP <- sub('DP=', '', variants$DP)
# variants$DP <- as.integer(variants$DP)
# variants <- variants %>% dplyr::filter(DP > 300)
#                    
# chrom <- masker(vcf, min_QUAL=0, min_DP=300, max_DP=1000, min_MQ=59.5, max_MQ=60.5)
# 
# 
# 
# funcotated <- data.frame()