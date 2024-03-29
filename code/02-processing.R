
# Packages ----------------------------------------------------------------

library(dplyr)

# Data --------------------------------------------------------------------
## Reading the mapping IDs
mapping_ids <- read.table(here::here("data/id_mapping_icgc2tcga.45_donors.txt"), header = TRUE)

## Reading the IDs of interest
list_ids_icgc <- scan(here::here("data/ids_of_interest_icgc.txt"), character(), quote = "")
list_ids_tcga <- scan(here::here("data/ids_of_interest_tcga.txt"), character(), quote = "")

## Reading overall data
id_with_pcawg_overview <- read.csv(here::here("data/id_and_response_with_pcawg_overview.csv"), row.names = 1)

## Reading CNV files for ids of interest
all_CNV_files <- list.files(path = "data/consensus_cnv/consensus.20170119.somatic.cna.annotated", pattern = "*.txt")
interest_CNV_files <- subset(all_CNV_files, grepl(paste0(list_ids_icgc, collapse = "|"),
                                                  all_CNV_files))

## Creating selected IDs to put back into dataframe
selected_icgc_ids <- sapply(strsplit(interest_CNV_files, split = ".", fixed = TRUE), function(x) (x[1]))


### For reading all the CNV files into a list (Use this)
CNV_files_list <- lapply(paste0(here::here("data/consensus_cnv/consensus.20170119.somatic.cna.annotated/"), interest_CNV_files), 
                         read.table, header = TRUE, stringsAsFactors = FALSE)

## Reading SV files for ids of interest

# all_SV_files <- list.files(path = "data/consensus_sv/tcga/open", pattern = "*.gz")
# interest_SV_files <- subset(all_SV_files, grepl(paste0(list_ids_icgc, collapse = "|"),
#                                                 all_SV_files))


### For reading all the SV files into a list (Use this)
# SV_files_list <- lapply(paste0(here::here("data/consensus_sv/tcga/open/"), interest_SV_files), 
#                         read.table, header = TRUE, stringsAsFactors = FALSE)


## Reading gene level calls in CNV
### CN gene level calls
cnv_cn_gene_level_calls <- read.table(here::here("data/consensus_cnv/gene_level_calls/all_samples.consensus_CN.by_gene.170214.txt"),
                                        fill = TRUE, header = TRUE)

### Consensus gene level calls
cnv_consensus_gene_level_calls <- read.table(here::here("data/consensus_cnv/gene_level_calls/all_samples.consensus_level_calls.by_gene.170214.txt"), 
                                             fill = TRUE, header = TRUE)

## Reading pan-cancer compendium CN signatures for PCAWG
pan_cancer_cn_signatures_pcawg <- read.csv(here::here("data/chromosomal_instability_data/PCAWG_activities_scaled.csv"), row.names = 1)

### Subset selected ids for pan-cancer compendium CN signatures for PCAWG
pan_cancer_cn_signatures_pcawg_selected <- subset(pan_cancer_cn_signatures_pcawg, grepl(paste0(list_ids_icgc, collapse = "|"),
                                                                            rownames(pan_cancer_cn_signatures_pcawg)))

## Reading pan-cancer compendium CN signautres for tcga
## Note: Use the PCAWG dataset as it uses WGS
# pan_cancer_cn_signatures_tcga <- read.csv(here::here("data/chromosomal_instability_data/tcga_activities_scaled.csv"), row.names = 1)

### Subset selected ids for pan-cancer compendium CN signatures for TCGA
# pan_cancer_cn_signatures_tcga_selected <- subset(pan_cancer_cn_signatures_tcga, grepl(paste0(list_ids_tcga, collapse = "|"),
#                                                                                       rownames(pan_cancer_cn_signatures_tcga)))

## Reading whole genome duplication for PCAWG from evolution and heterogeneity
wgd_pcawg_evolution <- read.table(here::here("data/evolution_and_heterogeneity/2018-07-24-wgdMrcaTiming.txt"), header = TRUE)

## Reading SNV gene level calls
utr3_matrix <- read.csv(here::here("data/snv_gene_level_calls/gene2sample_mutation_matrix.3UTR.csv"), row.names = 1)
utr5_matrix <- read.csv(here::here("data/snv_gene_level_calls/gene2sample_mutation_matrix.5UTR.csv"), row.names = 1)
high_moderate_matrix <- read.csv(here::here("data/snv_gene_level_calls/gene2sample_mutation_matrix.High_and_Moderate.csv"), row.names = 1)
ncTx_exon_matrix <- read.csv(here::here("data/snv_gene_level_calls/gene2sample_mutation_matrix.ncTx_exon.csv"), row.names = 1)
NMD_matrix <- read.csv(here::here("data/snv_gene_level_calls/gene2sample_mutation_matrix.NMD.csv"), row.names = 1)
upstream_2kb_matrix <- read.csv(here::here("data/snv_gene_level_calls/gene2sample_mutation_matrix.upstream_2kb.csv"), row.names = 1)


# Functions for other processing methods ---------------------------------------------------------------

filtering_snv_genes_interest <- function(snv_dataframe, ids_subset, genes_subset){
  # Function to filter out genes and patients of interest in a patient by gene matrix for SNV data. Also converts it to a gene by patient matrix.
  # snv_dataframe: Input your dataframe
  # ids_subset: Input the ids you want to subset
  # genes_subset: Input the genes you want to subset
  snv_selected_patients <- snv_dataframe[, grepl(paste(ids_subset, collapse = "|"),
                                                 names(snv_dataframe),)]
  snv_selected_genes_patients <- subset(snv_selected_patients, rownames(snv_selected_patients) %in% genes_subset)
  snv_selected_genes_patients_t <- data.frame(t(snv_selected_genes_patients))
  return(snv_selected_genes_patients_t)
}


## To generate total copy number for each sample, including major and minor DO NOT USE
# generate_total_cn_mut_load <- function(CNV_dataframe)
# {
#   total_cn <- sum(CNV_dataframe$total_cn, na.rm = TRUE)
#   major_cn <- sum(CNV_dataframe$major_cn, na.rm = TRUE)
#   minor_cn <- sum(CNV_dataframe$minor_cn, na.rm = TRUE)
#   out <- list("total_cn" = total_cn, "major_cn" = major_cn, "minor_cn" = minor_cn)
#   return(out)
# }


# SV Classes --------------------------------------------------------------

## To obtain all sv classes in dataset
# unique(unlist(lapply(SV_files_list, `[[`, "svclass")))

# CNV processing absolute copy numbers DO NOT USE --------------------------------------------------------------
# ## Adding up copy numbers for each sample
# cn_mut_load <- lapply(CNV_files_list, generate_total_cn_mut_load)
# # 
# # ## Converting the list of lists to one dataframe
# mut_load_dataframe <- data.frame(do.call(rbind, cn_mut_load))
# mut_load_dataframe <- data.frame(apply(mut_load_dataframe, 2, unlist))
# rownames(mut_load_dataframe) <- selected_icgc_ids # Adding id rownames


# Filtering Gene Level Calls for CNV data --------------------------------------------------------
## Columns to keep

list_ids_icgc_sub_sep <- gsub("-", "\\.", list_ids_icgc) # This replacement needs to be done due to different patient ID names
columns_keep_gene_level <- append(list_ids_icgc_sub_sep, c("Gene", "Symbol", "Locus", "ID", "Cytoband"))

## Processing CN by gene: Subsetting dataframe for patients of interest
cnv_cn_gene_level_calls_selected <- cnv_cn_gene_level_calls[, grepl(paste(columns_keep_gene_level, collapse = "|"),
                                                        names(cnv_cn_gene_level_calls),)]

cnv_cn_gene_level_calls_colnames_corrected <- gsub("\\.", "-", colnames(cnv_cn_gene_level_calls_selected)) # Replacement to revert back to original id names
colnames(cnv_cn_gene_level_calls_selected) <- cnv_cn_gene_level_calls_colnames_corrected

### Removing unnecessary columns and setting gene as rownames
drop_columns <- c("Symbol", "Locus", "ID", "Cytoband")
cnv_cn_gene_level_calls_selected <- cnv_cn_gene_level_calls_selected[, !names(cnv_cn_gene_level_calls_selected) %in% drop_columns]
rownames(cnv_cn_gene_level_calls_selected) <- cnv_cn_gene_level_calls_selected[, 1]
cnv_cn_gene_level_calls_selected[, 1] <- NULL

## Processing level calls by gene

calls_by_gene_selected <- cnv_consensus_gene_level_calls[, grepl(paste(columns_keep_gene_level, collapse = "|"),
                                                          names(cnv_consensus_gene_level_calls),)]

colnames(calls_by_gene_selected) <- cnv_cn_gene_level_calls_colnames_corrected

### Removing unnecessary columns and setting gene as rownames
calls_by_gene_selected<- calls_by_gene_selected[, !names(calls_by_gene_selected) %in% drop_columns]
rownames(calls_by_gene_selected) <- calls_by_gene_selected[, 1]
calls_by_gene_selected[, 1] <- NULL

## Output files for gene level calls
# write.csv(cnv_cn_gene_level_calls_selected, file = "data/consensus_cnv/gene_level_calls/selected_cn_by_gene.csv", 
#             quote = FALSE)

# write.csv(calls_by_gene_selected, file = "data/consensus_cnv/gene_level_calls/selected_calls_by_gene.csv",
#           quote = FALSE)


## Filtering out genes of interest
list_genes_selected <- c("BRCA1", "BRCA2", "ATM", "BARD1",
                         "BRIP1", "CHEK1", "CHEK2", "FAM175A",
                         "MRE11A", "NBN", "PALB2", "RAD51C", "RAD51D")

calls_by_gene_selected_genes_interest <- subset(calls_by_gene_selected, rownames(calls_by_gene_selected) %in% list_genes_selected)
calls_by_gene_selected_genes_interest_t <- data.frame(t(calls_by_gene_selected_genes_interest)) # Transposing rows and columns

## Removing X in rownames
rownames(calls_by_gene_selected_genes_interest_t) <- sub("X*", "", rownames(calls_by_gene_selected_genes_interest_t))



# Filtering Gene level calls for SNV data ---------------------------------
list_ids_tcga_sub_sep <- gsub("-", "\\.", list_ids_tcga) # This replacement needs to be done due to different patient ID names

## Subsetting relevant patient ids and genes for SNV data using function

utr3_selected_genes_interest <- filtering_snv_genes_interest(utr3_matrix, list_ids_tcga_sub_sep, list_genes_selected)
utr5_selected_genes_interest <- filtering_snv_genes_interest(utr5_matrix, list_ids_tcga_sub_sep, list_genes_selected)
high_moderate_selected_genes_interest <- filtering_snv_genes_interest(high_moderate_matrix, list_ids_tcga_sub_sep, list_genes_selected)
ncTx_exon_selected_genes_interest <- filtering_snv_genes_interest(ncTx_exon_matrix, list_ids_tcga_sub_sep, list_genes_selected)
NMD_selected_genes_interest <- filtering_snv_genes_interest(NMD_matrix, list_ids_tcga_sub_sep, list_genes_selected)
upstream_2kb_selected_genes_interest <- filtering_snv_genes_interest(upstream_2kb_matrix, list_ids_tcga_sub_sep, list_genes_selected)

## Converting ID separators to -
rownames(utr3_selected_genes_interest) <- gsub("\\.", "-", rownames(utr3_selected_genes_interest))
rownames(utr5_selected_genes_interest) <- gsub("\\.", "-", rownames(utr5_selected_genes_interest))
rownames(high_moderate_selected_genes_interest) <- gsub("\\.", "-", rownames(high_moderate_selected_genes_interest))
rownames(ncTx_exon_selected_genes_interest) <- gsub("\\.", "-", rownames(ncTx_exon_selected_genes_interest))
rownames(NMD_selected_genes_interest) <- gsub("\\.", "-", rownames(NMD_selected_genes_interest))
rownames(upstream_2kb_selected_genes_interest) <- gsub("\\.", "-", rownames(upstream_2kb_selected_genes_interest))


# Whole Genome Duplication ------------------------------------------------

wgd_pcawg_evolution_selected <- subset(wgd_pcawg_evolution, grepl(paste0(list_ids_icgc, collapse = "|"),
                                                                             wgd_pcawg_evolution$uuid))

# Combining into one dataframe** --------------------------------------------

final_dataframe <- id_with_pcawg_overview[, c("tumour_specimen_aliquot_id","Condition", "age")]

## Merging overview data
final_dataframe <- merge(final_dataframe, pan_cancer_cn_signatures_pcawg_selected, 
                         by.x = "tumour_specimen_aliquot_id", by.y = "row.names", all.x = TRUE)
final_dataframe[2, 3] <- 50 # Add age value to NA. This value was taken from evolution and heterogeneity dataset

## Merging gene level calls data
final_dataframe <- merge(final_dataframe, calls_by_gene_selected_genes_interest_t,
                         by.x = "tumour_specimen_aliquot_id", by.y = "row.names", all.x = TRUE)

## Merging WGD data
final_dataframe <- merge(final_dataframe, wgd_pcawg_evolution_selected[, c("uuid", "WGD")],
                         by.x = "tumour_specimen_aliquot_id", by.y = "uuid", all.x = TRUE)

## Adding TCGA IDs and reordering them
final_dataframe <- merge(final_dataframe, mapping_ids[, c("tumour_specimen_aliquot_id", "TCGA_id")],
                         by = "tumour_specimen_aliquot_id")

final_dataframe <- final_dataframe %>%
  relocate(TCGA_id, .before = Condition)

## Adding high and moderate impact variants to the dataframe
### Renaming columns for high and moderate impact variants
colnames(high_moderate_selected_genes_interest) <- paste("mh", colnames(high_moderate_selected_genes_interest), sep = "")

## Merging high and moderate impact variants
final_dataframe <- merge(final_dataframe, high_moderate_selected_genes_interest,
                         by.x = "TCGA_id", by.y = "row.names", all.x = TRUE)

## Finalising dataframe by making the IDs to row names
# rownames(final_dataframe) <- final_dataframe$tumour_specimen_aliquot_id
# final_dataframe <- subset(final_dataframe, select = -c(tumour_specimen_aliquot_id))

## Output final dataframe
# write.csv(final_dataframe, "data/final_dataframe.csv")

# Initial Analysis --------------------------------------------------------
# hist(final_dataframe$age)
# barplot(summary(final_dataframe$Condition))
# barplot(total_cn ~ Row.names, data = final_dataframe, las = 2, axisnames = FALSE)