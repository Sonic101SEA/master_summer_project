
# Data --------------------------------------------------------------------

## For one individual's copy number variation data
cna_annotated <- read.table(here::here("data/consensus_cnv/consensus.20170119.somatic.cna.annotated/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20170119.somatic.cna.annotated.txt"), header = TRUE)

## For gene data
### Major
consensus_cn_major_gene <- read.table(here::here("data/consensus_cnv/gene_level_calls/all_samples.consensus_CN.by_gene.170214.txt"), fill = TRUE, header = TRUE)

### Level calls
consensus_cn_level <- read.table(here::here("data/consensus_cnv/gene_level_calls/all_samples.consensus_level_calls.by_gene.170214.txt"), fill = TRUE, header = TRUE)

### GISTIC Analysis

all_data_by_genes <- read.table(here::here("data/consensus_cnv/GISTIC_analysis/all_data_by_genes.rmcnv.pt_170207.txt"), fill = TRUE, header = TRUE)

## Data mapping
platinum_response <- read.table(here::here("data/TCGA_OV_PlatinumResponse_PMID27526849.txt"), header = TRUE)
id_mapping <- read.table(here::here("data/id_mapping_icgc2tcga.45_donors.txt"), header = TRUE)

## Structural variant data
### For one individual's structural variant data
sv_data <- read.table(here::here("data/consensus_sv/icgc/open/0a6be23a-d5a0-4e95-ada2-a61b2b5d9485.pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz"), header = TRUE)

## Overiew Data
overview_data <- readxl::read_excel(here::here("data/pcawg_overview_data.xlsx"))



# Processing Steps --------------------------------------------------------
## Merging platinum response and the IDs of ICGC and TCGA 45 donors together
id_with_response <- merge(platinum_response, id_mapping, by.x = "SampleCode", by.y = "TCGA_id")

## Merging overview data with id_with_response
id_with_response_overview <- merge(id_with_response, overview_data, by = "tumour_specimen_aliquot_id", all.x = TRUE)

## Names of IDs with response
ids <- as.character(id_with_response_overview$tumour_specimen_aliquot_id)
ids_X <- paste("X", ids, sep = "")

## Processing consensus_cn_gene
### Keeping all columns seen in ids
consensus_cn_gene_selected <- consensus_cn_gene[, grepl(paste(ids_X, collapse = "|"),
                                                        names(consensus_cn_gene),)]
  # Work in progress
consensus_cn_gene_selected <- as.data.frame(consensus_cn_gene_selected)


# Output files ------------------------------------------------------------

## Output id and response with pcawg overview
# write.csv(id_with_response_overview, file = "data/id_and_response_with_pcawg_overview.csv")
## Output id and response only
# write.csv(id_with_response, file = "data/id_and_therapy_response.csv")
## Output id only
# write.table(ids, file = "data/ids.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
  # Note: The ids here are the 42 donors we are interested in


