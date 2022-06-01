
# Data --------------------------------------------------------------------

## For one individual's copy number variation data
cna_annotated <- read.table(here::here("data/consensus_cnv/consensus.20170119.somatic.cna.annotated/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20170119.somatic.cna.annotated.txt"), header = TRUE)

## For combined data
consensus_cn_gene <- read.table(here::here("data/consensus_cnv/gene_level_calls/all_samples.consensus_CN.by_gene.170214.txt"), fill = TRUE)

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

## Processing consensus_cn_gene
