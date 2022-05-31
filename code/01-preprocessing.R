
# Data --------------------------------------------------------------------

cna_annotated <- read.table(here::here("data/consensus_cnv/consensus.20170119.somatic.cna.annotated/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20170119.somatic.cna.annotated.txt"), header = TRUE)
  # For one individual's data
consensus_cn_gene <- read.table(here::here("data/consensus_cnv/gene_level_calls/all_samples.consensus_CN.by_gene.170214.txt"), fill = TRUE)
platinum_response <- read.table(here::here("data/TCGA_OV_PlatinumResponse_PMID27526849.txt"), header = TRUE)
id_mapping <- read.table(here::here("data/id_mapping_icgc2tcga.45_donors.txt"), header = TRUE)


# Processing Steps --------------------------------------------------------
## Merging platinum response and the IDs of ICGC and TCGA 45 donors together
id_with_response <- merge(platinum_response, id_mapping, by.x = "SampleCode", by.y = "TCGA_id")

