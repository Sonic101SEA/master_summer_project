
# Data --------------------------------------------------------------------

cna_annotated <- read.table(here::here("data/consensus_cnv/consensus.20170119.somatic.cna.annotated/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20170119.somatic.cna.annotated.txt"), header = TRUE)
  # For one individual's data
consensus_cn_gene <- read.table(here::here("data/consensus_cnv/gene_level_calls/all_samples.consensus_CN.by_gene.170214.txt"), fill = TRUE)
platinum_response <- read.table(here::here("data/TCGA_OV_PlatinumResponse_PMID27526849.txt"), header = TRUE)


# Processing Steps --------------------------------------------------------


