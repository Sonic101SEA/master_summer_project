
# Data --------------------------------------------------------------------
list_ids <- scan(here::here("data/ids.txt"), character(), quote = "")

## Reading CNV files for ids of interest
all_CNV_files <- list.files(path = "data/consensus_cnv/consensus.20170119.somatic.cna.annotated", pattern = "*.txt")
interest_CNV_files <- subset(all_CNV_files, grepl(paste0(list_ids, collapse = "|"),
                                                  all_CNV_files))

for (i in 1:length(interest_CNV_files)) {
  assign(interest_CNV_files[i],
         read.table(paste0(here::here("data/consensus_cnv/consensus.20170119.somatic.cna.annotated/"),
                      interest_CNV_files[i]), header = TRUE))
}

## Reading SV files for ids of interest

all_SV_files <- list.files(path = "data/consensus_sv/tcga/open", pattern = "*.gz")
interest_SV_files <- subset(all_SV_files, grepl(paste0(list_ids, collapse = "|"),
                                                  all_SV_files))

for (i in 1:length(interest_SV_files)) {
  assign(interest_SV_files[i],
         read.table(paste0(here::here("data/consensus_sv/tcga/open/"),
                           interest_SV_files[i]), header = TRUE))
}

# Processing --------------------------------------------------------------
## Processing consensus_cn_gene
### Keeping all columns seen in ids
consensus_cn_gene_selected <- consensus_cn_gene[, grepl(paste(ids_X, collapse = "|"),
                                                        names(consensus_cn_gene),)]
# Work in progress
consensus_cn_gene_selected <- as.data.frame(consensus_cn_gene_selected)



