
# Data --------------------------------------------------------------------
list_ids <- scan(here::here("data/ids_of_interest.txt"), character(), quote = "")

## Reading CNV files for ids of interest
all_CNV_files <- list.files(path = "data/consensus_cnv/consensus.20170119.somatic.cna.annotated", pattern = "*.txt")
interest_CNV_files <- subset(all_CNV_files, grepl(paste0(list_ids, collapse = "|"),
                                                  all_CNV_files))

### For creating multiple objects
for (i in 1:length(interest_CNV_files)) {
  assign(interest_CNV_files[i],
         read.table(paste0(here::here("data/consensus_cnv/consensus.20170119.somatic.cna.annotated/"),
                      interest_CNV_files[i]), header = TRUE))
}

### For reading all the CNV files into a list (Use this)
CNV_files_list <- lapply(paste0(here::here("data/consensus_cnv/consensus.20170119.somatic.cna.annotated/"), interest_CNV_files), 
                         read.table, header = TRUE, stringsAsFactors = FALSE)

## Reading SV files for ids of interest

all_SV_files <- list.files(path = "data/consensus_sv/tcga/open", pattern = "*.gz")
interest_SV_files <- subset(all_SV_files, grepl(paste0(list_ids, collapse = "|"),
                                                  all_SV_files))

### For creating multiple objects
for (i in 1:length(interest_SV_files)) {
  assign(interest_SV_files[i],
         read.table(paste0(here::here("data/consensus_sv/tcga/open/"),
                           interest_SV_files[i]), header = TRUE))
}

### For reading all the SV files into a list (Use this)
SV_files_list <- lapply(paste0(here::here("data/consensus_sv/tcga/open/"), interest_SV_files), 
                         read.table, header = TRUE, stringsAsFactors = FALSE)


## Reading gene level calls in CNV
### Major
cnv_major_gene_level_calls <- read.table(here::here("data/consensus_cnv/gene_level_calls/all_samples.consensus_CN.by_gene.170214.txt"), 
                                   fill = TRUE, header = TRUE)



# Functions ---------------------------------------------------------------



# Processing --------------------------------------------------------------
## Processing gene level calls
### Keeping all columns seen in ids
consensus_cn_gene_selected <- consensus_cn_gene[, grepl(paste(ids_X, collapse = "|"),
                                                        names(consensus_cn_gene),)]
# Work in progress
consensus_cn_gene_selected <- as.data.frame(consensus_cn_gene_selected)



