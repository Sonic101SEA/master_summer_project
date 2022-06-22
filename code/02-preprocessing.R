
# Data --------------------------------------------------------------------





# Processing --------------------------------------------------------------
## Processing consensus_cn_gene
### Keeping all columns seen in ids
consensus_cn_gene_selected <- consensus_cn_gene[, grepl(paste(ids_X, collapse = "|"),
                                                        names(consensus_cn_gene),)]
# Work in progress
consensus_cn_gene_selected <- as.data.frame(consensus_cn_gene_selected)



