
# Data --------------------------------------------------------------------
list_ids_icgc <- scan(here::here("data/ids_of_interest.txt"), character(), quote = "")

## Reading CNV files for ids of interest
all_CNV_files <- list.files(path = "data/consensus_cnv/consensus.20170119.somatic.cna.annotated", pattern = "*.txt")
interest_CNV_files <- subset(all_CNV_files, grepl(paste0(list_ids_icgc, collapse = "|"),
                                                  all_CNV_files))

### For creating multiple objects
# for (i in 1:length(interest_CNV_files)) {
#   assign(interest_CNV_files[i],
#          read.table(paste0(here::here("data/consensus_cnv/consensus.20170119.somatic.cna.annotated/"),
#                       interest_CNV_files[i]), header = TRUE))
# }

### For reading all the CNV files into a list (Use this)
CNV_files_list <- lapply(paste0(here::here("data/consensus_cnv/consensus.20170119.somatic.cna.annotated/"), interest_CNV_files), 
                         read.table, header = TRUE, stringsAsFactors = FALSE)

## Reading SV files for ids of interest

all_SV_files <- list.files(path = "data/consensus_sv/tcga/open", pattern = "*.gz")
interest_SV_files <- subset(all_SV_files, grepl(paste0(list_ids_icgc, collapse = "|"),
                                                all_SV_files))

### For creating multiple objects
# for (i in 1:length(interest_SV_files)) {
#   assign(interest_SV_files[i],
#          read.table(paste0(here::here("data/consensus_sv/tcga/open/"),
#                            interest_SV_files[i]), header = TRUE))
# }

### For reading all the SV files into a list (Use this)
SV_files_list <- lapply(paste0(here::here("data/consensus_sv/tcga/open/"), interest_SV_files), 
                        read.table, header = TRUE, stringsAsFactors = FALSE)


## Reading gene level calls in CNV
### Major
cnv_major_gene_level_calls <- read.table(here::here("data/consensus_cnv/gene_level_calls/all_samples.consensus_CN.by_gene.170214.txt"), 
                                         fill = TRUE, header = TRUE)

## Reading pan-cancer compendium CN signatures
pan_cancer_cn_signatures <- read.csv(here::here("data/chromosomal_instability_data/PCAWG_activities_scaled.csv"))

### Subset selected ids for pan-cancer compendium CN signatures
pan_cancer_cn_signatures_selected <- subset(pan_cancer_cn_signatures, grepl(paste0(list_ids_icgc, collapse = "|"),
                                                                            pan_cancer_cn_signatures$X))


# Functions for other processing methods ---------------------------------------------------------------

## To drop columns
drop_columns_CNV <- function(CNV_dataframe)
{
  CNV_dataframe <- CNV_dataframe[, -c(4:38)]
}

## To generate segment value from start and end of chromosome
generate_segValue <- function(CNV_dataframe)
{
  CNV_dataframe$segVal <- CNV_dataframe$end - CNV_dataframe$start
  CNV_dataframe
}

## To generate total copy number for each sample, including major and minor
generate_total_cn_mut_load <- function(CNV_dataframe)
{
  total_cn <- sum(CNV_dataframe$total_cn, na.rm = TRUE)
  major_cn <- sum(CNV_dataframe$major_cn, na.rm = TRUE)
  minor_cn <- sum(CNV_dataframe$minor_cn, na.rm = TRUE)
  out <- list("total_cn" = total_cn, "major_cn" = major_cn, "minor_cn" = minor_cn)
  return(out)
}


# Copy Number Functions ---------------------------------------------------------------

## Helper Functions
source(here::here("../cnsignatures/helper_functions.R"))

## Main Functions
## extractCopynumnberFeatures: This function takes as input a collection of absolute copy-number profiles and returns a list of copy-number features extracted from these samples. 
## Copy-number profiles can be input as either a QDNAseq object, or as a list of segment tables. 
##The segment tables (one for each sample) should have the following column headers: "chromosome", "start", "end", "segVal".  

extractCopynumberFeatures<-function(CN_data, cores = 1)
{
  #get chromosome lengths
  chrlen<-read.table(here::here("../cnsignatures/data/hg19.chrom.sizes.txt"), sep="\t",stringsAsFactors = F)[1:24,]
  
  #get centromere locations
  gaps<-read.table(here::here("../cnsignatures/data/gap_hg19.txt"),sep="\t",header=F,stringsAsFactors = F)
  centromeres<-gaps[gaps[,8]=="centromere",]
  
  if(cores > 1) {
    require(foreach)
    doMC::registerDoMC(cores)
    
    temp_list = foreach::foreach(i=1:6) %dopar% {
      if(i == 1){
        list(segsize = getSegsize(CN_data) )
      } else if (i == 2) {
        list(bp10MB = getBPnum(CN_data,chrlen) )
      } else if (i == 3) {
        list(osCN = getOscilation(CN_data,chrlen) )
      } else if (i == 4) {
        list(bpchrarm = getCentromereDistCounts(CN_data,centromeres,chrlen) )
      } else if (i == 5) {
        list(changepoint = getChangepointCN(CN_data) )
      } else {
        list(copynumber = getCN(CN_data) )
      }
      
    }
    unlist( temp_list, recursive = FALSE )
  } else {  
    
    segsize<-getSegsize(CN_data)
    bp10MB<-getBPnum(CN_data,chrlen)
    osCN<-getOscilation(CN_data,chrlen)
    bpchrarm<-getCentromereDistCounts(CN_data,centromeres,chrlen)
    changepoint<-getChangepointCN(CN_data)
    copynumber<-getCN(CN_data)
    
    list(segsize=segsize,bp10MB=bp10MB,osCN=osCN,bpchrarm=bpchrarm,changepoint=changepoint,copynumber=copynumber)
  }
  
}

## generateSampleByComponentMatrix: Given a set of extracted copy number features profiles this function returns a sum-of-posterior sample-by-component matrix. 
## If the all_components argument is specified, then the sum-of-posteriors is calculated using these components, otherwise the component definitions from the manuscript are used.

generateSampleByComponentMatrix<-function(CN_features, all_components=NULL, cores = 1, rowIter = 1000, subcores = 2)
{
  if(is.null(all_components))
  {
    all_components<-readRDS(here::here("../cnsignatures/data/component_parameters.rds"))
  }
  
  if(cores > 1){
    require(foreach)
    
    feats = c( "segsize", "bp10MB", "osCN", "changepoint", "copynumber", "bpchrarm" )
    doMC::registerDoMC(cores)
    
    full_mat = foreach(feat=feats, .combine=cbind) %dopar% {
      calculateSumOfPosteriors(CN_features[[feat]],all_components[[feat]], 
                               feat, rowIter = rowIter, cores = subcores)
    }
  } else {
    full_mat<-cbind(
      calculateSumOfPosteriors(CN_features[["segsize"]],all_components[["segsize"]],"segsize"),
      calculateSumOfPosteriors(CN_features[["bp10MB"]],all_components[["bp10MB"]],"bp10MB"),
      calculateSumOfPosteriors(CN_features[["osCN"]],all_components[["osCN"]],"osCN"),
      calculateSumOfPosteriors(CN_features[["changepoint"]],all_components[["changepoint"]],"changepoint"),
      calculateSumOfPosteriors(CN_features[["copynumber"]],all_components[["copynumber"]],"copynumber"),
      calculateSumOfPosteriors(CN_features[["bpchrarm"]],all_components[["bpchrarm"]],"bpchrarm"))
  }
  
  rownames(full_mat)<-unique(CN_features[["segsize"]][,1])
  full_mat[is.na(full_mat)]<-0
  full_mat
}

## quantifySignautres: Given a sample-by-component matrix this function quantifies signature exposures using the LCD function from the YAPSA package, returning a normalised signature-by-sample matrix. 
## If the component_by_signature matrix is specified then this matrix is used to define the signatures otherwise the signature definitions from the manuscript are used.

quantifySignatures<-function(sample_by_component,component_by_signature=NULL)
{
  if(is.null(component_by_signature))
  {
    component_by_signature<-readRDS(here::here("data/feat_sig_mat.rds"))
  }
  signature_by_sample<-YAPSA::LCD(t(sample_by_component),
                                  YAPSA:::normalize_df_per_dim(component_by_signature,2))
  signature_by_sample<-normaliseMatrix(signature_by_sample)
  signature_by_sample
}


# SV Classes --------------------------------------------------------------

## To obtain all sv classes in dataset
unique(unlist(lapply(SV_files_list, `[[`, "svclass")))

# CNV processing --------------------------------------------------------------
## Processing CNV data
### Dropping columns
drop_CNV_files <- lapply(CNV_files_list, drop_columns_CNV)

### Generating segment values
segVal_CNV_files <- lapply(drop_CNV_files, generate_segValue)

### Adding up copy numbers for each sample
cn_mut_load <- lapply(CNV_files_list, generate_total_cn_mut_load)


# Gene Level Calls --------------------------------------------------------

## Processing gene level calls for CNV
### Keeping all columns seen in ids
consensus_cn_gene_selected <- consensus_cn_gene[, grepl(paste(ids_X, collapse = "|"),
                                                        names(consensus_cn_gene),)]
# Work in progress
consensus_cn_gene_selected <- as.data.frame(consensus_cn_gene_selected)

## Copy Number Signatures
extractCopynumberFeatures(segVal_CNV_files)

