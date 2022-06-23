
# Data --------------------------------------------------------------------





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

# Processing --------------------------------------------------------------

## Copy Number Signatures
extractCopynumberFeatures(segVal_CNV_files)

