###############################################################################
###               1. REDUCE INCIDENCE MATRICES                             ####
###############################################################################
# Function to reduce size of incidence matrix by removing unique (not multi-mapping)
# subsequences and origin elements with only unique subsequences

reduceIncM <- function(incM){
  
  ### Remove unique (not multi-mapping) subsequences
  incM_RowFilter <- incM[-which(rowSums(incM)==1),]
  dim_RowFilter <- dim(incM_RowFilter)
  
  ### Remove origin elements with only 0 values: they only have connections from
  ### unique subsequences which were removed in previous step
  OriginsWithOnlyUniqueSubseqs <- which(colSums(incM_RowFilter)==0)
  if(length(OriginsWithOnlyUniqueSubseqs)!=0){
    incM_RowColFilter <- incM_RowFilter[,-OriginsWithOnlyUniqueSubseqs] 
  }else{
    incM_RowColFilter <- incM_RowFilter
  }
  
  result <- incM_RowColFilter
  return(result)
}