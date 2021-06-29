###############################################################################
####     5. EXTRACT SUBSEQS AND SUB-INCIDENCE MATRICES FROM EACH CC        ####
###############################################################################
### Define function to get subsequences, origin elements and sub-incidence matrices 
### from multiple-origin connected components

get.subseq.subIncM <- function(cc.origins, incMlist){
  
  print("Generating list of subsequences for multi-origin connected components")
  cc.subseq <- list()
  for(o in 1:length(cc.origins)){
    subseqlistAll <- vector()
    originlist <- cc.origins[[o]]
    for(i in 1:length(incMlist)){
      subincM <- incMlist[[i]][,which(colnames(incMlist[[i]]) %in% originlist)]
      if(is.vector(subincM)){
        subseqlist <- which(subincM!=0)
      }else{
        subseqlist <- which(rowSums(subincM)!=0)
      }
      subseqlistAll <- c(subseqlistAll, names(subseqlist))
    }
    cc.subseq[[i]] <- subseqlistAll
  }
  rm(subseqlist,originlist,subincM,i) # clean memory
  gc()
  
  print("Generating list of sub-incidence matrices for multi-origin connected components")
  cc.subincM <- list()
  for(o in 1:length(cc.origins)){
    subincMAll <- list()
    originlist <- cc.origins[[o]]
    for(i in 1:length(incMlist)){
      subincM <- incMlist[[i]][,which(colnames(incMlist[[i]]) %in% originlist)]
      if(is.vector(subincM)){
        cc.subincM[[i]] <- subincM[which(subincM)!=0]
      }else{
        cc.subincM[[i]] <- subincM[which(rowSums(subincM)!=0),]
      }
      subincMAll[[i]] <- 
    }
  }
  rm(originlist,subincM,i) # clean memory
  gc()
  
  result <- list(cc.subseq=cc.subseq, cc.subincM=cc.subincM)
  return(result)
  
}