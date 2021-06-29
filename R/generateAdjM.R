###############################################################################
###                2. GENERATE ADJACENCY MATRIX                            ####
###############################################################################

# Define function to generate adjacency matrix of origin-to-origin connections,
# based on multi-mapping subsequences) by cross product of the incidence matrix
# describing origin-to-subsequence assignments

generateAdjM <- function(incM){
  
  # Transform dataframe object into logical matrix object (to reduce its size in memory)
  SplitSizeCols <- round(dim(incM)[2]/5)
  incM1 <- as.matrix(incM[,1:SplitSizeCols])
  incM2 <- as.matrix(incM[,(SplitSizeCols+1):(SplitSizeCols*2)])
  incM3 <- as.matrix(incM[,((SplitSizeCols*2)+1):(SplitSizeCols*3)])
  incM4 <- as.matrix(incM[,((SplitSizeCols*3)+1):(SplitSizeCols*4)])
  incM5 <- as.matrix(incM[,((SplitSizeCols*4)+1):dim(incM)[2]])
  incM1 <- apply(incM1, 2, as.logical)
  incM2 <- apply(incM2, 2, as.logical)
  incM3 <- apply(incM3, 2, as.logical)
  incM4 <- apply(incM4, 2, as.logical)
  incM5 <- apply(incM5, 2, as.logical)
  rm(incM)
  gc()
  incM <- cbind(incM1, incM2, incM3, incM4, incM5)
  rm(incM1, incM2, incM3, incM4, incM5, SplitSizeCols) # clean memory
  gc()
  
  # Convert to adjacency matrix
  adjM=t(incM)%&%incM
  print(paste0("Dimensions adjacency matrix ", paste(dim(adjM)[1], ", ", dim(adjM)[2]), "\n"))
  #rm(incM)
  #gc()
  
  # Transform ngCMatrix object into matrix object
  SplitSizeCols <- round(dim(adjM)[2]/5)
  adjM1 <- as.matrix(adjM[,1:SplitSizeCols])
  adjM2 <- as.matrix(adjM[,(SplitSizeCols+1):(SplitSizeCols*2)])
  adjM3 <- as.matrix(adjM[,((SplitSizeCols*2)+1):(SplitSizeCols*3)])
  adjM4 <- as.matrix(adjM[,((SplitSizeCols*3)+1):(SplitSizeCols*4)])
  adjM5 <- as.matrix(adjM[,((SplitSizeCols*4)+1):dim(adjM)[2]])
  #rm(adjM)
  #gc()
  adjM <- cbind(adjM1, adjM2, adjM3, adjM4, adjM5)
  rm(incM1, incM2, incM3, incM4, incM5, SplitSizeCols) # clean memory
  gc()
  
  # Remove self-connecting edges
  diag(adjM) <- rep(FALSE,length(diag(adjM))) 
  
  return(adjM)
}