###############################################################################
####           3. MERGE MULTIPLE ADJACENCY MATRICES                      ####
###############################################################################
# Merge multiple adjacency matrices and generate a unique adjacency matrix describing
# if each pair of origin elements is sharing (value 1) or not (value 0) at least one
# subsequence. 
# corresponding to transcripts observed in the other omic and missing
# Step 2. order columns and rows in the same way for RNA and PROT adjacency matrices
# Step 3. sum the two matrices
# Step 4. convert all values >1 to exactly 1 in the resulting matrix so that it becomes a binary matrix of 0s and 1s

completeAdjMforMerge <- function(adjMlist){
  
  ## Make all adjacency matrices of the same size, adding to each one separately missing rows and columns
  
  # Extract complete, non-redundant list of origin elements from all adjacency matrices
  allOrigins <- unique(unlist(lapply(adjMlist, function(x) colnames(x))))
  # Find origins missing in each adjacency matrix
  toAdd <- lapply(adjMlist, function(x) setdiff(allOrigins, colnames(x)))
  inputAdj_addedColsRows_ordered <- list()
  for(i in 1:length(toAdd)){
    if(length(toAdd[[i]])>0){
      # Generate columns with missing origins and add them to each adjacency matrix
      colsToAdd <- matrix(nrow = dim(adjMlist[[i]])[1], ncol = length(toAdd[[i]]), dimnames = list(rownames(adjMlist[[i]]), toAdd[[i]]))
      colsToAdd <- replace(colsToAdd, is.na(colsToAdd), 0) # replace NAs by 0s
      colsToAdd <- apply(colsToAdd, 2, as.logical) # convert to logical (to reduce object size)
      rownames(colsToAdd) <- rownames(adjMlist[[i]]) # put back rownames (lost in logical conversion)
      rownames_adjM <- rownames(adjMlist[[i]])
      adjMlist[[i]] <- apply(adjMlist[[i]], 2, as.logical) # convert to logical (to reduce object size)
      rownames(adjMlist[[i]]) <- rownames_adjM # put back rownames (lost in logical conversion)
      inputAdj_addedCols <- cbind(adjMlist[[i]], colsToAdd) # add columns
      rm(colsToAdd)
      gc()
      # Generate rows with origins to add to each adjacency matrix and add them to each adjacency matrix
      rowsToAdd <- matrix(nrow = length(toAdd[[i]]), ncol = dim(adjMlist[[i]])[2]+length(toAdd[[i]]), dimnames = list(toAdd[[i]], colnames(inputAdj_addedCols)))
      rowsToAdd <- replace(rowsToAdd, is.na(rowsToAdd), 0) # replace NAs by 0s
      rowsToAdd <- apply(rowsToAdd, 2, as.logical) # convert to logical (to reduce object size)
      rownames(rowsToAdd) <- toAdd[[i]] # put back rownames (lost in logical conversion)
      inputAdj_addedColsRows <- rbind(inputAdj_addedCols, rowsToAdd) # add rows
      rm(rowsToAdd, inputAdj_addedCols)
      gc()
    }else{
      inputAdj_addedColsRows <- adjMlist[[i]]
    }
    # Order
    inputAdj_addedColsRows_ordered[[i]] <- inputAdj_addedColsRows[order(as.character(rownames(inputAdj_addedColsRows))), order(as.character(colnames(inputAdj_addedColsRows)))]
  }
  ## Merge into a unique adjacency matrix by sum
  adjMfinal <- Reduce("+", inputAdj_addedColsRows_ordered)
  diag(adjMfinal) <- rep(0,length(diag(adjMfinal))) # remove self-loops (origins connected with themselves)
  adjMfinal[adjMfinal > 1] <- 1 # convert all values >1 to 1
  
  return(adjMfinal)
}