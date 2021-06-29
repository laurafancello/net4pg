###############################################################################
####   TEST FOR FUNCTION 3. TO MERGE MULTIPLE ADJACENCY MATRICES            ####
###############################################################################

## Generate input for test
adjM1=cbind(c(1,1,0,0),c(1,0,0,0),c(1,1,0,1),c(1,0,0,1))
colnames(adjM1) <- c("a","b","c","d")
rownames(adjM1) <- colnames(adjM1)
adjM2=cbind(c(0,1,0,0,1),c(1,1,1,1,0),c(1,1,1,0,0),c(1,0,0,1,1),c(1,1,1,0,0))
colnames(adjM2) <- c("a","c","e","f","g")
rownames(adjM2) <- colnames(adjM2)
adjM3=cbind(c(1,1,0),c(1,0,0),c(1,0,1))
colnames(adjM3) <- c("a","d","e")
rownames(adjM3) <- colnames(adjM3)
adjMlist <- list(adjM1=adjM1, adjM2=adjM2, adjM3=adjM3)
adjMlist <- lapply(adjMlist, function(x) apply(x, 2, as.logical)) # convert to logical (to reduce object size)

## Extract complete, non-redundant list of origin elements from all adjacency matrices
allOrigins <- unique(unlist(lapply(adjMlist, function(x) colnames(x))))

## Find origins missing in each adjacency matrix
toAdd <- lapply(adjMlist, function(x) setdiff(allOrigins, colnames(x)))

## Generate columns with missing origins and add them to each adjacency matrix
colsToAdd <- mapply(function(x, y) matrix(nrow = dim(x)[1], ncol = length(y), dimnames = list(rownames(x), y)), adjMlist, toAdd)
colsToAdd <- lapply(colsToAdd, function(x) replace(x, is.na(x), 0)) # replace NAs by 0s
colsToAdd <- lapply(colsToAdd, function(x) apply(x, 2, as.logical)) # convert to logical (to reduce object size)
inputAdj_addedCols <- mapply(function(x,y) cbind(x,y), adjMlist, colsToAdd) # add columns
rm(colsToAdd)
gc()

## Generate rows with origins to add to each adjacency matrix and add them to each adjacency matrix
rowsToAdd <- mapply(function(x, y) matrix(nrow = length(y), ncol = dim(x)[2], dimnames = list(y, colnames(x))), inputAdj_addedCols_logical, toAdd)
rowsToAdd <- lapply(rowsToAdd, function(x) replace(x, is.na(x), 0)) # replace NAs by 0s
rowsToAdd <- lapply(rowsToAdd, function(x) apply(x, 2, as.logical)) # convert to logical (to reduce object size)
inputAdj_addedColsRows <- mapply(function(x,y) rbind(x,y), inputAdj_addedCols, rowsToAdd_ok, SIMPLIFY = F) # add rows
rm(rowsToAdd, inputAdj_addedCols)
gc()

## Order
inputAdj_addedColsRows_ordered <- lapply(inputAdj_addedColsRows, function(x) x[order(as.character(rownames(x))), order(as.character(colnames(x)))])

## Merge into a unique adjacency matrix by sum
adjMfinal <- Reduce("+", inputAdj_addedColsRows_ordered)
diag(adjMfinal) <- rep(0,length(diag(adjMfinal))) # remove self-loops (origins connected with themselves)
adjMfinal[adjMfinal > 1] <- 1


#####################################################################################################
#### TEST FOR FUNCTION 4. TO EXTRACT SUBSEQS AND SUB-INCIDENCE MATRICES FOR SUB-GRAPH PLOTTING   ####
#####################################################################################################

## Generate input for test
# Read incidence matrices
incM1=cbind(c(0,1,0,0,1),c(0,1,0,0,0),c(0,1,0,0,0),c(0,0,0,0,1),c(1,0,1,0,0),c(0,0,1,0,0),c(0,1,0,0,1),c(0,0,0,0,1))
colnames(incM1) <- c("a","b","c","d","e","f","g","h")
rownames(incM1) <- c("pep1","pep2","pep3","pep4","pep5")
incM2=cbind(c(0,0,0,0,1),c(0,0,0,1,0),c(0,0,1,0,0),c(0,0,0,0,1),c(1,0,1,0,0),c(0,0,0,1,0),c(0,0,1,0,0),c(0,1,0,0,0))
colnames(incM2) <- c("a","b","c","d","e","f","g","h")
rownames(incM2) <- c("RNAss1", "RNAss2", "RNAss3", "RNAss4", "RNAss5")
incM3=cbind(c(1,0,0),c(0,1,0),c(1,0,0),c(0,0,1),c(0,0,1),c(1,0,0))
colnames(incM3) <- c("a","b","c","d","e","f")
rownames(incM3) <- c("Episs1", "Episs2", "Episs3")
incMlist <- list(incM1=incM1, incM2=incM2, incM3=incM3)
# Reduce incidence matrices
reduced_incMlist <- lapply(incMlist, reduceincidenceMatrix)
# Generate adjacency matrices
adjM1=t(incM1)%&%incM1
adjM1=as.matrix(adjM1)
diag(adjM1) <- rep(FALSE,length(diag(adjM1)))
adjM2=t(incM2)%&%incM2
adjM2=as.matrix(adjM2)
diag(adjM2) <- rep(FALSE,length(diag(adjM2)))
adjM3=t(incM3)%&%incM3
adjM3=as.matrix(adjM3)
diag(adjM3) <- rep(FALSE,length(diag(adjM3)))
# Merge in unique adjacency matrix
adjMlist <- list(adjM1=adjM1, adjM2=adjM2, adjM3=adjM3)
adjMfinal <- completeAdjMforMerge(adjMlist)
g <- graphAM(adjMfinal, edgemode='undirected', values=NA)
cc.origins <- connComp(as(g, 'graphNEL'))

cc.subseq.subincM <- get.subseq.subincM(cc.origins, incMlist)
  









##### DEPRECATED (APPLY REPLACED BY FOR LOOPS) #####
colsToAdd <- mapply(function(x, y) matrix(nrow = dim(x)[1], ncol = length(y), dimnames = list(rownames(x), y)), adjMlist, toAdd)
colsToAdd <- lapply(colsToAdd, function(x) replace(x, is.na(x), 0)) # replace NAs by 0s
colsToAdd <- lapply(colsToAdd, function(x) apply(x, 2, as.logical)) # convert to logical (to reduce object size)
for(i in 1:length(colsToAdd)){ # put back rownames (lost in logical conversion)
  rownames(colsToAdd[[i]]) <- rownames(adjMlist[[i]])
}
adjMlist <- lapply(adjMlist, function(x) apply(x, 2, as.logical)) # convert to logical (to reduce object size)
inputAdj_addedCols <- mapply(function(x,y) cbind(x,y), adjMlist, colsToAdd) # add columns
rm(colsToAdd)
gc()
# Generate rows with origins to add to each adjacency matrix and add them to each adjacency matrix
rowsToAdd <- mapply(function(x, y) matrix(nrow = length(y), ncol = dim(x)[2], dimnames = list(y, colnames(x))), inputAdj_addedCols_logical, toAdd)
rowsToAdd <- lapply(rowsToAdd, function(x) replace(x, is.na(x), 0)) # replace NAs by 0s
rowsToAdd <- lapply(rowsToAdd, function(x) apply(x, 2, as.logical)) # convert to logical (to reduce object size)
for(i in 1:length(rowsToAdd)){ # put back rownames (lost in logical conversion)
  rownames(rowsToAdd[[i]]) <- toAdd[[i]]
}
inputAdj_addedColsRows <- mapply(function(x,y) rbind(x,y), inputAdj_addedCols, rowsToAdd, SIMPLIFY = F) # add rows
rm(rowsToAdd, inputAdj_addedCols)
gc()
# Order
inputAdj_addedColsRows_ordered <- lapply(inputAdj_addedColsRows, function(x) x[order(as.character(rownames(x))), order(as.character(colnames(x)))])

## Merge into a unique adjacency matrix by sum
adjMfinal <- Reduce("+", inputAdj_addedColsRows_ordered)
diag(adjMfinal) <- rep(0,length(diag(adjMfinal))) # remove self-loops (origins connected with themselves)
adjMfinal[adjMfinal > 1] <- 1 # convert all values >1 to 1

return(adjMfinal)



