####### NOTE: CODE MODIFIED FROM "Generate_read_peptide_transcript_graph_PSM_and_RNA_22042020_ImprovedCode.R" ######

library("data.table") # Read big input files containing incidence matrices with freqd function
library("graph") # Generate graph from adjacency matrix and calculate connected components
library("dplyr") # Merge RNA and PROT incidence matrices by bind_rows() into a unique incidence matrix for connected components subgraphs visualization
library("igraph") # Plot connecetd component subgraphs
library("ggplot2") # Generate scatterplots section 11 ("further analyses")
library("pdp") # for grid.arrange() in plots section 11 ("further analyses")
library("stringr") # for str_replace (to conevrt simplified ENST transcript ids into ENST00000 and viceversa)
library("tidyr") # for replace_na() into completeAdjMforMerge() in-house function
library("Matrix") # for %&% function in adjacency matrix generation
library("rlist") # to filter list of max scores in compariosn between SingleCC_NoRNAsubseq and SingleCC_RNAanadPeptides

path=getwd()
source(paste0(path, "/Functions_GraphAndCCs.R"))


# ########################################################################################################
# ####                    1. READ RNA INCIDENCE MATRIX                                                ####
# ########################################################################################################
# # Issue to upload entire input matrix in RStudio by graphical user interface (scp transfer not easily available 
# # on IFB RStudio virtual machine. Therefore, input matrix (subseq \t transcriptID) was split by bash into submatrices
# # and submatrices were uploaded and read one by one and concatenated
# 
# #Read header (colnames)
# transcriptIDs <- fread(file=paste0(path, "/Header_CrossTable_RNAsubsequences"), sep=" ", header=F, keepLeadingZeros = TRUE, data.table=FALSE)
# transcriptIDs <- unlist(transcriptIDs[1,])
# 
# for (i in 1:10){
#   #Read rownames
#   subseqIDs <- fread(file=paste0(path, "/CrossTable_RNAsubsequences_NoHeader", i), sep=" ", header=F, data.table=F, select=1)
#   subseqIDs <- subseqIDs$V1
#   #Read data.frame
#   incidenceMatrix <- fread(file=paste0(path, "/CrossTable_RNAsubsequences_NoHeader", i), sep=" ", header=F, data.table=F, drop=1)
#   #Set col and row names
#   colnames(incidenceMatrix) <- as.character(as.vector(transcriptIDs))
#   rownames(incidenceMatrix) <- as.character(as.vector(subseqIDs))
#   assign(paste("incidenceMatrix",i, sep=""), (incidenceMatrix))
#   }
# 
# #Concatenate and save in RDS object
# incidenceMatrixRNA <- rbind(incidenceMatrix1, incidenceMatrix2, incidenceMatrix3, incidenceMatrix4,incidenceMatrix5,incidenceMatrix6,incidenceMatrix7,incidenceMatrix8,incidenceMatrix9,incidenceMatrix10)
# saveRDS(incidenceMatrixRNA, file=paste0(path, "/incidenceMatrix_RNA.RDS" ), compress="gzip")
# rm(incidenceMatrix1, incidenceMatrix2, incidenceMatrix3, incidenceMatrix4,incidenceMatrix5,incidenceMatrix6,incidenceMatrix7,incidenceMatrix8,incidenceMatrix9,incidenceMatrix10, subseqIDs1, subseqIDs2, subseqIDs3, subseqIDs4, subseqIDs5, subseqIDs6, subseqIDs7, subseqIDs8, subseqIDs9, subseqIDs10)
# gc()
# 
# RowNamesRNASimplified <- cbind(as.character(as.vector(paste0("R_",1:dim(incidenceMatrixRNA)[1]))), rownames(incidenceMatrixRNA))
# colnames(RowNamesRNASimplified)=c("subsequence_id", "subsequence_coord")
# rownames(incidenceMatrixRNA) <- RowNamesRNASimplified[,1]
# saveRDS(incidenceMatrixRNA, file=paste0(path, "/incidenceMatrix_RNA_Rownamessimplified.RDS" ), compress="gzip")
# 
# write.table(file=paste0(path,"/RowNamesRNASimplified.txt"), RowNamesRNASimplified, col.names = T, row.names = F, sep="\t", quote=F)


# ########################################################################################################
# ####                    2. READ PROTEOMICS INCIDENCE MATRIX                                         ####
# ########################################################################################################
# # Issue to upload entire input matrix in RStudio by graphical user interface (scp transfer not easily available 
# # on IFB RStudio virtual machine. Therefore, input matrix (peptide \t transcriptID) was split by bash into submatrices
# # and submatrices were uploaded and read one by one and concatenated
# 
# # Read header (colnames)
# transcriptIDs <- fread(file=paste0(path, "/Header_CrossTable_Proteomics_ProlineIDs"), sep=" ", header=F, keepLeadingZeros = TRUE, data.table=FALSE)
# transcriptIDs <- unlist(transcriptIDs[1,])
# 
# #Read rownames
# subseqIDs <- fread(file=paste0(path, "/Rownames_CrossTable_Proteomics_ProlineIDs"), sep=" ", header=F, data.table=F, select=1)
# subseqIDs <- subseqIDs$V1
# 
# #Read data.frame
# for (i in 1:10){
#   incidenceMatrix <- fread(file=paste0(path, "/CrossTable_Proteomics_ProlineIDs_", i), sep=" ", header=F, data.table=F)
#   assign(paste("incidenceMatrix",i, sep=""), (incidenceMatrix))
# }
# incidenceMatrix_PSM <- rbind(incidenceMatrix1, incidenceMatrix2, incidenceMatrix3, incidenceMatrix4, incidenceMatrix5, incidenceMatrix6, incidenceMatrix7, incidenceMatrix8, incidenceMatrix9, incidenceMatrix10)
# 
# #Set col and row names
# colnames(incidenceMatrix_PSM)=as.character(as.vector(transcriptIDs))
# rownames(incidenceMatrix_PSM)=as.character(as.vector(subseqIDs))
# colnames(incidenceMatrix_PSM) <- str_replace(colnames(incidenceMatrix_PSM), "ENST00000", "ENST")
# dim(incidenceMatrix_PSM) # 97351 53031
# rm(incidenceMatrix1, incidenceMatrix2, incidenceMatrix3, incidenceMatrix4, incidenceMatrix5, incidenceMatrix6, incidenceMatrix7, incidenceMatrix8, incidenceMatrix9, incidenceMatrix10)  
# saveRDS(incidenceMatrix_PSM, file=paste0(path, "/incidenceMatrix_PSM.RDS" ), compress="gzip")
# 
# RowNamesPROTSimplified <- cbind(as.character(as.vector(paste0("P_",1:dim(incidenceMatrix_PSM)[1]))), rownames(incidenceMatrix_PSM))
# colnames(RowNamesPROTSimplified)=c("peptide_id", "peptide_seq")
# rownames(incidenceMatrix_PSM) <- RowNamesPROTSimplified[,1]
# saveRDS(incidenceMatrix_PSM, file=paste0(path, "/incidenceMatrix_PSM_Rownamessimplified.RDS" ), compress="gzip")
# write.table(file=paste0(path,"/RowNamesPROTSimplified.txt"), RowNamesPROTSimplified, col.names = T, row.names = F, sep="\t", quote=F)

incMlist <- list(incM1=incM1, incM2=incM2, incM3=incM3)


#######################################################################################################
####                        3. REDUCE INCIDENCE MATRICES                                           ####
#######################################################################################################
# Reduce size of incidence matrice by removing unique (not shared) RNA subsequences/peptides and
# transcripts with only unique RNA subsequences/peptides and no shared ones

##### 3.1 DEFINE FUNCTION TO REDUCE INCIDENCE MATRICES
reduceIncidenceMatrix <- function(incidenceMatrix, RDSfile){
  ### Remove unique (not shared) RNA subsequences/peptides 
  incidenceMatrix_RowFilter <- incidenceMatrix[-which(rowSums(incidenceMatrix)==1),]
  dim_RowFilter <- dim(incidenceMatrix_RowFilter)
  rm(incidenceMatrix)
  
  ### Remove transcripts with only 0 values: they only have unique RNA subsequences/peptides which were removed in previous step
  incidenceMatrix_RowColFilter <- incidenceMatrix_RowFilter[,-which(colSums(incidenceMatrix_RowFilter)==0)] 
  dim_RowColFilter <- dim(incidenceMatrix_RowColFilter)
  #Save in RDS object and clean
  saveRDS(incidenceMatrix_RowColFilter, file=paste0(path, "/", RDSfile, ".RDS"), compress="xz")
  result <- list(dim_RowFilter=dim_RowFilter, dim_RowColFilter=dim_RowColFilter, incidenceMatrix_RowColFilter=incidenceMatrix_RowColFilter)
  return(result)
}

##### 3.2 APPLY FUNCTION TO REDUCE RNA and PROT INCIDENCE MATRICES
filterRNA <- reduceIncidenceMatrix(incidenceMatrixRNA, "incidenceMatrixRNA_RowColFilter")
dim_RowFilterRNA <- filterRNA[[1]] # 28661 57247
dim_RowColFilterRNA <- filterRNA[[2]] # 28661 46057
incidenceMatrixRNA_RowColFilter <- filterRNA[[3]]
saveRDS(file=paste0(path, "/incidenceMatrixRNA_RowColFilter.RDS"), incidenceMatrixRNA_RowColFilter)
filterPROT <- reduceIncidenceMatrix(incidenceMatrixPSM, "incidenceMatrixPROT_RowColFilter")
dim_RowFilterPROT <- filterPROT[[1]] # 88893 53031
dim_RowColFilterPROT <- filterPROT[[2]] # 88893 51677
incidenceMatrixPROT_RowColFilter <- filterPROT[[3]]
saveRDS(file=paste0(path, "/incidenceMatrixPROT_RowColFilter.RDS"), incidenceMatrixPROT_RowColFilter)


########################################################################################################
###                                          4. GENERATE ADJACENCY MATRICES                         ####
########################################################################################################

##### 4.1 DEFINE FUNCTION TO GENERATE RNA and PROT ADJACENCY MATRICES FROM RNA and PROT INCDENCE MATRICES
generateAdjMatrix <- function(incM, adjacencyFilename){
  # Read incidence matrix
  print(paste0("Dimensions incidence matrix ", dim(incM), "\n"))
  print(paste0("is(incidence matrix) ", is(incM), "\n"))
  print(paste0("object.size(incidence matrix) ", object.size(incM), "\n"))
  
  # Convert to adjacency matrix
  adjM=t(incM)%&%incM
  print(paste0("Dimensions adjacency matrix ", paste(dim(adjM)[1], ", ", dim(adjM)[2]), "\n"))
  print(paste0("object.size(adjacency matrix) ", object.size(adjM), "\n"))
  rm(incM)
  gc()
  # Transform ngCMatrix object into dataframe object (required to merge PROT - RNA)
  SplitSizeCols <- round(dim(adjM)[2]/5)
  adjM1 <- as.matrix(adjM[,1:SplitSizeCols])
  adjM2 <- as.matrix(adjM[,(SplitSizeCols+1):(SplitSizeCols*2)])
  adjM3 <- as.matrix(adjM[,((SplitSizeCols*2)+1):(SplitSizeCols*3)])
  adjM4 <- as.matrix(adjM[,((SplitSizeCols*3)+1):(SplitSizeCols*4)])
  adjM5 <- as.matrix(adjM[,((SplitSizeCols*4)+1):dim(adjM)[2]])
  rm(adjM)
  gc()
  adjM1=as.data.frame(adjM1)
  adjM2=as.data.frame(adjM2)
  adjM3=as.data.frame(adjM3)
  adjM4=as.data.frame(adjM4)
  adjM5=as.data.frame(adjM5)
  adjM <- cbind(adjM1, adjM2, adjM3, adjM4, adjM5)
  print(paste0("Dimensions adjacency matrix converted to dataframe", dim(adjM), "\n")) 
  object.size(adjM)
  return(adjM)
}

##### 4.2 APPLY FUNCTION TO GENERATE ADJACENCY MATRICES
adjRNA <- generateAdjMatrix(incidenceMatrixRNA_RowColFilter, "adjRNA") 
# Incidence matrix
# dimensions: 28661, 46057
# is: matrix array mMatrix structure vector
# object.size: 5,285,309,608 bytes
# Adjacency matrix
# dimensions. 46057 46057
# is: "data.frame" "list"       "oldClass"   "vector"    
# object.size: 8562376 bytes
saveRDS(file=paste0(path,"/adjRNA.RDS"), adjRNA, compress="xz")
rm(incidenceMatrixRNA_RowColFilter)
adjPROT <- generateAdjMatrix(incidenceMatrixPROT_RowColFilter, "adjPROT")
# Incidence matrix
# dimensions: 88893 51677
# is: matrix array mMatrix structure vector
# object.size: 18387405360 bytes
# Adjacency matrix
# dimensions. 51677 51677
# is: "data.frame" "list"       "oldClass"   "vector"    
# object.size: 10,961,336 bytes
saveRDS(file=paste0(path,"/adjPROT.RDS"), adjPROT, compress="xz")
rm(incidenceMatrixPROT_RowColFilter)

##### 4.3 REPLACE ALL VALUES ON DIAGONAL BY "FALSE", BECASUE THEY ONLY REPRESENT THE CONNECTION OF A TRANSCRIPT WITH ITSELF (TRANSCRIPT i VS TRANSCRIPT i)
adjRNA <-as.matrix(adjRNA)
adjPROT <-as.matrix(adjPROT)
diag(adjRNA) <- rep(FALSE,length(diag(adjRNA))) # remove self-connecting edges
diag(adjPROT) <- rep(FALSE,length(diag(adjPROT))) # remove self-connecting edges
saveRDS(adjRNA, file=paste0(path, "/adjRNA_noSelfLoops.RDS"), compress="xz")
saveRDS(adjPROT, file=paste0(path, "/adjPROT_noSelfLoops.RDS"), compress="xz")


########################################################################################################
###                                  5. MERGE RNA - PROT ADJACENCY MATRICES                         ####
########################################################################################################
# Generate a unique adjacency matrix describing if each pair of transcripts is sharing (value 1)
# or not (value 0) at least one RNA subsequence or peptide, following these steps:
# 1. make RNA and PROT adjacency matrices of the same size, adding to each separately rows and columns
# corresponding to transcripts observed in the other omic and missing
# 2. order columns and rows in the same way for RNA and PROT adjacency matrices
# 3. sum the two matrices
# 4. convert all values >1 to exactly 1 in the resulting matrix so that it becomes a binary matrix of 0s and 1s

##### 5.1 DEFINE FUNCTION TO COMPLETE EACH OMIC (RNA AND PROT) ADJACENCY MATRIX WITH TRANSCRIPTS IDENTIFIED ONLY IN THE OTHER OMIC
completeAdjMforMerge <- function(inputAdj, otherAdj){
  ## Find transcripts which are in otherAdj but not in inputAdj adjacency matrix
  transcrToAdd <- setdiff(colnames(otherAdj),colnames(inputAdj))
  print(paste0("Length transcrToAdd ", length(transcrToAdd), "\n"))
  
  ## Generate columns with transcripts to add to adjacency matrix
  colsToAdd <- matrix(nrow = dim(inputAdj)[1], ncol = length(transcrToAdd))
  colnames(colsToAdd) <- transcrToAdd
  rownames(colsToAdd) <- c(rownames(inputAdj))
  colsToAdd[(colsToAdd)] <-0
  ## Add columns to inutAdj
  inputAdj_addedCols <- cbind(inputAdj, colsToAdd)
  rm(colsToAdd)
  gc()
  ## Transform ngcMatrix object in logical matrix (subsetting for memory issues)
  inputAdj_addedCols <- as.matrix(inputAdj_addedCols)
  inputAdj_addedCols_logical <- apply(inputAdj_addedCols, 2, as.logical)
  rm(inputAdj_addedCols)
  gc()
  
  ## Generate rows with transcripts to add to RNA adjacency matrix
  rowsToAdd <- matrix(nrow = length(transcrToAdd), ncol = dim(inputAdj)[2]+length(transcrToAdd))
  colnames(rowsToAdd) <- c(colnames(inputAdj), transcrToAdd)
  rownames(rowsToAdd) <- transcrToAdd
  rowsToAdd[is.na(rowsToAdd)] <- 0
  rowsToAdd_logical <- apply(rowsToAdd, 2, as.logical)
  rm(rowsToAdd)
  gc()
  inputAdj_addedColsRows <- rbind(inputAdj_addedCols_logical, rowsToAdd_logical)
  rownames(inputAdj_addedColsRows) <- colnames(inputAdj_addedColsRows)
  rm(rowsToAdd_logical,inputAdj_addedCols_logical)
  gc()
  return(inputAdj_addedColsRows)
}

##### 5.2 APPLY FUNCTION TO COMPLETE EACH OMIC (RNA AND PROT) ADJACENCY MATRIX WITH TRANSCRIPTS IDENTIFIED ONLY IN THE OTHER OMIC
adjRNA_addedColsRows <- completeAdjMforMerge(adjRNA,adjPROT)
saveRDS(file=paste0(path,"/adjRNA_addedColsRows.RDS"), adjRNA_addedColsRows, compress="xz")
#object.size: 24,238,753,464 bytes
#dim: 77835 77835

adjPROT_addedColsRows <- completeAdjMforMerge(adjPROT,adjRNA)
saveRDS(file=paste0(path,"/adjPROT_addedColsRows.RDS"), adjPROT_addedColsRows, compress="xz")
#object.size: 24,238,753,464 bytes
#dim: 77835 77835
# matrix, logical

##### 5.3 ORDER IN THE SAME WAY ROWS/COLUMNS OF COMPLETE RNA AND PROT ADJACENCY MATRICES
adjRNA_addedColsRows_ordered <- adjRNA_addedColsRows[order(rownames(adjRNA_addedColsRows)),order(colnames(adjRNA_addedColsRows))]
adjPROT_addedColsRows_ordered <- adjPROT_addedColsRows[order(rownames(adjPROT_addedColsRows)),order(colnames(adjPROT_addedColsRows))]
rm(adjRNA_addedColsRows, adjPROT_addedColsRows)
gc()
saveRDS(file=paste0(path, "/adjPROT_addedColsRows_ordered.RDS"), adjPROT_addedColsRows_ordered)
saveRDS(file=paste0(path, "/adjRNA_addedColsRows_ordered.RDS"), adjRNA_addedColsRows_ordered)

##### 5.4 PUT TOGETHER RNA AND PROT ADJACENCY MATRICES
adjPROTandRNA=adjPROT_addedColsRows_ordered+adjRNA_addedColsRows_ordered
rm(adjPROT_addedColsRows_ordered,adjRNA_addedColsRows_ordered)
gc()
dim(adjPROTandRNA) # 77835 77835
object.size(adjPROTandRNA) # 24,244,357,632 bytes
is(adjPROTandRNA) # "matrix"    "array"     "mMatrix"   "structure" "vector"   
is(adjPROTandRNA[,1]) # "integer" "double" "numeric" "vector" "data.frameRowLabels" "atomicVector"  "numericVector" "index" "replValue" "numLike" "number" "Mnumeric" "replValueSp"
saveRDS(file=paste0(path,"/adjPROTandRNA.RDS"), adjPROTandRNA, compress="xz")
#Remove self-loops (transcripts connected with themselves)
diag(adjPROTandRNA) <- rep(0,length(diag(adjPROTandRNA)))
dim(adjPROTandRNA) # 77835 77835
object.size(adjPROTandRNA) # 48,477,506,528 bytes
is(adjPROTandRNA) # matrix"    "array"     "mMatrix"   "structure" "vector" 
is(adjPROTandRNA[,1]) # "numeric" "vector"
gc()
saveRDS(file=paste0(path,"/adjPROTandRNA_NoSelfLoops.RDS"), adjPROTandRNA, compress="xz")

##### 5.5 TRANSFORM IN BINARY MATRIX: VALUES >1 REPLACED BY 1
SplitSizeCols <- round(dim(adjPROTandRNA)[2]/5)
adjPROTandRNA1 <- as.matrix(adjPROTandRNA[,1:SplitSizeCols])
adjPROTandRNA1[adjPROTandRNA1 > 1] <- 1
adjPROTandRNA2 <- as.matrix(adjPROTandRNA[,(SplitSizeCols+1):(SplitSizeCols*2)])
adjPROTandRNA2[adjPROTandRNA2 > 1] <- 1
adjPROTandRNA3 <- as.matrix(adjPROTandRNA[,((SplitSizeCols*2)+1):(SplitSizeCols*3)])
adjPROTandRNA3[adjPROTandRNA3 > 1] <- 1
adjPROTandRNA4 <- as.matrix(adjPROTandRNA[,((SplitSizeCols*3)+1):(SplitSizeCols*4)])
adjPROTandRNA4[adjPROTandRNA4 > 1] <- 1
adjPROTandRNA5 <- as.matrix(adjPROTandRNA[,((SplitSizeCols*4)+1):dim(adjPROTandRNA)[2]])
adjPROTandRNA5[adjPROTandRNA5 > 1] <- 1
adjPROTandRNA_Binary <- cbind(adjPROTandRNA1,adjPROTandRNA2,adjPROTandRNA3,adjPROTandRNA4,adjPROTandRNA5)
dim(adjPROTandRNA_Binary) # 77835 77835
object.size(adjPROTandRNA_Binary) # 48,477,506,528 bytes
saveRDS(file=paste0(path,"/adjPROTandRNA_Binary.RDS"), adjPROTandRNA_Binary, compress="xz")


######################################################################################################################
####                  6. CALCULATE CONNECTED COMPONENTS CONTAINING >1 TRANSCRIPT (SHARED SUBSEQUENCES)          ####
######################################################################################################################

##### 6.1 GENERATE GRAPH
print("Building graph object")
ptm <- proc.time() # Start the clock
g <- graph::graphAM(adjPROTandRNA_Binary, edgemode='undirected', values=NA)
proc.time() - ptm # Stop the clock
# user  system elapsed 
#243.624  83.369 326.999 
saveRDS(file=paste0(path, "/g.RDS"), g)

##### DEPRECATED: DEGREES OF GRAPH.
# Since this graph is only built upon transcripts with at least one shared RNA_subseq/peptide, its degrees are not that informative: all nodes (transcripts) with unique
# RNA_subseqs/peptides are missing
# graphDegrees <- graph::degree(g)
# write.table(file=paste0(path, "/graphDegrees.txt"), graphDegrees, sep="\t", row.names=F, col.names=F, quote=F)
# tableDegreeNnodes <- as.data.frame(table(graphDegrees))
# colnames(tableDegreeNnodes) <- c("degree", "N_nodes")
# write.table(file=paste0(path, "/tableDegreeNnodes.txt"), tableDegreeNnodes, sep="\t", row.names=F, col.names=F, quote=F)
# pdf(paste0(path, "/Degree_Nnodes.pdf"))
# plot(tableDegreeNnodes)
# dev.off()
# pdf(paste0(path, "/Degree_Nnodes_log10.pdf"))
# plot(log10(as.numeric(as.vector(tableDegreeNnodes$N_nodes))), log10(as.numeric(as.vector(tableDegreeNnodes$degree))))
# dev.off()

##### 6.2 CALCULATE CONNECTED COMPONENTS
print("Calculating connected components")
ptm <- proc.time() # Start the clock
multTranscript.cc.transcripts <- graph::connComp(as(g, 'graphNEL')) # graphNEL representation (= nodes and an edge list) suitable when large number of nodes and relatively few edges
proc.time() - ptm # Stop the clock
#    user    system   elapsed 
#12599.468    16.691 12614.158 
is(multTranscript.cc.transcripts) # "list"   "vector"
length(multTranscript.cc.transcripts) # 11133
head(multTranscript.cc.transcripts[[1]]) # "ENST000233" "ENST256680" "ENST256682" "ENST272102" "ENST303436" "ENST398092"
saveRDS(file=paste0(path, "/multTranscript.cc.transcripts.RDS"), multTranscript.cc.transcripts)


######################################################################################################################################################################
####      7. READ AND MERGE RNA AND PROT INCIDENCE MATRICES (required for step 8. where subsequences/peptides are extracted for each connected component)        ####
######################################################################################################################################################################

##### 7.1 READ BACK FULL PROT and RNA INCIDENCE MATRICES
incidenceMatrixPSM=readRDS(file=paste0(path, "/incidenceMatrix_PSM_Rownamessimplified.RDS"))
is(incidenceMatrixPSM) # "data.frame" "list"       "oldClass"   "vector"    
dim(incidenceMatrixPSM) # 97351 53031
object.size(incidenceMatrixPSM) #20663714680bytes

incidenceMatrixRNA=readRDS(file=paste0(path, "/incidenceMatrix_RNA_Rownamessimplified.RDS"))
is(incidenceMatrixRNA) # "data.frame" "list"       "oldClass"   "vector"    
dim(incidenceMatrixRNA) # 68595 57247
object.size(incidenceMatrixRNA) # 15719379144 bytes

##### 7.2 MERGE PROT and RNA INCIDENCE MATRICES
incidenceMatrixPROTandRNA <- dplyr::bind_rows(incidenceMatrixPSM, incidenceMatrixRNA)
rownames(incidenceMatrixPROTandRNA) <- c(as.character(as.vector(rownames(incidenceMatrixPSM))),as.character(as.vector(rownames(incidenceMatrixRNA))))
saveRDS(incidenceMatrixPROTandRNA, file=paste0(path, "/incidenceMatrixPROTandRNA.RDS" ), compress="xz")
is(incidenceMatrixPROTandRNA) # "data.frame" "list"       "oldClass"   "vector"    
dim(incidenceMatrixPROTandRNA) # 165946 88073
object.size(incidenceMatrixPROTandRNA) # 58483342720 bytes

##### 7.3 REPLACE NAs BY 0s (NAs WERE INTRODUCES BY dplyr::bind_rows FUNCTION FOR MERGE) AND CONVERT TO LOGICLA (SMALLER OBJECT)
### 7.3.1 Split to replace NAs by 0s and convert to logical
SplitSizeCols <- round(dim(incidenceMatrixPROTandRNA)[2]/5)
incidenceMatrix1 <- as.matrix(incidenceMatrixPROTandRNA[,1:SplitSizeCols])
incidenceMatrix1[is.na(incidenceMatrix1)] <-0
incidenceMatrix1_logical <- apply(incidenceMatrix1,2,as.logical)
rm(incidenceMatrix1)
gc()
incidenceMatrix2 <- as.matrix(incidenceMatrixPROTandRNA[,(SplitSizeCols+1):(SplitSizeCols*2)])
incidenceMatrix2[is.na(incidenceMatrix2)] <-0
incidenceMatrix2_logical <- apply(incidenceMatrix2,2,as.logical)
rm(incidenceMatrix2)
gc()
incidenceMatrix3 <- as.matrix(incidenceMatrixPROTandRNA[,((SplitSizeCols*2)+1):(SplitSizeCols*3)])
incidenceMatrix3[is.na(incidenceMatrix3)] <-0
incidenceMatrix3_logical <- apply(incidenceMatrix3,2,as.logical)
rm(incidenceMatrix3)
gc()
incidenceMatrix4 <- as.matrix(incidenceMatrixPROTandRNA[,((SplitSizeCols*3)+1):(SplitSizeCols*4)])
incidenceMatrix4[is.na(incidenceMatrix4)] <-0
incidenceMatrix4_logical <- apply(incidenceMatrix4,2,as.logical)
rm(incidenceMatrix4)
gc()
incidenceMatrix5 <- as.matrix(incidenceMatrixPROTandRNA[,((SplitSizeCols*4)+1):dim(incidenceMatrixPROTandRNA)[2]])
incidenceMatrix5[is.na(incidenceMatrix5)] <-0
incidenceMatrix5_logical <- apply(incidenceMatrix5,2,as.logical)
rm(incidenceMatrix5,incidenceMatrixPROTandRNA, SplitSizeCols)
gc()
### 7.3.2 Paste back chunks in unique object
incidenceMatrixPROTandRNA_NoNA_logical <-cbind(incidenceMatrix1_logical, incidenceMatrix2_logical, incidenceMatrix3_logical, incidenceMatrix4_logical, incidenceMatrix5_logical)
dim(incidenceMatrixPROTandRNA_NoNA_logical) # 165946  88073
is(incidenceMatrixPROTandRNA_NoNA_logical) # "matrix"    "array"     "mMatrix"   "structure" "vector"
is(incidenceMatrixPROTandRNA_NoNA_logical[1,]) #"logical"       "vector"        "index"         "replValue"     "numLike"       "atomicVector"  "numericVector" "replValueSp"   "Mnumeric"
object.size(incidenceMatrixPROTandRNA_NoNA_logical) # 58467789928 bytes
rm(incidenceMatrix1_logical, incidenceMatrix2_logical, incidenceMatrix3_logical, incidenceMatrix4_logical, incidenceMatrix5_logical)
gc()
### 7.3.3 Put back rownames
incidenceMatrixPROTandRNA <- readRDS(paste0(path, "/incidenceMatrixPROTandRNA.RDS"))
rownames(incidenceMatrixPROTandRNA_NoNA_logical) <- rownames(incidenceMatrixPROTandRNA)
saveRDS(file=paste0(path, "/incidenceMatrixPROTandRNA_NoNA_logical.RDS"), incidenceMatrixPROTandRNA_NoNA_logical)


################################################################################################################
####      8. GENERATE LIST OF RNA_SUBSEQUENCES/PEPTIDES AND TRANSCRIPTS FOR EACH CONNECTED COMPONENT        ####
################################################################################################################

##### 8.1 FOR SINGLE TRANSCRIPT CONNECTED COMPONENTS
# Find number of transcript NOT in multiple-transcript connected components
allTranscripts <- colnames(incidenceMatrixPROTandRNA_NoNA_logical)
length(allTranscripts) # 88073
singleTranscript.cc.transcripts=setdiff(allTranscripts,unlist(multTranscript.cc.transcripts))
length(singleTranscript.cc.transcripts) # 10,238
saveRDS(file=paste0(path, "/singleTranscript.cc.transcripts.RDS"), singleTranscript.cc.transcripts)

print("Generating list of subsequences/peptides for each single transcript CC")
ptm <- proc.time() # Start the clock
singleTranscript.cc.subseq <- list()
for(i in 1:length(singleTranscript.cc.transcripts)){
  transcriptlist <- singleTranscript.cc.transcripts[[i]]
  subX <- as.matrix(incidenceMatrixPROTandRNA_NoNA_logical[,transcriptlist])
  subseqlist <- which(rowSums(subX)!=0)
  singleTranscript.cc.subseq[[i]] <- names(subseqlist)
}
proc.time() - ptm # Stop the clock
# user  system elapsed 
#52.415  27.758  80.152 
saveRDS(file=paste0(path, "/singleTranscript.cc.subseq.RDS"), singleTranscript.cc.subseq) 

print("Generating list of sub-incidence-matrices for each single transcript CC")
ptm <- proc.time() # Start the clock
singleTranscript.cc.incMatr <- list()
for(i in 1:length(singleTranscript.cc.transcripts)){
  transcriptlist <- singleTranscript.cc.transcripts[[i]]
  subX <- as.matrix(incidenceMatrixPROTandRNA_NoNA_logical[,transcriptlist])
  singleTranscript.cc.incMatr[[i]] <- subX[rowSums(subX)!=0,]
}
proc.time() - ptm # Stop the clock
#user  system elapsed 
#49.572   0.021  49.591 
saveRDS(file=paste0(path, "/singleTranscript.cc.incMatr.RDS"), singleTranscript.cc.incMatr) 

##### 8.2 FOR MULTIPLE TRANSCRIPT CONNECTED COMPONENTS
print("Generating list of subsequences/peptides for each multiple transcript CC")
ptm <- proc.time() # Start the clock
multTranscript.cc.subseq <- list()
for(i in 1:length(multTranscript.cc.transcripts)){
  transcriptlist <- multTranscript.cc.transcripts[[i]]
  subX <- incidenceMatrixPROTandRNA_NoNA_logical[,transcriptlist]
  subseqlist <- which(rowSums(subX)!=0)
  multTranscript.cc.subseq[[i]] <- names(subseqlist)
}
proc.time() - ptm # Stop the clock
#user  system elapsed 
#160.109   4.672 164.764 
saveRDS(file=paste0(path, "/multTranscript.cc.subseq.RDS"), multTranscript.cc.subseq) 

print("Generating list of sub-incidence-matrices for each multiple transcript CC")
ptm <- proc.time() # Start the clock
multTranscript.cc.incMatr <- list()
for(i in 1:length(multTranscript.cc.transcripts)){
  transcriptlist <- multTranscript.cc.transcripts[[i]]
  subX <- incidenceMatrixPROTandRNA_NoNA_logical[,transcriptlist]
  multTranscript.cc.incMatr[[i]] <- subX[rowSums(subX)!=0,]
}
proc.time() - ptm # Stop the clock
#   user  system elapsed 
#155.576   0.031 155.623 
saveRDS(file=paste0(path, "/multTranscript.cc.incMatr.RDS"), multTranscript.cc.incMatr) 
rm(i,g,subX,subseqlist,transcriptlist) # Clean memory
gc()

##### 8.3 PASTE TOGETHER RNA_SUBSEQUENCES/PEPTIDES AND TRANSCRIPTS EXTRACTED FROM SINGLE AND MULTI CCs
all.cc.subseq <- c(multTranscript.cc.subseq, singleTranscript.cc.subseq)
all.cc.transcripts <- c(multTranscript.cc.transcripts, singleTranscript.cc.transcripts)
saveRDS(file=paste0(path, "/all.cc.subseq.RDS"), all.cc.subseq)
saveRDS(file=paste0(path, "/all.cc.transcripts.RDS"), all.cc.transcripts)


####################################################################################################################
####                                  9. STATISTICS ON ALL CONNECTED COMPONENTS                                    ####
####################################################################################################################

##### 9.1 CCs SIZE (=NB TRANSCRIPTS) DISTRIBUTION: WRITE AND PLOT NB OF CCs WITH 1, 2, ..., 10, >10 TRANSCRIPTS
N_transcripts <- unlist(lapply(all.cc.transcripts, length))
N_subseqpeptides <- unlist(lapply(all.cc.subseq, length))
all.cc.summary <- as.data.frame(cbind(N_transcripts, N_subseqpeptides))
colnames(all.cc.summary) <- c("Nb_transcripts","Nb_subseqpeptides")
write.csv2(file=paste0(path, "/AllConnectedComponents_NTranscripts_Nsubseqpeptides.csv"), all.cc.summary, row.names = F)
a1 <- dim(all.cc.summary[all.cc.summary$Nb_transcripts==1,])[1]
a2 <- dim(all.cc.summary[all.cc.summary$Nb_transcripts==2,])[1]
a3 <- dim(all.cc.summary[all.cc.summary$Nb_transcripts==3,])[1]
a4 <- dim(all.cc.summary[all.cc.summary$Nb_transcripts==4,])[1]
a5 <- dim(all.cc.summary[all.cc.summary$Nb_transcripts==5,])[1]
a6 <- dim(all.cc.summary[all.cc.summary$Nb_transcripts==6,])[1]
a7 <- dim(all.cc.summary[all.cc.summary$Nb_transcripts==7,])[1]
a8 <- dim(all.cc.summary[all.cc.summary$Nb_transcripts==8,])[1]
a9 <- dim(all.cc.summary[all.cc.summary$Nb_transcripts==9,])[1]
a10 <- dim(all.cc.summary[all.cc.summary$Nb_transcripts==10,])[1]
a10plus <- dim(all.cc.summary[all.cc.summary$Nb_transcripts>10,])[1]
NtranscriptsDistribution <- as.data.frame(cbind(c("1","2","3","4","5","6","7","8","9","10",">10"),as.numeric(as.vector(c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a10plus)))))
colnames(NtranscriptsDistribution)=c("Nb_transcripts", "Nb_CC")
write.csv2(file=paste0(path, "/AllConnectedComponents_NTranscriptsDistribution.csv"), NtranscriptsDistribution, row.names = F)
pdf(paste0(path, "/AllConnectedComponents_NTranscriptsDistribution.pdf"))
plot(factor(NtranscriptsDistribution$Nb_transcripts, levels=c("1","2", "3", "4", "5", "6", "7", "8", "9", "10", ">10")), as.numeric(as.vector(NtranscriptsDistribution$Nb_CC)), type="s", xlab="Nb_transcripts", ylab="N_CCs")
dev.off()

##### DEPRECATED (NB RNA_SUBSEQUENCES NOT COMPARABLE TO NB_PEPTIDES BECAUSE MUCH LONGER AS THEY ARE THE MERGE OF ADJACENT READS)
# Plot nb transcripts vs nb RNA_subseqs+peptides
# pdf(paste0(path, "/AllConnectedComponents_NTranscripts_Nsubseqpeptides.pdf"))
# plot(all.cc.summary, pch=20)
# dev.off()
# # Zoom plot nb transcripts vs nb RNA_subseqs+peptides on CCs with 1-100 Transcripts vs 1-100 subseq/peptides 
# pdf(paste0(path, "/AllConnectedComponents_NTranscripts_Nsubseqpeptides_ZOOM_100transcripts_100subseqpeptides.pdf"))
# plot(all.cc.summary, xlim=c(1,100), ylim=c(1,100), pch=20)
# dev.off()

##### 9.2 COUNT NB OF CCs WITH ONLY RNA_SUBSEQ, ONLY PEPTIDES OR RNA_SUBSEQ+PEPTIDES 
m <- as.data.frame(matrix(nrow = length(all.cc.transcripts), ncol = 4))
colnames(m)=c("ConnCompID", "N_transcripts", "N_peptides", "N_RNAsubsequences")
for(i in 1:length(all.cc.transcripts)){
   m[i,1]=i
   m[i,2]=length(unlist(all.cc.transcripts[[i]]))
   m[i,3]=length(grep("P", all.cc.subseq[[i]]))
   m[i,4]=length(grep("R", all.cc.subseq[[i]]))
}
saveRDS(file=paste0(path, "/ConnCompID_Ntranscripts_Npeptides_Nsubseq.RDS"), m)
write.table(file=paste0(path, "/ConnCompID_Ntranscripts_Npeptides_Nsubseq.txt"), m, sep="\t", col.names = T, row.names = F, quote = F)

CC_NoRNAsubseq <- m[m$N_RNAsubsequences==0,]
dim(CC_NoRNAsubseq[CC_NoRNAsubseq$N_transcripts==1,])[1] # 770 single transcript CCs with peptides but no RNA subsequences
dim(CC_NoRNAsubseq[CC_NoRNAsubseq$N_transcripts>1,])[1] # 3224 multiple transcript CCs with peptides but no RNA subsequences
write.table(file=paste0(path, "/CC_NoRNAsubseq.txt"), CC_NoRNAsubseq, sep="\t", row.names=F, col.names=F, quote=F)

CC_RNAandPeptides <- m[((m$N_RNAsubsequences>0)&(m$N_peptides>0)),]
dim(CC_RNAandPeptides[CC_RNAandPeptides$N_transcripts==1,])[1] # 343 single transcript CCs with peptides and RNA subsequences
dim(CC_RNAandPeptides[CC_RNAandPeptides$N_transcripts>1,])[1] # 4790 multiple transcript CCs with peptides and RNA subsequences
write.table(file=paste0(path, "/CC_RNAandPeptides.txt"), CC_RNAandPeptides, sep="\t", row.names=F, col.names=F, quote=F)

CC_NoPeptides <- m[m$N_peptides==0,]
dim(CC_NoPeptides[CC_NoPeptides$N_transcripts==1,])[1] #  9125 single transcript CCs with only RNA subsequences and no peptides
dim(CC_NoPeptides[CC_NoPeptides$N_transcripts>1,])[1] #  3119 multiple transcript CCs with only RNA subsequences and no peptides
write.table(file=paste0(path, "/CC_NoPeptides.txt"), CC_NoPeptides, sep="\t", row.names=F, col.names=F, quote=F)


####################################################################################################################
####                                              10. FURTHER ANALYSES                                          ####
####################################################################################################################
PSMscore <- read.table(paste0(path, "/AllValidatedPSMsFromProteinSets_Ensembl_03122019_ContamHitReadableByR.txt"), sep="\t", header=T, quote="")
PSMscore$sequence <- as.character(as.vector(PSMscore$sequence))

RowNamesPROTSimplified <- as.data.frame(read.table(file=paste0(path, "/RowNamesPROTSimplified.txt"), sep="\t", header=T))

proteinID_geneID_transcriptID <- as.data.frame(readRDS(paste0(path, "/proteinID_geneID_transcriptID.RDS")))
proteinID_geneID_transcriptID[,1] <- stringr::str_replace(proteinID_geneID_transcriptID[,1], "ENSP00000", "ENSP")
proteinID_geneID_transcriptID[,2] <- stringr::str_replace(proteinID_geneID_transcriptID[,2], "ENST00000", "ENST")
proteinID_geneID_transcriptID[,3] <- stringr::str_replace(proteinID_geneID_transcriptID[,3], "ENSG00000", "ENSG")
proteinID_geneID_transcriptID=as.data.frame(proteinID_geneID_transcriptID)
#Error found in Proline output generating the proteinID_geneID_transcriptID object: ENSP00000375660 annotated to gene 
#ENG00000204604 instead of ENSG00000182986. Correct manually
proteinID_geneID_transcriptID[proteinID_geneID_transcriptID$ProteinID=="ENSP375660",]$GeneID="ENSG182986"

TranscriptToGeneIDsOriginalIncM_forRNA <- readRDS(paste0(path, "/TranscriptToGeneIDsOriginalIncM.RDS"))
colnames(TranscriptToGeneIDsOriginalIncM_forRNA) <- c("TranscriptID", "GeneID")

TranscriptID_GeneID_RNAandPROT <- unique(rbind(TranscriptToGeneIDsOriginalIncM_forRNA, proteinID_geneID_transcriptID[,c(2,3)]))
saveRDS(file=paste0(path, "/TranscriptID_GeneID_RNAandPROT.RDS"), TranscriptID_GeneID_RNAandPROT)

##### 10.1 DEFINE FUNCTION TO CREATE DATAFRAME WITH CONNCOMPID, N_PEPTIDES, N_RNASUBSEQS, N_TRANSCRIPTS, PSM SCORE FOR EACH CC
extract_CCs_features <- function(CC_df, all.cc.transcripts, all.cc.subseq, RowNamesPROTSimplified, PSMscore, TranscriptID_GeneID_RNAandPROT){
  CC_details <- matrix(nrow = length(CC_df$ConnCompID), ncol = 7)
  colnames(CC_details) <- c("ConnCompID", "N_peptides", "N_transcripts", "avg_PSMscore", "avg_maxPSMscorePerPeptide", "N_PSMs", "N_genes")
  AllMaxPSMscores <- list()
  for(c in 1:length(CC_df$ConnCompID)){
    ConnCompID <- CC_df[c,]$ConnCompID
    transcript_ids <- all.cc.transcripts[[ConnCompID]]
    N_transcripts <- length(transcript_ids)
    peptide_ids <- all.cc.subseq[[ConnCompID]][grep("P_", all.cc.subseq[[ConnCompID]])]
    N_peptides <- length(peptide_ids)
    peptides <- vector()
    for(p_index in 1:length(peptide_ids)){
      p <- as.character(as.vector(RowNamesPROTSimplified[RowNamesPROTSimplified$peptide_id==peptide_ids[p_index],]$peptide_seq))
      peptides <- c(peptides, p)
    }
    psm_scores <- vector()
    max_psm_scores <- vector()
    for(p_index in 1:length(peptides)){
      p_score <- as.numeric(as.vector(PSMscore[PSMscore$sequence==peptides[p_index],]$psm_score))
      psm_scores <- c(psm_scores, p_score)
      max_psm_scores <- c(max_psm_scores, max(p_score))
    }
    AllMaxPSMscores[[c]] <- list(max_psm_scores=max_psm_scores, N_transcripts=N_transcripts)
    names(AllMaxPSMscores)[c] <- ConnCompID
    N_PSMs <- length(psm_scores)
    avg_PSMscore <- mean(psm_scores)
    max_avg_PSMscore <- mean(max_psm_scores)
    genes <- vector()
    for(t_index in 1:length(transcript_ids)){
    t <- stringr::str_replace(transcript_ids[t_index], "ENST", "ENST00000")
    gene <- as.character(as.vector(TranscriptID_GeneID_RNAandPROT[TranscriptID_GeneID_RNAandPROT$TranscriptID==t,]$GeneID))
    genes <- c(genes, gene)
    }
    N_genes <- length(unique(genes))
    CC_details[c,] <- c(ConnCompID, N_peptides, N_transcripts, avg_PSMscore, max_avg_PSMscore, N_PSMs, N_genes)
  }  
  CC_details <- as.data.frame(CC_details)
  return(list(CC_details=CC_details,  AllMaxPSMscores=AllMaxPSMscores))
}

##### 10.2 APPLY FUNCTION TO CREATE DATAFRAME WITH CONNCOMPID, N_PEPTIDES, N_RNASUBSEQS, N_TRANSCRIPTS, PSM SCORE FOR:
#####     a. CCs WITH NO RNA_SUBSEQS
#####     b. CCs WITH BOTH RNA_SUBSEQs AND PEPTIDES
CC_NoRNAsubseq_res <- extract_CCs_features(CC_NoRNAsubseq, all.cc.transcripts, all.cc.subseq, RowNamesPROTSimplified, PSMscore, proteinID_geneID_transcriptID)
CC_NoRNAsubseq_details <- CC_NoRNAsubseq_res$CC_details
saveRDS(file=paste0(path, "/CC_NoRNAsubseq_details.RDS"), CC_NoRNAsubseq_details)
AllMaxPSMscores_NoRNAsubseq <- CC_NoRNAsubseq_res$AllMaxPSMscores
saveRDS(file=paste0(path, "/AllMaxPSMscores_NoRNAsubseq.RDS"), AllMaxPSMscores_NoRNAsubseq)
rm(CC_details, CC_df, ConnCompID, c, AllMaxPSMscores, transcript_ids, N_transcripts, peptide_ids, p, N_peptides, peptides, psm_scores, max_psm_scores, p_score, N_PSMs, max_avg_PSMscore, avg_PSMscore, genes, gene, N_genes, p_index, t, t_index, CC_NoRNAsubseq_res)
gc()

CC_RNAandPeptides_res <- extract_CCs_features(CC_RNAandPeptides, all.cc.transcripts, all.cc.subseq, RowNamesPROTSimplified, PSMscore, proteinID_geneID_transcriptID)
CC_RNAandPeptides_details <- CC_RNAandPeptides_res$CC_details
saveRDS(file=paste0(path, "/CC_RNAandPeptides_details.RDS"), CC_RNAandPeptides_details)
AllMaxPSMscores_RNAandPeptides <- CC_RNAandPeptides_res$AllMaxPSMscores
saveRDS(file=paste0(path, "/AllMaxPSMscores_RNAandPeptides.RDS"), AllMaxPSMscores_RNAandPeptides)
rm(CC_details, CC_df, ConnCompID, c, AllMaxPSMscores,transcript_ids, N_transcripts, peptide_ids, p, N_peptides, peptides, psm_scores, max_psm_scores, p_score, N_PSMs, max_avg_PSMscore, avg_PSMscore, genes, gene, N_genes, p_index, t, t_index, CC_RNAandPeptides_res)
gc()


##### 10.3 COMPARE SCORES CCs WITH NO RNA_SUBSEQS VS CCs WITH BOTH RNA_SUBSEQs AND PEPTIDES
# Define function to compare average PSM score (all PSMs for all peptides) or average maximum PSM score per peptide
averagePSMscore <-function(df1, df1_name, df2, df2_name, scoreType){
  BoxplotPSMscores <- as.data.frame(rbind(cbind(rep(df1_name, length(df1[[scoreType]])), df1[[scoreType]]), cbind(rep(df2_name, length(df2[[scoreType]])), df2[[scoreType]])))
  colnames(BoxplotPSMscores) <-c("CC_type", scoreType)
  title <- paste0("Nb CCs ", df1_name, "=", length(unique(df1$ConnCompID)), "; Nb CCs ", df2_name, "=", length(unique(df2$ConnCompID)))
  subtitle <- paste0("Nb peptides ", df1_name, "=", sum(df1$N_peptides), "; Nb peptides ", df2_name, "=", sum(df2$N_peptides))
  pdf(paste0("Boxplot_", scoreType, "_", df1_name, "_", df2_name, ".pdf"))
  boxplot(as.numeric(as.vector(BoxplotPSMscores[[scoreType]]))~BoxplotPSMscores$CC_type, data=BoxplotPSMscores, xlab="CC type", ylab=scoreType, main=title, cex.main=1, sub=subtitle, cex.sub=0.7)
  dev.off()
}

# Define function to compare maximum PSM score per peptide (one PSM score per peptide, all peptides for each CC)
maxPSMscoreAllPeptides <- function(df1, df1_name, AllMaxPSMscores_df1, df2, df2_name, AllMaxPSMscores_df2){
  #AllMaxPSMscores objects are list with the best PSMs for each of their peptides, separately for each CC. Concatenate them all together in a dataframe and indicate only 
  #if they come from NoRNAsubseq type CCs or RNAandPeptides type CCs.
  df1_MaxPSMscores <- vector()
  for(i in 1:length(AllMaxPSMscores_df1)){
    df1_MaxPSMscores <- c(df1_MaxPSMscores, as.numeric(as.vector(unlist(AllMaxPSMscores_df1[[i]]["max_psm_scores"]))))
  }
  df1_MaxPSMscores <- as.data.frame(cbind(as.character(as.vector(rep(df1_name, length(df1_MaxPSMscores)))), as.numeric(as.vector(df1_MaxPSMscores))))
  colnames(df1_MaxPSMscores) <- c("CC_type", "MaxPSMscorePerPeptide")
  df2_MaxPSMscores <- vector()
  for(i in 1:length(AllMaxPSMscores_df2)){
    df2_MaxPSMscores <- c(df2_MaxPSMscores, as.numeric(as.vector(unlist(AllMaxPSMscores_df2[[i]]["max_psm_scores"]))))
  }
  df2_MaxPSMscores <- as.data.frame(cbind(as.character(as.vector(rep(df2_name, length(df2_MaxPSMscores)))), as.numeric(as.vector(unlist(df2_MaxPSMscores)))))
  colnames(df2_MaxPSMscores)  <- c("CC_type", "MaxPSMscorePerPeptide")
  #Generate boxplot
  BoxplotAllPeptidesMaxPSMscore <- rbind(df1_MaxPSMscores, df2_MaxPSMscores)
  title <- paste0("Nb CCs ", df1_name, "=", length(unique(df1$ConnCompID)), "; Nb CCs ", df2_name, "=", length(unique(df2$ConnCompID)))
  subtitle <- paste0("Nb peptides ", df1_name, "=", sum(df1$N_peptides), "; Nb peptides ", df2_name, "=", sum(df2$N_peptides))
  pdf(paste0("BoxplotAllPeptidesMaxPSMscore_", df1_name, "_", df2_name, ".pdf"))
  boxplot(as.numeric(as.vector(BoxplotAllPeptidesMaxPSMscore$MaxPSMscorePerPeptide))~BoxplotAllPeptidesMaxPSMscore$CC_type, data=BoxplotAllPeptidesMaxPSMscore, xlab="CC type", ylab="Max PSM score Per Peptide", main=title, cex.main=1, sub=subtitle, cex.sub=0.7)
  dev.off()
  pdf(paste0("BoxplotAllPeptidesMaxPSMscore_", df1_name, "_", df2_name, "_log10.pdf"))
  boxplot(log10(as.numeric(as.vector(BoxplotAllPeptidesMaxPSMscore$MaxPSMscorePerPeptide)))~BoxplotAllPeptidesMaxPSMscore$CC_type, data=BoxplotAllPeptidesMaxPSMscore, xlab="CC type", ylab="log10 max PSM score per peptide", main=title, cex.main=1, sub=subtitle, cex.sub=0.7)
  dev.off()
}

### 10.3.1 Compare CCs_NoRNAsubseqs versus CCs_RNAsubseqAndPeptides using all CCs (single+multi) 
# Average PSM score (all PSMs for all peptides)
averagePSMscore(CC_NoRNAsubseq_details, "CCNoRNAsubseq", CC_RNAandPeptides_details, "CCRNAandPeptides", "avg_PSMscore")
# Average maximum PSM score per peptide (best PSM score per peptide)
averagePSMscore(CC_NoRNAsubseq_details, "CCNoRNAsubseq", CC_RNAandPeptides_details, "CCRNAandPeptides", "avg_maxPSMscorePerPeptide")
# Maximum PSM score per peptide (one PSM score per peptide, all peptides for each CC)
maxPSMscoreAllPeptides(CC_NoRNAsubseq_details, "CCNoRNAsubseq", AllMaxPSMscores_NoRNAsubseq, CC_RNAandPeptides_details, "CCRNAandPeptides", AllMaxPSMscores_RNAandPeptides)
# Maximum PSM score per peptide (one PSM score per peptide, all peptides for each CC) in CC_NoRNAsubseq with only 1 peptide, 2, 3 or >3 peptides
CC_NoRNAsubseq_MaxPSMscores_pep1 <- vector()
CC_NoRNAsubseq_MaxPSMscores_pep2 <- vector()
CC_NoRNAsubseq_MaxPSMscores_pep3 <- vector()
CC_NoRNAsubseq_MaxPSMscores_pep3More <- vector()
CC_NoRNAsubseq_MaxPSMscores_pep4 <- vector()
CC_NoRNAsubseq_MaxPSMscores_pep5 <- vector()
CC_NoRNAsubseq_MaxPSMscores_pep6 <- vector()
CC_NoRNAsubseq_MaxPSMscores_pep7 <- vector()
CC_NoRNAsubseq_MaxPSMscores_pep8 <- vector()
CC_NoRNAsubseq_MaxPSMscores_pep9 <- vector()
CC_NoRNAsubseq_MaxPSMscores_pep10 <- vector()
CC_NoRNAsubseq_MaxPSMscores_pep10More <- vector()
for(i in 1:length(AllMaxPSMscores_NoRNAsubseq)){
  if(length(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==1){
    CC_NoRNAsubseq_MaxPSMscores_pep1 <- c(CC_NoRNAsubseq_MaxPSMscores_pep1, as.numeric(as.vector(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==2){
    CC_NoRNAsubseq_MaxPSMscores_pep2 <- c(CC_NoRNAsubseq_MaxPSMscores_pep2, as.numeric(as.vector(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==3){
    CC_NoRNAsubseq_MaxPSMscores_pep3 <- c(CC_NoRNAsubseq_MaxPSMscores_pep3, as.numeric(as.vector(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)>3){
    CC_NoRNAsubseq_MaxPSMscores_pep3More <- c(CC_NoRNAsubseq_MaxPSMscores_pep3More, as.numeric(as.vector(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==4){
    CC_NoRNAsubseq_MaxPSMscores_pep4 <- c(CC_NoRNAsubseq_MaxPSMscores_pep4, as.numeric(as.vector(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==5){
    CC_NoRNAsubseq_MaxPSMscores_pep5 <- c(CC_NoRNAsubseq_MaxPSMscores_pep5, as.numeric(as.vector(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==6){
    CC_NoRNAsubseq_MaxPSMscores_pep6 <- c(CC_NoRNAsubseq_MaxPSMscores_pep6, as.numeric(as.vector(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==7){
    CC_NoRNAsubseq_MaxPSMscores_pep7 <- c(CC_NoRNAsubseq_MaxPSMscores_pep7, as.numeric(as.vector(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==8){
    CC_NoRNAsubseq_MaxPSMscores_pep8 <- c(CC_NoRNAsubseq_MaxPSMscores_pep8, as.numeric(as.vector(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==9){
    CC_NoRNAsubseq_MaxPSMscores_pep9 <- c(CC_NoRNAsubseq_MaxPSMscores_pep9, as.numeric(as.vector(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==10){
    CC_NoRNAsubseq_MaxPSMscores_pep10 <- c(CC_NoRNAsubseq_MaxPSMscores_pep10, as.numeric(as.vector(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)>10){
    CC_NoRNAsubseq_MaxPSMscores_pep10More <- c(CC_NoRNAsubseq_MaxPSMscores_pep10More, as.numeric(as.vector(AllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
}
CC_NoRNAsubseq_MaxPSMscores_pep1 <- as.data.frame(cbind(as.character(as.vector(rep("CC_NoRNAsubseq", length(CC_NoRNAsubseq_MaxPSMscores_pep1)))), as.numeric(as.vector(CC_NoRNAsubseq_MaxPSMscores_pep1)), as.character(as.vector(rep("1", length(CC_NoRNAsubseq_MaxPSMscores_pep1))))))
CC_NoRNAsubseq_MaxPSMscores_pep2 <- as.data.frame(cbind(as.character(as.vector(rep("CC_NoRNAsubseq", length(CC_NoRNAsubseq_MaxPSMscores_pep2)))), as.numeric(as.vector(CC_NoRNAsubseq_MaxPSMscores_pep2)), as.character(as.vector(rep("2", length(CC_NoRNAsubseq_MaxPSMscores_pep2))))))
CC_NoRNAsubseq_MaxPSMscores_pep3 <- as.data.frame(cbind(as.character(as.vector(rep("CC_NoRNAsubseq", length(CC_NoRNAsubseq_MaxPSMscores_pep3)))), as.numeric(as.vector(CC_NoRNAsubseq_MaxPSMscores_pep3)), as.character(as.vector(rep("3", length(CC_NoRNAsubseq_MaxPSMscores_pep3))))))
CC_NoRNAsubseq_MaxPSMscores_pep3More <- as.data.frame(cbind(as.character(as.vector(rep("CC_NoRNAsubseq", length(CC_NoRNAsubseq_MaxPSMscores_pep3More)))), as.numeric(as.vector(CC_NoRNAsubseq_MaxPSMscores_pep3More)), as.character(as.vector(rep(">3", length(CC_NoRNAsubseq_MaxPSMscores_pep3More))))))
CC_NoRNAsubseq_MaxPSMscores1to3 <- rbind(CC_NoRNAsubseq_MaxPSMscores_pep1, CC_NoRNAsubseq_MaxPSMscores_pep2, CC_NoRNAsubseq_MaxPSMscores_pep3, CC_NoRNAsubseq_MaxPSMscores_pep3More)
colnames(CC_NoRNAsubseq_MaxPSMscores1to3) <- c("CC_type", "MaxPSMscorePerPeptide", "N_peptides")

CC_NoRNAsubseq_MaxPSMscores_pep4 <- as.data.frame(cbind(as.character(as.vector(rep("CC_NoRNAsubseq", length(CC_NoRNAsubseq_MaxPSMscores_pep4)))), as.numeric(as.vector(CC_NoRNAsubseq_MaxPSMscores_pep4)), as.character(as.vector(rep("4", length(CC_NoRNAsubseq_MaxPSMscores_pep4))))))
CC_NoRNAsubseq_MaxPSMscores_pep5 <- as.data.frame(cbind(as.character(as.vector(rep("CC_NoRNAsubseq", length(CC_NoRNAsubseq_MaxPSMscores_pep5)))), as.numeric(as.vector(CC_NoRNAsubseq_MaxPSMscores_pep5)), as.character(as.vector(rep("5", length(CC_NoRNAsubseq_MaxPSMscores_pep5))))))
CC_NoRNAsubseq_MaxPSMscores_pep6 <- as.data.frame(cbind(as.character(as.vector(rep("CC_NoRNAsubseq", length(CC_NoRNAsubseq_MaxPSMscores_pep6)))), as.numeric(as.vector(CC_NoRNAsubseq_MaxPSMscores_pep6)), as.character(as.vector(rep("6", length(CC_NoRNAsubseq_MaxPSMscores_pep6))))))
CC_NoRNAsubseq_MaxPSMscores_pep7 <- as.data.frame(cbind(as.character(as.vector(rep("CC_NoRNAsubseq", length(CC_NoRNAsubseq_MaxPSMscores_pep7)))), as.numeric(as.vector(CC_NoRNAsubseq_MaxPSMscores_pep7)), as.character(as.vector(rep("7", length(CC_NoRNAsubseq_MaxPSMscores_pep7))))))
CC_NoRNAsubseq_MaxPSMscores_pep8 <- as.data.frame(cbind(as.character(as.vector(rep("CC_NoRNAsubseq", length(CC_NoRNAsubseq_MaxPSMscores_pep8)))), as.numeric(as.vector(CC_NoRNAsubseq_MaxPSMscores_pep8)), as.character(as.vector(rep("8", length(CC_NoRNAsubseq_MaxPSMscores_pep8))))))
CC_NoRNAsubseq_MaxPSMscores_pep9 <- as.data.frame(cbind(as.character(as.vector(rep("CC_NoRNAsubseq", length(CC_NoRNAsubseq_MaxPSMscores_pep9)))), as.numeric(as.vector(CC_NoRNAsubseq_MaxPSMscores_pep9)), as.character(as.vector(rep("9", length(CC_NoRNAsubseq_MaxPSMscores_pep9))))))
CC_NoRNAsubseq_MaxPSMscores_pep10 <- as.data.frame(cbind(as.character(as.vector(rep("CC_NoRNAsubseq", length(CC_NoRNAsubseq_MaxPSMscores_pep10)))), as.numeric(as.vector(CC_NoRNAsubseq_MaxPSMscores_pep10)), as.character(as.vector(rep("10", length(CC_NoRNAsubseq_MaxPSMscores_pep10))))))
CC_NoRNAsubseq_MaxPSMscores_pep10More <- as.data.frame(cbind(as.character(as.vector(rep("CC_NoRNAsubseq", length(CC_NoRNAsubseq_MaxPSMscores_pep10More)))), as.numeric(as.vector(CC_NoRNAsubseq_MaxPSMscores_pep10More)), as.character(as.vector(rep(">10", length(CC_NoRNAsubseq_MaxPSMscores_pep10More))))))
CC_NoRNAsubseq_MaxPSMscores1to10 <- rbind(CC_NoRNAsubseq_MaxPSMscores_pep1, CC_NoRNAsubseq_MaxPSMscores_pep2, CC_NoRNAsubseq_MaxPSMscores_pep3, CC_NoRNAsubseq_MaxPSMscores_pep4, CC_NoRNAsubseq_MaxPSMscores_pep5, CC_NoRNAsubseq_MaxPSMscores_pep6, CC_NoRNAsubseq_MaxPSMscores_pep7, CC_NoRNAsubseq_MaxPSMscores_pep8, CC_NoRNAsubseq_MaxPSMscores_pep9, CC_NoRNAsubseq_MaxPSMscores_pep10, CC_NoRNAsubseq_MaxPSMscores_pep10More)
colnames(CC_NoRNAsubseq_MaxPSMscores1to10) <- c("CC_type", "MaxPSMscorePerPeptide", "N_peptides")

CC_RNAandPeptides_MaxPSMscores_pep1 <- vector()
CC_RNAandPeptides_MaxPSMscores_pep2 <- vector()
CC_RNAandPeptides_MaxPSMscores_pep3 <- vector()
CC_RNAandPeptides_MaxPSMscores_pep3More <- vector()
CC_RNAandPeptides_MaxPSMscores_pep4 <- vector()
CC_RNAandPeptides_MaxPSMscores_pep5 <- vector()
CC_RNAandPeptides_MaxPSMscores_pep6 <- vector()
CC_RNAandPeptides_MaxPSMscores_pep7 <- vector()
CC_RNAandPeptides_MaxPSMscores_pep8 <- vector()
CC_RNAandPeptides_MaxPSMscores_pep9 <- vector()
CC_RNAandPeptides_MaxPSMscores_pep10 <- vector()
CC_RNAandPeptides_MaxPSMscores_pep10More <- vector()
for(i in 1:length(AllMaxPSMscores_RNAandPeptides)){
  if(length(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==1){
    CC_RNAandPeptides_MaxPSMscores_pep1 <- c(CC_RNAandPeptides_MaxPSMscores_pep1, as.numeric(as.vector(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==2){
    CC_RNAandPeptides_MaxPSMscores_pep2 <- c(CC_RNAandPeptides_MaxPSMscores_pep2, as.numeric(as.vector(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==3){
    CC_RNAandPeptides_MaxPSMscores_pep3 <- c(CC_RNAandPeptides_MaxPSMscores_pep3, as.numeric(as.vector(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)>3){
    CC_RNAandPeptides_MaxPSMscores_pep3More <- c(CC_RNAandPeptides_MaxPSMscores_pep3More, as.numeric(as.vector(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==4){
    CC_RNAandPeptides_MaxPSMscores_pep4 <- c(CC_RNAandPeptides_MaxPSMscores_pep4, as.numeric(as.vector(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==5){
    CC_RNAandPeptides_MaxPSMscores_pep5 <- c(CC_RNAandPeptides_MaxPSMscores_pep5, as.numeric(as.vector(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==6){
    CC_RNAandPeptides_MaxPSMscores_pep6 <- c(CC_RNAandPeptides_MaxPSMscores_pep6, as.numeric(as.vector(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==7){
    CC_RNAandPeptides_MaxPSMscores_pep7 <- c(CC_RNAandPeptides_MaxPSMscores_pep7, as.numeric(as.vector(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==8){
    CC_RNAandPeptides_MaxPSMscores_pep8 <- c(CC_RNAandPeptides_MaxPSMscores_pep8, as.numeric(as.vector(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==9){
    CC_RNAandPeptides_MaxPSMscores_pep9 <- c(CC_RNAandPeptides_MaxPSMscores_pep9, as.numeric(as.vector(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==10){
    CC_RNAandPeptides_MaxPSMscores_pep10 <- c(CC_RNAandPeptides_MaxPSMscores_pep10, as.numeric(as.vector(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)>10){
    CC_RNAandPeptides_MaxPSMscores_pep10More <- c(CC_RNAandPeptides_MaxPSMscores_pep10More, as.numeric(as.vector(AllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
}
CC_RNAandPeptides_MaxPSMscores_pep1 <- as.data.frame(cbind(as.character(as.vector(rep("CC_RNAandPeptides", length(CC_RNAandPeptides_MaxPSMscores_pep1)))), as.numeric(as.vector(CC_RNAandPeptides_MaxPSMscores_pep1)), as.character(as.vector(rep("1", length(CC_RNAandPeptides_MaxPSMscores_pep1))))))
CC_RNAandPeptides_MaxPSMscores_pep2 <- as.data.frame(cbind(as.character(as.vector(rep("CC_RNAandPeptides", length(CC_RNAandPeptides_MaxPSMscores_pep2)))), as.numeric(as.vector(CC_RNAandPeptides_MaxPSMscores_pep2)), as.character(as.vector(rep("2", length(CC_RNAandPeptides_MaxPSMscores_pep2))))))
CC_RNAandPeptides_MaxPSMscores_pep3 <- as.data.frame(cbind(as.character(as.vector(rep("CC_RNAandPeptides", length(CC_RNAandPeptides_MaxPSMscores_pep3)))), as.numeric(as.vector(CC_RNAandPeptides_MaxPSMscores_pep3)), as.character(as.vector(rep("3", length(CC_RNAandPeptides_MaxPSMscores_pep3))))))
CC_RNAandPeptides_MaxPSMscores_pep3More <- as.data.frame(cbind(as.character(as.vector(rep("CC_RNAandPeptides", length(CC_RNAandPeptides_MaxPSMscores_pep3More)))), as.numeric(as.vector(CC_RNAandPeptides_MaxPSMscores_pep3More)), as.character(as.vector(rep(">3", length(CC_RNAandPeptides_MaxPSMscores_pep3More))))))
CC_RNAandPeptides_MaxPSMscores1to3 <- rbind(CC_RNAandPeptides_MaxPSMscores_pep1, CC_RNAandPeptides_MaxPSMscores_pep2, CC_RNAandPeptides_MaxPSMscores_pep3, CC_RNAandPeptides_MaxPSMscores_pep3More)
colnames(CC_RNAandPeptides_MaxPSMscores1to3) <- c("CC_type", "MaxPSMscorePerPeptide", "N_peptides")

CC_RNAandPeptides_MaxPSMscores_pep4 <- as.data.frame(cbind(as.character(as.vector(rep("CC_RNAandPeptides", length(CC_RNAandPeptides_MaxPSMscores_pep4)))), as.numeric(as.vector(CC_RNAandPeptides_MaxPSMscores_pep4)), as.character(as.vector(rep("4", length(CC_RNAandPeptides_MaxPSMscores_pep4))))))
CC_RNAandPeptides_MaxPSMscores_pep5 <- as.data.frame(cbind(as.character(as.vector(rep("CC_RNAandPeptides", length(CC_RNAandPeptides_MaxPSMscores_pep5)))), as.numeric(as.vector(CC_RNAandPeptides_MaxPSMscores_pep5)), as.character(as.vector(rep("5", length(CC_RNAandPeptides_MaxPSMscores_pep5))))))
CC_RNAandPeptides_MaxPSMscores_pep6 <- as.data.frame(cbind(as.character(as.vector(rep("CC_RNAandPeptides", length(CC_RNAandPeptides_MaxPSMscores_pep6)))), as.numeric(as.vector(CC_RNAandPeptides_MaxPSMscores_pep6)), as.character(as.vector(rep("6", length(CC_RNAandPeptides_MaxPSMscores_pep6))))))
CC_RNAandPeptides_MaxPSMscores_pep7 <- as.data.frame(cbind(as.character(as.vector(rep("CC_RNAandPeptides", length(CC_RNAandPeptides_MaxPSMscores_pep7)))), as.numeric(as.vector(CC_RNAandPeptides_MaxPSMscores_pep7)), as.character(as.vector(rep("7", length(CC_RNAandPeptides_MaxPSMscores_pep7))))))
CC_RNAandPeptides_MaxPSMscores_pep8 <- as.data.frame(cbind(as.character(as.vector(rep("CC_RNAandPeptides", length(CC_RNAandPeptides_MaxPSMscores_pep8)))), as.numeric(as.vector(CC_RNAandPeptides_MaxPSMscores_pep8)), as.character(as.vector(rep("8", length(CC_RNAandPeptides_MaxPSMscores_pep8))))))
CC_RNAandPeptides_MaxPSMscores_pep9 <- as.data.frame(cbind(as.character(as.vector(rep("CC_RNAandPeptides", length(CC_RNAandPeptides_MaxPSMscores_pep9)))), as.numeric(as.vector(CC_RNAandPeptides_MaxPSMscores_pep9)), as.character(as.vector(rep("9", length(CC_RNAandPeptides_MaxPSMscores_pep9))))))
CC_RNAandPeptides_MaxPSMscores_pep10 <- as.data.frame(cbind(as.character(as.vector(rep("CC_RNAandPeptides", length(CC_RNAandPeptides_MaxPSMscores_pep10)))), as.numeric(as.vector(CC_RNAandPeptides_MaxPSMscores_pep10)), as.character(as.vector(rep("10", length(CC_RNAandPeptides_MaxPSMscores_pep10))))))
CC_RNAandPeptides_MaxPSMscores_pep10More <- as.data.frame(cbind(as.character(as.vector(rep("CC_RNAandPeptides", length(CC_RNAandPeptides_MaxPSMscores_pep10More)))), as.numeric(as.vector(CC_RNAandPeptides_MaxPSMscores_pep10More)), as.character(as.vector(rep(">10", length(CC_RNAandPeptides_MaxPSMscores_pep10More))))))
CC_RNAandPeptides_MaxPSMscores1to10 <- rbind(CC_RNAandPeptides_MaxPSMscores_pep1, CC_RNAandPeptides_MaxPSMscores_pep2, CC_RNAandPeptides_MaxPSMscores_pep3, CC_RNAandPeptides_MaxPSMscores_pep4, CC_RNAandPeptides_MaxPSMscores_pep5, CC_RNAandPeptides_MaxPSMscores_pep6, CC_RNAandPeptides_MaxPSMscores_pep7, CC_RNAandPeptides_MaxPSMscores_pep8, CC_RNAandPeptides_MaxPSMscores_pep9, CC_RNAandPeptides_MaxPSMscores_pep10, CC_RNAandPeptides_MaxPSMscores_pep10More)
colnames(CC_RNAandPeptides_MaxPSMscores1to10) <- c("CC_type", "MaxPSMscorePerPeptide", "N_peptides")

BoxplotAllPeptidesMaxPSMscore1to3 <- rbind(CC_NoRNAsubseq_MaxPSMscores1to3, CC_RNAandPeptides_MaxPSMscores1to3)
BoxplotAllPeptidesMaxPSMscore1to3$N_peptides <-factor(BoxplotAllPeptidesMaxPSMscore1to3$N_peptides, levels=c("1","2","3",">3"))
pdf("BoxplotAllPeptidesMaxPSMscore_CCNoRNAsubseq_CCRNAandPeptides_ByNb1to3Peptides.pdf")
ggplot(BoxplotAllPeptidesMaxPSMscore1to3, aes(x=N_peptides, y=as.numeric(as.vector(MaxPSMscorePerPeptide)), fill=CC_type)) + 
  geom_boxplot() +
  labs(x="CC_type", y="max PSM score per peptide")
dev.off()

BoxplotAllPeptidesMaxPSMscore1to10 <- rbind(CC_NoRNAsubseq_MaxPSMscores1to10, CC_RNAandPeptides_MaxPSMscores1to10)
BoxplotAllPeptidesMaxPSMscore1to10$N_peptides <-factor(BoxplotAllPeptidesMaxPSMscore1to10$N_peptides, levels=c("1","2","3","4","5","6","7","8","9","10",">10"))
pdf("BoxplotAllPeptidesMaxPSMscore_CCNoRNAsubseq_CCRNAandPeptides_ByNb1to10Peptides.pdf")
ggplot(BoxplotAllPeptidesMaxPSMscore1to10, aes(x=N_peptides, y=as.numeric(as.vector(MaxPSMscorePerPeptide)), fill=CC_type)) + 
  geom_boxplot() +
  labs(x="CC_type", y="max PSM score per peptide")
dev.off()


### 10.3.2 Compare CCs_NoRNAsubseqs versus CCs_RNAsubseqAndPeptides using only single transcript CCs:
SingleCC_NoRNAsubseq_details <- CC_NoRNAsubseq_details[as.numeric(as.vector(CC_NoRNAsubseq_details$N_transcripts))==1,] # 770 6
SingleCC_RNAandPeptides_details <- CC_RNAandPeptides_details[as.numeric(as.vector(CC_RNAandPeptides_details$N_transcripts))==1,] # 343 6
# Average PSM score (all PSMs for all peptides)
averagePSMscore(SingleCC_NoRNAsubseq_details, "SingleCCNoRNAsubseq", SingleCC_RNAandPeptides_details, "SingleCCRNAandPeptides", "avg_PSMscore")
# Maximum PSM score on average (best PSM score per peptide)
averagePSMscore(SingleCC_NoRNAsubseq_details, "SingleCCNoRNAsubseq", SingleCC_RNAandPeptides_details, "SingleCCRNAandPeptides", "avg_maxPSMscorePerPeptide")
# Maximum PSM score per peptide (one PSM score per peptide, all peptides for each CC)
SingleAllMaxPSMscores_NoRNAsubseq <- rlist::list.filter(AllMaxPSMscores_NoRNAsubseq, N_transcripts==1)
SingleAllMaxPSMscores_RNAandPeptides <- rlist::list.filter(AllMaxPSMscores_RNAandPeptides, N_transcripts==1)
maxPSMscoreAllPeptides(SingleCC_NoRNAsubseq_details, "SingleCCNoRNAsubseq", SingleAllMaxPSMscores_NoRNAsubseq, SingleCC_RNAandPeptides_details, "SingleCCRNAandPeptides", SingleAllMaxPSMscores_RNAandPeptides)
# Maximum PSM score per peptide (one PSM score per peptide, all peptides for each CC) in singleCC_NoRNAsubseq with only 1 peptide, 2, 3 or >3 peptides
SingleCC_NoRNAsubseq_MaxPSMscores_pep1 <- vector()
SingleCC_NoRNAsubseq_MaxPSMscores_pep1More <- vector()
SingleCC_NoRNAsubseq_MaxPSMscores_pep2 <- vector()
SingleCC_NoRNAsubseq_MaxPSMscores_pep3 <- vector()
SingleCC_NoRNAsubseq_MaxPSMscores_pep3More <- vector()
SingleCC_NoRNAsubseq_MaxPSMscores_pep4 <- vector()
SingleCC_NoRNAsubseq_MaxPSMscores_pep5 <- vector()
SingleCC_NoRNAsubseq_MaxPSMscores_pep6 <- vector()
SingleCC_NoRNAsubseq_MaxPSMscores_pep7 <- vector()
SingleCC_NoRNAsubseq_MaxPSMscores_pep8 <- vector()
SingleCC_NoRNAsubseq_MaxPSMscores_pep9 <- vector()
SingleCC_NoRNAsubseq_MaxPSMscores_pep10 <- vector()
SingleCC_NoRNAsubseq_MaxPSMscores_pep10More <- vector()
for(i in 1:length(SingleAllMaxPSMscores_NoRNAsubseq)){
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==1){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep1 <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep1, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)>1){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep1More <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep1More, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==2){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep2 <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep2, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==3){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep3 <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep3, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)>3){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep3More <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep3More, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==4){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep4 <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep4, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==5){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep5 <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep5, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==6){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep6 <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep6, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==7){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep7 <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep7, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==8){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep8 <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep8, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==9){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep9 <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep9, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)==10){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep10 <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep10, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)>10){
    SingleCC_NoRNAsubseq_MaxPSMscores_pep10More <- c(SingleCC_NoRNAsubseq_MaxPSMscores_pep10More, as.numeric(as.vector(SingleAllMaxPSMscores_NoRNAsubseq[[i]]$max_psm_scores)))
  }
}
SingleCC_NoRNAsubseq_MaxPSMscores_pep1 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep1)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep1)), as.character(as.vector(rep("1", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep1))))))
SingleCC_NoRNAsubseq_MaxPSMscores_pep1More <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep1More)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep1More)), as.character(as.vector(rep(">1", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep1More))))))
SingleCC_NoRNAsubseq_MaxPSMscores1HitWonderOrNot <- rbind(SingleCC_NoRNAsubseq_MaxPSMscores_pep1, SingleCC_NoRNAsubseq_MaxPSMscores_pep1More)
colnames(SingleCC_NoRNAsubseq_MaxPSMscores1HitWonderOrNot) <- c("CC_type", "MaxPSMscorePerPeptide", "N_peptides")

SingleCC_NoRNAsubseq_MaxPSMscores_pep2 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep2)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep2)), as.character(as.vector(rep("2", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep2))))))
SingleCC_NoRNAsubseq_MaxPSMscores_pep3 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep3)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep3)), as.character(as.vector(rep("3", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep3))))))
SingleCC_NoRNAsubseq_MaxPSMscores_pep3More <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep3More)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep3More)), as.character(as.vector(rep(">3", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep3More))))))
SingleCC_NoRNAsubseq_MaxPSMscores1to3 <- rbind(SingleCC_NoRNAsubseq_MaxPSMscores_pep1, SingleCC_NoRNAsubseq_MaxPSMscores_pep2, SingleCC_NoRNAsubseq_MaxPSMscores_pep3, SingleCC_NoRNAsubseq_MaxPSMscores_pep3More)
colnames(SingleCC_NoRNAsubseq_MaxPSMscores1to3) <- c("CC_type", "MaxPSMscorePerPeptide", "N_peptides")

SingleCC_NoRNAsubseq_MaxPSMscores_pep4 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep4)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep4)), as.character(as.vector(rep("4", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep4))))))
SingleCC_NoRNAsubseq_MaxPSMscores_pep5 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep5)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep5)), as.character(as.vector(rep("5", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep5))))))
SingleCC_NoRNAsubseq_MaxPSMscores_pep6 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep6)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep6)), as.character(as.vector(rep("6", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep6))))))
SingleCC_NoRNAsubseq_MaxPSMscores_pep7 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep7)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep7)), as.character(as.vector(rep("7", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep7))))))
SingleCC_NoRNAsubseq_MaxPSMscores_pep8 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep8)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep8)), as.character(as.vector(rep("8", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep8))))))
SingleCC_NoRNAsubseq_MaxPSMscores_pep9 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep9)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep9)), as.character(as.vector(rep("9", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep9))))))
SingleCC_NoRNAsubseq_MaxPSMscores_pep10 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep10)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep10)), as.character(as.vector(rep("10", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep10))))))
SingleCC_NoRNAsubseq_MaxPSMscores_pep10More <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep10More)))), as.numeric(as.vector(SingleCC_NoRNAsubseq_MaxPSMscores_pep10More)), as.character(as.vector(rep(">10", length(SingleCC_NoRNAsubseq_MaxPSMscores_pep10More))))))
SingleCC_NoRNAsubseq_MaxPSMscores1to10 <- rbind(SingleCC_NoRNAsubseq_MaxPSMscores_pep1, SingleCC_NoRNAsubseq_MaxPSMscores_pep2, SingleCC_NoRNAsubseq_MaxPSMscores_pep3, SingleCC_NoRNAsubseq_MaxPSMscores_pep4, SingleCC_NoRNAsubseq_MaxPSMscores_pep5, SingleCC_NoRNAsubseq_MaxPSMscores_pep6, SingleCC_NoRNAsubseq_MaxPSMscores_pep7, SingleCC_NoRNAsubseq_MaxPSMscores_pep8, SingleCC_NoRNAsubseq_MaxPSMscores_pep9, SingleCC_NoRNAsubseq_MaxPSMscores_pep10, SingleCC_NoRNAsubseq_MaxPSMscores_pep10More)
colnames(SingleCC_NoRNAsubseq_MaxPSMscores1to10) <- c("CC_type", "MaxPSMscorePerPeptide", "N_peptides")

SingleCC_RNAandPeptides_MaxPSMscores_pep1 <- vector()
SingleCC_RNAandPeptides_MaxPSMscores_pep1More <- vector()
SingleCC_RNAandPeptides_MaxPSMscores_pep2 <- vector()
SingleCC_RNAandPeptides_MaxPSMscores_pep3 <- vector()
SingleCC_RNAandPeptides_MaxPSMscores_pep3More <- vector()
SingleCC_RNAandPeptides_MaxPSMscores_pep4 <- vector()
SingleCC_RNAandPeptides_MaxPSMscores_pep5 <- vector()
SingleCC_RNAandPeptides_MaxPSMscores_pep6 <- vector()
SingleCC_RNAandPeptides_MaxPSMscores_pep7 <- vector()
SingleCC_RNAandPeptides_MaxPSMscores_pep8 <- vector()
SingleCC_RNAandPeptides_MaxPSMscores_pep9 <- vector()
SingleCC_RNAandPeptides_MaxPSMscores_pep10 <- vector()
SingleCC_RNAandPeptides_MaxPSMscores_pep10More <- vector()
for(i in 1:length(SingleAllMaxPSMscores_RNAandPeptides)){
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==1){
    SingleCC_RNAandPeptides_MaxPSMscores_pep1 <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep1, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)>1){
    SingleCC_RNAandPeptides_MaxPSMscores_pep1More <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep1More, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==2){
    SingleCC_RNAandPeptides_MaxPSMscores_pep2 <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep2, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==3){
    SingleCC_RNAandPeptides_MaxPSMscores_pep3 <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep3, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)>3){
    SingleCC_RNAandPeptides_MaxPSMscores_pep3More <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep3More, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==4){
    SingleCC_RNAandPeptides_MaxPSMscores_pep4 <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep4, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==5){
    SingleCC_RNAandPeptides_MaxPSMscores_pep5 <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep5, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==6){
    SingleCC_RNAandPeptides_MaxPSMscores_pep6 <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep6, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==7){
    SingleCC_RNAandPeptides_MaxPSMscores_pep7 <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep7, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==8){
    SingleCC_RNAandPeptides_MaxPSMscores_pep8 <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep8, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==9){
    SingleCC_RNAandPeptides_MaxPSMscores_pep9 <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep9, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)==10){
    SingleCC_RNAandPeptides_MaxPSMscores_pep10 <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep10, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
  if(length(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)>10){
    SingleCC_RNAandPeptides_MaxPSMscores_pep10More <- c(SingleCC_RNAandPeptides_MaxPSMscores_pep10More, as.numeric(as.vector(SingleAllMaxPSMscores_RNAandPeptides[[i]]$max_psm_scores)))
  }
}
SingleCC_RNAandPeptides_MaxPSMscores_pep1 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep1)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep1)), as.character(as.vector(rep("1", length(SingleCC_RNAandPeptides_MaxPSMscores_pep1))))))
SingleCC_RNAandPeptides_MaxPSMscores_pep1More <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep1More)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep1More)), as.character(as.vector(rep(">1", length(SingleCC_RNAandPeptides_MaxPSMscores_pep1More))))))
SingleCC_RNAandPeptides_MaxPSMscores1HitWonderOrNot <- rbind(SingleCC_RNAandPeptides_MaxPSMscores_pep1, SingleCC_RNAandPeptides_MaxPSMscores_pep1More)
colnames(SingleCC_RNAandPeptides_MaxPSMscores1HitWonderOrNot) <- c("CC_type", "MaxPSMscorePerPeptide", "N_peptides")

SingleCC_RNAandPeptides_MaxPSMscores_pep2 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep2)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep2)), as.character(as.vector(rep("2", length(SingleCC_RNAandPeptides_MaxPSMscores_pep2))))))
SingleCC_RNAandPeptides_MaxPSMscores_pep3 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep3)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep3)), as.character(as.vector(rep("3", length(SingleCC_RNAandPeptides_MaxPSMscores_pep3))))))
SingleCC_RNAandPeptides_MaxPSMscores_pep3More <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep3More)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep3More)), as.character(as.vector(rep(">3", length(SingleCC_RNAandPeptides_MaxPSMscores_pep3More))))))
SingleCC_RNAandPeptides_MaxPSMscores1to3 <- rbind(SingleCC_RNAandPeptides_MaxPSMscores_pep1, SingleCC_RNAandPeptides_MaxPSMscores_pep2, SingleCC_RNAandPeptides_MaxPSMscores_pep3, SingleCC_RNAandPeptides_MaxPSMscores_pep3More)
colnames(SingleCC_RNAandPeptides_MaxPSMscores1to3) <- c("CC_type", "MaxPSMscorePerPeptide", "N_peptides")

SingleCC_RNAandPeptides_MaxPSMscores_pep4 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep4)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep4)), as.character(as.vector(rep("4", length(SingleCC_RNAandPeptides_MaxPSMscores_pep4))))))
SingleCC_RNAandPeptides_MaxPSMscores_pep5 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep5)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep5)), as.character(as.vector(rep("5", length(SingleCC_RNAandPeptides_MaxPSMscores_pep5))))))
SingleCC_RNAandPeptides_MaxPSMscores_pep6 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep6)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep6)), as.character(as.vector(rep("6", length(SingleCC_RNAandPeptides_MaxPSMscores_pep6))))))
SingleCC_RNAandPeptides_MaxPSMscores_pep7 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep7)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep7)), as.character(as.vector(rep("7", length(SingleCC_RNAandPeptides_MaxPSMscores_pep7))))))
SingleCC_RNAandPeptides_MaxPSMscores_pep8 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep8)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep8)), as.character(as.vector(rep("8", length(SingleCC_RNAandPeptides_MaxPSMscores_pep8))))))
SingleCC_RNAandPeptides_MaxPSMscores_pep9 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep9)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep9)), as.character(as.vector(rep("9", length(SingleCC_RNAandPeptides_MaxPSMscores_pep9))))))
SingleCC_RNAandPeptides_MaxPSMscores_pep10 <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep10)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep10)), as.character(as.vector(rep("10", length(SingleCC_RNAandPeptides_MaxPSMscores_pep10))))))
SingleCC_RNAandPeptides_MaxPSMscores_pep10More <- as.data.frame(cbind(as.character(as.vector(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_MaxPSMscores_pep10More)))), as.numeric(as.vector(SingleCC_RNAandPeptides_MaxPSMscores_pep10More)), as.character(as.vector(rep(">10", length(SingleCC_RNAandPeptides_MaxPSMscores_pep10More))))))
SingleCC_RNAandPeptides_MaxPSMscores1to10 <- rbind(SingleCC_RNAandPeptides_MaxPSMscores_pep1, SingleCC_RNAandPeptides_MaxPSMscores_pep2, SingleCC_RNAandPeptides_MaxPSMscores_pep3, SingleCC_RNAandPeptides_MaxPSMscores_pep4, SingleCC_RNAandPeptides_MaxPSMscores_pep5, SingleCC_RNAandPeptides_MaxPSMscores_pep6, SingleCC_RNAandPeptides_MaxPSMscores_pep7, SingleCC_RNAandPeptides_MaxPSMscores_pep8, SingleCC_RNAandPeptides_MaxPSMscores_pep9, SingleCC_RNAandPeptides_MaxPSMscores_pep10, SingleCC_RNAandPeptides_MaxPSMscores_pep10More)
colnames(SingleCC_RNAandPeptides_MaxPSMscores1to10) <- c("CC_type", "MaxPSMscorePerPeptide", "N_peptides")

BoxplotAllPeptidesMaxPSMscore1HitWonderOrNot_SingleCCs <- rbind(SingleCC_NoRNAsubseq_MaxPSMscores1HitWonderOrNot, SingleCC_RNAandPeptides_MaxPSMscores1HitWonderOrNot)
BoxplotAllPeptidesMaxPSMscore1HitWonderOrNot_SingleCCs$N_peptides <-factor(BoxplotAllPeptidesMaxPSMscore1HitWonderOrNot_SingleCCs$N_peptides, levels=c("1",">1"))
pdf("BoxplotAllPeptidesMaxPSMscore_SingleCCNoRNAsubseq_SingleCCRNAandPeptides_1HitWonderOrNot.pdf")
ggplot(BoxplotAllPeptidesMaxPSMscore1HitWonderOrNot_SingleCCs, aes(x=N_peptides, y=as.numeric(as.vector(MaxPSMscorePerPeptide)), fill=CC_type)) + 
  geom_boxplot() +
  labs(x="CC_type", y="max PSM score per peptide")
dev.off()

BoxplotAllPeptidesMaxPSMscore1to3_SingleCCs <- rbind(SingleCC_NoRNAsubseq_MaxPSMscores1to3, SingleCC_RNAandPeptides_MaxPSMscores1to3)
BoxplotAllPeptidesMaxPSMscore1to3_SingleCCs$N_peptides <-factor(BoxplotAllPeptidesMaxPSMscore1to3_SingleCCs$N_peptides, levels=c("1","2",">3"))
pdf("BoxplotAllPeptidesMaxPSMscore_SingleCCNoRNAsubseq_SingleCCRNAandPeptides_ByNb1to3Peptides.pdf")
ggplot(BoxplotAllPeptidesMaxPSMscore1to3_SingleCCs, aes(x=N_peptides, y=as.numeric(as.vector(MaxPSMscorePerPeptide)), fill=CC_type)) + 
  geom_boxplot() +
  labs(x="CC_type", y="max PSM score per peptide")
dev.off()

BoxplotAllPeptidesMaxPSMscore1to10_SingleCCs <- rbind(SingleCC_NoRNAsubseq_MaxPSMscores1to10, SingleCC_RNAandPeptides_MaxPSMscores1to10)
BoxplotAllPeptidesMaxPSMscore1to10_SingleCCs$N_peptides <-factor(BoxplotAllPeptidesMaxPSMscore1to10_SingleCCs$N_peptides, levels=c("1","2","3","4","5","6","7","8","9","10",">10"))
pdf("BoxplotAllPeptidesMaxPSMscore_SingleCCNoRNAsubseq_SingleCCRNAandPeptides_ByNb1to10Peptides.pdf")
ggplot(BoxplotAllPeptidesMaxPSMscore1to10_SingleCCs, aes(x=N_peptides, y=as.numeric(as.vector(MaxPSMscorePerPeptide)), fill=CC_type)) + 
  geom_boxplot() +
  labs(x="CC_type", y="max PSM score per peptide")
dev.off()


##### 10.4 COMPARE NB PEPTIDES CCs WITH NO RNA_SUBSEQS VS CCs WITH BOTH RNA_SUBSEQs AND PEPTIDES
### 10.4.1 Compare number of peptides per connected component
title <- paste0("Nb SingleCC_NoRNAsubseq=", length(unique(SingleCC_NoRNAsubseq_details$ConnCompID)), "; Nb SingleCC_RNAandPeptides=",length(unique(SingleCC_RNAandPeptides_details$ConnCompID)))
subtitle <- paste0("Nb peptides SingleCC_NoRNAsubseq=", sum(SingleCC_NoRNAsubseq_details$N_peptides), "; Nb peptides SingleCC_RNAandPeptides=",sum(SingleCC_RNAandPeptides_details$N_peptides))
BoxplotPeptidesSingleCC <- as.data.frame(rbind(cbind(rep("SingleCC_NoRNAsubseq", length(SingleCC_NoRNAsubseq_details$N_peptides)), SingleCC_NoRNAsubseq_details$N_peptides), cbind(rep("SingleCC_RNAandPeptides", length(SingleCC_RNAandPeptides_details$N_peptides)), SingleCC_RNAandPeptides_details$N_peptides)))
colnames(BoxplotPeptidesSingleCC) <-c("CC_type", "N_peptides")
pdf("BoxplotNbPeptides_SingleCCNoRNAsubseq_SingleCCRNAandPeptides.pdf")
boxplot(as.numeric(as.vector(BoxplotPeptidesSingleCC$N_peptides))~BoxplotPeptidesSingleCC$CC_type, data=BoxplotPeptidesSingleCC, xlab="CC type", ylab="Nb peptides (per CC)", main=title, cex.main=1, sub=subtitle, cex.sub=0.8)
dev.off()
pdf("BoxplotNbPeptides_SingleCCNoRNAsubseq_SingleCCRNAandPeptides_log10.pdf")
boxplot(log10(as.numeric(as.vector(BoxplotPeptidesSingleCC$N_peptides)))~BoxplotPeptidesSingleCC$CC_type, data=BoxplotPeptidesSingleCC, xlab="CC type", ylab="Log10 nb peptides (per CC)", main=title, cex.main=1, sub=subtitle, cex.sub=0.8)
dev.off()

### 10.4.2 Compare nb singleCCs with 1-to-10 or >10 peptides 
Npeptides_SingleCC_NoRNAsubseq <- vector()
for(i in 1:10){
  Npeptides_SingleCC_NoRNAsubseq <- c(Npeptides_SingleCC_NoRNAsubseq, dim(CC_NoRNAsubseq[((CC_NoRNAsubseq$N_transcripts==1)&(CC_NoRNAsubseq$N_peptides==i)),])[1])
}
Npeptides_SingleCC_NoRNAsubseq <- c(Npeptides_SingleCC_NoRNAsubseq, dim(CC_NoRNAsubseq[((CC_NoRNAsubseq$N_transcripts==1)&(CC_NoRNAsubseq$N_peptides>10)),])[1])
Npeptides_SingleCC_NoRNAsubseq <- as.data.frame(cbind(c("1","2","3","4","5","6","7","8","9","10",">10"), Npeptides_SingleCC_NoRNAsubseq))
colnames(Npeptides_SingleCC_NoRNAsubseq) <- c("N_peptides", "N_CCs")
write.table(file=paste0(path, "/Npeptides_SingleCC_NoRNAsubseq.txt"), Npeptides_SingleCC_NoRNAsubseq, sep="\t", col.names = T, row.names = F, quote=F)

Npeptides_SingleCC_RNAandPeptides <- vector()
for(i in 1:10){
  Npeptides_SingleCC_RNAandPeptides <- c(Npeptides_SingleCC_RNAandPeptides, dim(CC_RNAandPeptides[((CC_RNAandPeptides$N_transcripts==1)&(CC_RNAandPeptides$N_peptides==i)),])[1])
}
Npeptides_SingleCC_RNAandPeptides <- c(Npeptides_SingleCC_RNAandPeptides, dim(CC_RNAandPeptides[((CC_RNAandPeptides$N_transcripts==1)&(CC_RNAandPeptides$N_peptides>10)),])[1])
Npeptides_SingleCC_RNAandPeptides <- as.data.frame(cbind(c("1","2","3","4","5","6","7","8","9","10",">10"), Npeptides_SingleCC_RNAandPeptides))
colnames(Npeptides_SingleCC_RNAandPeptides) <- c("N_peptides", "N_CCs")
write.table(file=paste0(path, "/Npeptides_SingleCC_RNAandPeptides.txt"), Npeptides_SingleCC_RNAandPeptides, sep="\t", col.names = T, row.names = F, quote=F)

### 10.4.3 Compare nb peptides in single CCs with or without RNA: one-hit wonders?
# Only single-transcript CCs
Npeptides_CC_NoRNAsubseq_details_Single <- as.data.frame(table(SingleCC_NoRNAsubseq_details$N_peptides))
colnames(Npeptides_CC_NoRNAsubseq_details_Single) <- c("N_peptides", "N_CCs")
write.table(file=paste0(path, "/Npeptides_CC_NoRNAsubseq_details_Single.txt"), Npeptides_CC_NoRNAsubseq_details_Single, col.names = T, sep="\t")
pdf("Npeptides_CC_NoRNAsubseq_details_Single.pdf")
plot(as.numeric(as.vector(Npeptides_CC_NoRNAsubseq_details_Single$N_peptides)), as.numeric(as.vector(Npeptides_CC_NoRNAsubseq_details_Single$N_CCs)), type='h', xlab ="N peptides", ylab="N CCs", main="Single-transcript CC with no RNAsubseqs")
dev.off()
Npeptides_CC_PeptidesAndRNAsubseq_details_Single <- as.data.frame(table(SingleCC_RNAandPeptides_details$N_peptides))
colnames(Npeptides_CC_PeptidesAndRNAsubseq_details_Single) <- c("N_peptides", "N_CCs")
write.table(file=paste0(path, "/Npeptides_CC_PeptidesAndRNAsubseq_details_Single.txt"), Npeptides_CC_PeptidesAndRNAsubseq_details_Single, col.names = T, sep="\t")
pdf("Npeptides_CC_PeptidesAndRNAsubseq_details_Single.pdf")
plot(as.numeric(as.vector(Npeptides_CC_PeptidesAndRNAsubseq_details_Single$N_peptides)), as.numeric(as.vector(Npeptides_CC_PeptidesAndRNAsubseq_details_Single$N_CCs)), type='h', xlab ="N peptides", ylab="N CCs", main="Single-transcript CC with RNAsubseq+peptides", xlim=c(0,120), ylim=c(0,350))
dev.off()
# Multi-transcript CCs
OneHitWonder_NoRNAsubseq_peptides <- vector()
OneHitWonder_NoRNAsubseq <- as.data.frame(matrix(nrow = length(CC_NoRNAsubseq_details$ConnCompID), ncol = 3))
colnames(OneHitWonder_NoRNAsubseq) <- c("ConnCompID", "N_OneHitWonders", "Tot_Ntranscripts")
for(cc in 1:length(CC_NoRNAsubseq_details$ConnCompID)){
  ConnCompID=CC_NoRNAsubseq_details[cc,]$ConnCompID
  # if single-transcript CC, if only one peptide it is one-hit wonder
  if((CC_NoRNAsubseq_details[cc,]$N_transcripts==1)&(CC_NoRNAsubseq_details[cc,]$N_peptides==1)){
    OneHitWonder_NoRNAsubseq[cc,]$ConnCompID <- ConnCompID
    OneHitWonder_NoRNAsubseq[cc,]$N_OneHitWonders <- 1
    OneHitWonder_NoRNAsubseq[cc,]$Tot_Ntranscripts <- 1
    OneHitWonder_NoRNAsubseq_peptides <- c(OneHitWonder_NoRNAsubseq_peptides, all.cc.subseq[[ConnCompID]])
  }
  # if multi-transcript CC extract one-hit wonders using original PROT incidence matrix
  else{
    OneHitWonders <- 0
    for(t in 1:length(all.cc.transcripts[[ConnCompID]])){
      transcript=all.cc.transcripts[[ConnCompID]][t]
      index_col=which(colnames(incidenceMatrix_PSM)==transcript)
      if(sum(incidenceMatrix_PSM[,index_col])==1){
        OneHitWonders <- OneHitWonders+1
        index_row <- which(incidenceMatrix_PSM[,index_col]==1)
        OneHitWonder_NoRNAsubseq_peptides <- c(OneHitWonder_NoRNAsubseq_peptides, rownames(incidenceMatrix_PSM)[index_row])
      }
    }
    OneHitWonder_NoRNAsubseq[cc,]$ConnCompID <- ConnCompID
    OneHitWonder_NoRNAsubseq[cc,]$N_OneHitWonders <- OneHitWonders
    OneHitWonder_NoRNAsubseq[cc,]$Tot_Ntranscripts <- CC_NoRNAsubseq_details[cc,]$N_transcripts
  }
}


OneHitWonder_RNAandPeptides_peptides <- vector()
OneHitWonder_RNAandPeptides <- as.data.frame(matrix(nrow = length(CC_RNAandPeptides_details$ConnCompID), ncol = 3))
colnames(OneHitWonder_RNAandPeptides) <- c("ConnCompID", "N_OneHitWonders", "Tot_Ntranscripts")
for(cc in 1:length(CC_RNAandPeptides_details$ConnCompID)){
  ConnCompID=CC_RNAandPeptides_details[cc,]$ConnCompID
  # if single-transcript CC, if only one peptide it is one-hit wonder
  if((CC_RNAandPeptides_details[cc,]$N_transcripts==1)&(CC_RNAandPeptides_details[cc,]$N_peptides==1)){
    OneHitWonder_RNAandPeptides[cc,]$ConnCompID <- ConnCompID
    OneHitWonder_RNAandPeptides[cc,]$N_OneHitWonders <- 1
    OneHitWonder_RNAandPeptides[cc,]$Tot_Ntranscripts <- 1
    OneHitWonder_RNAandPeptides_peptides <- c(OneHitWonder_RNAandPeptides_peptides, all.cc.subseq[[ConnCompID]][grep("P_", all.cc.subseq[[ConnCompID]])])
  }
  # if multi-transcript CC extract one-hit wonders using original PROT incidence matrix
  else{
    OneHitWonders <- 0
    # Extract transcripts with peptide hits (in case of multiCCs there may be some transcripts with only RNA subseqs)
    transcripts <- all.cc.transcripts[[ConnCompID]][all.cc.transcripts[[ConnCompID]] %in% colnames(incidenceMatrix_PSM)]
    for(t in 1:length(transcripts)){
      transcript=transcripts[t]
      index_col=which(colnames(incidenceMatrix_PSM)==transcript)
      if(sum(incidenceMatrix_PSM[,index_col])==1){
        OneHitWonders <- OneHitWonders+1
        index_row <- which(incidenceMatrix_PSM[,index_col]==1)
        OneHitWonder_RNAandPeptides_peptides <- c(OneHitWonder_RNAandPeptides_peptides, rownames(incidenceMatrix_PSM)[index_row])
      }
    }
    OneHitWonder_RNAandPeptides[cc,]$ConnCompID <- ConnCompID
    OneHitWonder_RNAandPeptides[cc,]$N_OneHitWonders <- OneHitWonders
    OneHitWonder_RNAandPeptides[cc,]$Tot_Ntranscripts <- length(transcripts) # I count only transcripts with at least one peptide
  }
}


sum(OneHitWonder_NoRNAsubseq$N_OneHitWonders)/sum(OneHitWonder_NoRNAsubseq$Tot_Ntranscripts)*100  # 14.22%
sum(OneHitWonder_RNAandPeptides$N_OneHitWonders)/sum(OneHitWonder_RNAandPeptides$Tot_Ntranscripts)*100 # 6.38%

PSMscore <- read.table(paste0(path, "/AllValidatedPSMsFromProteinSets_Ensembl_03122019_ContamHitReadableByR.txt"), sep="\t", header=T, quote="")
PSMscore$sequence <- as.character(as.vector(PSMscore$sequence))
RowNamesPROTSimplified <- as.data.frame(read.table(file=paste0(path, "/RowNamesPROTSimplified.txt"), sep="\t", header=T))
OneHitWonder_NoRNAsubseq_PSMscores <- vector()
OneHitWonder_NoRNAsubseq_MaxPSMscores <- vector()
for(p in 1:length(OneHitWonder_NoRNAsubseq_peptides)){
  peptide_seq <- RowNamesPROTSimplified[RowNamesPROTSimplified$peptide_id==OneHitWonder_NoRNAsubseq_peptides[p],]$peptide_seq
  score <- PSMscore[PSMscore$sequence==peptide_seq,]$psm_score
  maxscore <- max(score)
  OneHitWonder_NoRNAsubseq_PSMscores <- c(OneHitWonder_NoRNAsubseq_PSMscores, score)
  OneHitWonder_NoRNAsubseq_MaxPSMscores <- c(OneHitWonder_NoRNAsubseq_MaxPSMscores, maxscore)
}
OneHitWonder_RNAandPeptides_PSMscores <- vector()
OneHitWonder_RNAandPeptides_MaxPSMscores <- vector()
for(p in 1:length(OneHitWonder_RNAandPeptides_peptides)){
  peptide_seq <- RowNamesPROTSimplified[RowNamesPROTSimplified$peptide_id==OneHitWonder_RNAandPeptides_peptides[p],]$peptide_seq
  score <- PSMscore[PSMscore$sequence==peptide_seq,]$psm_score
  maxscore <- max(score)
  OneHitWonder_RNAandPeptides_PSMscores <- c(OneHitWonder_RNAandPeptides_PSMscores, score)
  OneHitWonder_RNAandPeptides_MaxPSMscores <- c(OneHitWonder_RNAandPeptides_MaxPSMscores, maxscore)
}

OneHitWonder_NoRNAsubseq_PSMscores <- as.data.frame(cbind(rep("CC_NoRNAsubseq", length(OneHitWonder_NoRNAsubseq_PSMscores)), as.numeric(as.vector(OneHitWonder_NoRNAsubseq_PSMscores))))
colnames(OneHitWonder_NoRNAsubseq_PSMscores) <- c("CC_type", "OneHitWonder_PeptideScores")
OneHitWonder_RNAandPeptides_PSMscores <-as.data.frame(cbind(rep("CC_RNAandPeptides", length(OneHitWonder_RNAandPeptides_PSMscores)), as.numeric(as.vector(OneHitWonder_RNAandPeptides_PSMscores))))
colnames(OneHitWonder_RNAandPeptides_PSMscores) <- c("CC_type", "OneHitWonder_PeptideScores")
Boxplot_OneHitWonderScores <- rbind(OneHitWonder_NoRNAsubseq_PSMscores, OneHitWonder_RNAandPeptides_PSMscores)
pdf(paste0(path, "/Boxplot_OneHitWonderScores_CCNoRNAsubseq_vs_CCRNAandPeptides.pdf"))
boxplot(as.numeric(as.vector(Boxplot_OneHitWonderScores$OneHitWonder_PeptideScores))~as.factor(Boxplot_OneHitWonderScores$CC_type), xlab="CC_type", ylab="One-hit wonder peptide scores")
dev.off()

OneHitWonder_NoRNAsubseq_MaxPSMscores <- as.data.frame(cbind(rep("CC_NoRNAsubseq", length(OneHitWonder_NoRNAsubseq_MaxPSMscores)), as.numeric(as.vector(OneHitWonder_NoRNAsubseq_MaxPSMscores))))
colnames(OneHitWonder_NoRNAsubseq_MaxPSMscores) <- c("CC_type", "OneHitWonder_PeptideScores")
OneHitWonder_RNAandPeptides_MaxPSMscores <-as.data.frame(cbind(rep("CC_RNAandPeptides", length(OneHitWonder_RNAandPeptides_MaxPSMscores)), as.numeric(as.vector(OneHitWonder_RNAandPeptides_MaxPSMscores))))
colnames(OneHitWonder_RNAandPeptides_MaxPSMscores) <- c("CC_type", "OneHitWonder_PeptideScores")
Boxplot_OneHitWonderMaxScores <- rbind(OneHitWonder_NoRNAsubseq_MaxPSMscores, OneHitWonder_RNAandPeptides_MaxPSMscores)
pdf(paste0(path, "/Boxplot_OneHitWonderBestScorePerPeptide_CCNoRNAsubseq_vs_CCRNAandPeptides.pdf"))
boxplot(as.numeric(as.vector(Boxplot_OneHitWonderMaxScores$OneHitWonder_PeptideScores))~as.factor(Boxplot_OneHitWonderMaxScores$CC_type), xlab="CC_type", ylab="One-hit wonder scores (best PSM score per peptide)")
dev.off()

##### 10.5 DEPRECATED: SINCE READS WERE MERGED INTO RNA_SUBSEQUENCES THEIR NUMBER IS UNINFORMATIVE, COVERAGE INSTEAD WOULD BE
#####  SINGLE TRANSCRIPT CCs HAVE MORE RNA SUBSEQUENCES AND PEPTIDES THAN MULTI-TRANSCRIPT CCs?     ######################  
# Define function to plot, for each CC, nb RNA_subsequences vs nb peptides and indicate nb transcripts as point size
plot_nbRNA_nbPeptide_nbTranscript <- function(df_all, df_singleTrCC, df_multiTrCC, CCs_toUse){
  pdf(paste0("Scatterplot_Nsubseqspeptides_Ntranscriptsize_allCCs", CCs_toUse, ".pdf"), width = 18, height = 12)
  ggplot(df_all, aes(x=N_RNAsubsequences, y=N_peptides)) + geom_point(aes(size=N_transcripts)) +
    theme_bw() + xlim(0, 500) + ylim(0,700) + facet_wrap(~ CCtype) + scale_size(range = c(1,10), breaks=c(1,2,10,50,100,200)) +
    theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), strip.text = element_text(size = 20))
  dev.off()
  pdf(paste0("Scatterplot_Nsubseqspeptides_Ntranscriptsize_allCCs", CCs_toUse, "_log10.pdf"), width = 18, height = 12)
  ggplot(df_all, aes(x= log10(as.numeric(as.vector(df_all$N_RNAsubsequences))+1), y=log10(as.numeric(as.vector(df_all$N_peptides))+1))) + geom_point(aes(size=N_transcripts)) +
    theme_bw() + facet_wrap(~ CCtype) + scale_size(range = c(1,10), breaks=c(1,2,10,50,100,200)) +
    labs(x="log10 (# RNAsubseqs +1)", y="log10 (# peptides +1)") +
    theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), strip.text = element_text(size = 20))
  dev.off()
  
  df_singleTrCC_melt <- data.table::melt(data=as.data.table(df_singleTrCC), id.vars = c("ConnCompID"), measure.vars = c("N_peptides","N_RNAsubsequences"))
  df_singleTrCC_melt <- cbind(as.data.table(df_singleTrCC_melt), rep("singleTransc_CC",dim(df_singleTrCC_melt)[1]))
  colnames(df_singleTrCC_melt) <- c("ConnCompID", "SubseqsOrPeptides", "N", "CCtype")
  df_multiTrCC_melt <- data.table::melt(data=as.data.table(df_multiTrCC), id.vars = c("ConnCompID"), measure.vars = c("N_peptides","N_RNAsubsequences"))
  df_multiTrCC_melt <- cbind(df_multiTrCC_melt, rep("multiTransc_CC",dim(df_multiTrCC_melt)[1]))
  colnames(df_multiTrCC_melt) <- c("ConnCompID", "SubseqsOrPeptides", "N", "CCtype")
  Boxplot_Nsubseqspeptides_singlemultiTrCCs <- rbind(df_singleTrCC_melt, df_multiTrCC_melt)
  pdf(paste0(path, "/Boxplot_Nsubseqspeptides_singlemultiTrCCs", CCs_toUse, ".pdf"))
  ggplot(Boxplot_Nsubseqspeptides_singlemultiTrCCs, aes(x=CCtype, y=N, fill=SubseqsOrPeptides)) + geom_boxplot() + theme_bw()
  dev.off()
  pdf(paste0(path, "/Boxplot_Nsubseqspeptides_singlemultiTrCCs", CCs_toUse, "_Log10.pdf"))
  ggplot(Boxplot_Nsubseqspeptides_singlemultiTrCCs, aes(x=CCtype, y=log10(N+1), fill=SubseqsOrPeptides)) + geom_boxplot() + theme_bw()
  dev.off()
}

### 10.5.1 Apply on all CCs function to plot, for each CC, nb RNA_subsequences vs nb peptides and indicate nb transcripts as point size
m_singleTrCC <- m[m$N_transcripts==1,]  # dim: 10238     4
m_multiTrCC <- m[m$N_transcripts>1,]  # dim: 11133     4
m_allCC <- rbind(m_singleTrCC, m_multiTrCC)
m_allCC$CCtype <- ifelse(m_allCC$N_transcripts==1, "singleTransc_CC", "multiTransc_CC")
m_allCC$color <- ifelse(m_allCC$CCtype=="singleTransc_CC", "Orange", "Blue")
m_allCC$color <- as.factor(m_allCC$color)
plot_nbRNA_nbPeptide_nbTranscript(m_allCC, m_singleTrCC, m_multiTrCC, "")
### 10.5.2 Apply on only CCs which contain at least one peptide function to plot, for each CC, nb RNA_subsequences vs nb peptides and indicate nb transcripts as point size
m_singleTrCC_2 <- m[((m$N_transcripts==1)&(m$N_peptides>0)),] # dim: 1113    4
m_multiTrCC_2 <- m[((m$N_transcripts>1)&(m$N_peptides>0)),] # dim: 8014     4
m_allCC_2 <- rbind(m_singleTrCC_2, m_multiTrCC_2)
m_allCC_2$CCtype <- ifelse(m_allCC_2$N_transcripts==1, "singleTransc_CC", "multiTransc_CC")
m_allCC_2$color <- ifelse(m_allCC_2$CCtype=="singleTransc_CC", "Orange", "Blue")
m_allCC_2$color <- as.factor(m_allCC_2$color)
plot_nbRNA_nbPeptide_nbTranscript(m_allCC_2, m_singleTrCC_2, m_multiTrCC_2, "ContainingAtLeastOnePeptide")
### 10.4.3 Apply on only CCs which contain at least one peptide and one RNAsubseq function to plot, for each CC, nb RNA_subsequences vs nb peptides and indicate nb transcripts as point size
m_singleTrCC_3 <- m[((m$N_transcripts==1)&(m$N_peptides>0)&(m$N_RNAsubsequences>0)),] # dm: 343     4
m_multiTrCC_3 <- m[((m$N_transcripts>1)&(m$N_peptides>0)&(m$N_RNAsubsequences>0)),] # dim: 4790     4
m_allCC_3 <- rbind(m_singleTrCC_3, m_multiTrCC_3)
m_allCC_3$CCtype <- ifelse(m_allCC_3$N_transcripts==1, "singleTransc_CC", "multiTransc_CC")
m_allCC_3$color <- ifelse(m_allCC_3$CCtype=="singleTransc_CC", "Orange", "Blue")
m_allCC_3$color <- as.factor(m_allCC_3$color)
plot_nbRNA_nbPeptide_nbTranscript(m_allCC_3, m_singleTrCC_3, m_multiTrCC_3, "ContainingAtLeastOnePeptideAndOneRNAsubseq")


####################################################################################################################
####                                  11. PLOT CONNECTED COMPONENTS                                             ####
####################################################################################################################

### Using igraph package
# Plot 10 randomly selected CCs (multitranscripts but 2-4 transcripts only to be able to visualize them; reads and peptides)
CCs_toPlot <- c(157, 356, 915, 4146, 2151, 2557, 778, 2361, 2223, 2326)
pdf(paste0(path, "/Plot_ExamplesCCsMultiTranscr_SelectedCCs.pdf"))
for (CC_id in 1:length(CCs_toPlot)){
  g <- igraph::graph_from_incidence_matrix(multTranscript.cc.incMatr[[CCs_toPlot[CC_id]]])
  V(g)$label.cex <- 0.5
  V(g)$label.color <- "black"
  V(g)$color <- rep("#0072B2", length(names(as.list(V(g))))) # Blue for RNAseq subseq nodes
  V(g)$color[grep("ENST",names(as.list(V(g))))] <- "#D55E00" # Orange for transcript nodes
  V(g)$color[grep("P_",names(as.list(V(g))))] <- "#009E73" # Green for peptide nodes
  #E(g)$color <- rep("#0072B2", length(names(as.list(E(g))))) # Blue for edges ro RNAseq subseqs
  #E(g)$color[grep("P_",names(as.list(E(g))))] <- "#009E73" # Green for peptide nodes
  #E(g)$width <- rep("2", length(names(as.list(E(g)))))
  opar <- par()$mar; par(mar=rep(0, 4)) #Give the graph lots of room
  #igraph::plot.igraph(g, edge.color="gray30", cex=.5,layout=layout_as_bipartite)
  igraph::plot.igraph(g, edge.color="gray30", cex=.5, vertex.size=5)
  #title(paste0("CC_", CC_id), line = -1)
  title(paste0("CC_", CCs_toPlot[CC_id]), line = -1)
  par(mar=opar)
}
dev.off()
