#' Read incidence matrix of proteomic identifications
#'
#' Read in input the incidence matrix of peptide-to-protein mappings generated
#' from valid proteomic identifications.
#' @param incM_filename the name of the tab-delimited file containing incidence
#' matrix values; the input incidence matrix must contain along the columns
#' protein identifiers and along the rows peptides; each cell must contain a 1
#' or 0 value indicating whether or not the peptide maps on the corresponding
#' protein.
#' @param colnames_filename name of the file containing incidence matrix column
#' names, which are protein identifiers. The file must contain one identifier
#' per line. Protein identifiers must be in Ensembl format (i.e.,
#' ENSPXXXXXXXXXXX for human).
#' @param rownames_filename name of the file containing incidence matrix row
#' names, which are peptide identifiers. The file must contain one identifier
#' per line. Peptide identifiers can be in any format (i.e. numeric identifiers,
#' amino acid sequences, ...)
#' @return a \code{logical} \code{matrix} containing the incidence matrix with
#' its column and row names (respectively, protein and peptide identifiers)
#' and 0 or 1 values indicating whether or not the peptide maps on the
#' corresponding protein.
#' @examples
#' # Read the tab-delimited file containing he proteome incidence matrix
#' incM_filename <- system.file( "extdata"
#'                              , "incM_example"
#'                              , package = "net4pg"
#'                              , mustWork = TRUE)
#' rownames_filename <- system.file( "extdata"
#'                                   , "peptideIDs_incM_example"
#'                                   , package = "net4pg"
#'                                   , mustWork = TRUE)
#' colnames_filename <- system.file( "extdata"
#'                                  , "proteinIDs_incM_example"
#'                                  , package = "net4pg"
#'                                  , mustWork = TRUE)
#' incM <- read_inc_matrix(incM_filename = incM_filename
#'                  , colnames_filename = colnames_filename
#'                  , rownames_filename = rownames_filename)
#'
#' @author Laura Fancello
#'
#' @export

read_inc_matrix <- function(incM_filename, colnames_filename, rownames_filename) {

  # Sanity Checks  ----------------------------------------------------------
  ## Check input arguments
  if (is.null(incM_filename)) {
    stop("argument 'incM_filename' is missing, with no default")
  }
  if (!((methods::is(incM_filename)[1] == "character") | (methods::is(incM_filename)[2] == "vector"))) {
    stop("argument 'incM_filename' is not a character vector")
  }
  if (is.null(colnames_filename)) {
    stop("argument 'colnames_filename' is missing, with no default")
  }
  if (!((methods::is(colnames_filename)[1] == "character") | (methods::is(colnames_filename)[2] == "vector"))) {
    stop("argument 'colnames_filename' is not a character vector")
  }
  if (is.null(rownames_filename)) {
    stop("argument 'rownames_filename' is missing, with no default")
  }
  if (!((methods::is(rownames_filename)[1] == "character") | (methods::is(rownames_filename)[2] == "vector"))) {
    stop("argument 'rownames_filename' is not a character vector")
  }

  # Read incidence matrix  --------------------------------------------------
  ## Read column names (protein IDs)
  proteinIDs <- data.table::fread(file = colnames_filename
                                  , sep = "\n"
                                  , header = FALSE
                                  , data.table = FALSE)
  proteinIDs <- proteinIDs$V1

  ## Read row names (peptide IDs)
  peptideIDs <- data.table::fread(file = rownames_filename
                                  , sep = "\n"
                                  , header = FALSE
                                  , data.table = FALSE)
  peptideIDs <- peptideIDs$V1

  ## Read incidence matrix values
  nbLines <- length(peptideIDs) # nb rows of incidence matrix

  ## If big matrix (more than 1000 rows) read it chunk by chunk
  if (nbLines > 10000) {

    chunk <- 10000
    skip <- 0
    cycles <- ceiling(nbLines / chunk) # set the number of chunks to read matrix

    for (i in 1:cycles) {

      if (i == cycles) {
        chunk <- nbLines - (chunk * (cycles - 1))
      }

      print(paste0("cycle: ", i))
      print(paste0("skip: ", skip))
      print(paste0("chunk: ", chunk))

      # Read sub-matrix (chunk)
      incM_tmp <- data.table::fread(file = incM_filename
                        , sep = "\t"
                        , header = FALSE
                        , data.table = FALSE
                        , skip = skip
                        , nrows = chunk)
      assign(paste("incM", skip, sep = ""), (incM_tmp))
      skip <- skip + chunk

    }

    rm(incM_tmp) # clean memory
    gc()

    ## Put together all chunks in a unique matrix object
    incM <- incM0
    rm(incM0)
    gc()
    incM <- as.matrix(incM)
    ## Convert to logical to save memory space
    incM <- apply(incM, 2, as.logical)

    chunk <- 10000
    for (i in 1:(cycles - 1)) {
      incM_toadd <- get(paste0("incM", (i * chunk)))
      rm(list = paste0("incM", (i * chunk)))
      incM_toadd <- as.matrix(incM_toadd)
      ## Convert to logical to save memory space
      incM_toadd <- apply(incM_toadd, 2, as.logical)
      incM <- rbind(incM, incM_toadd)
    }

    ## Clean memory
    rm(incM_toadd, cycles, chunk, skip, i, nbLines)
    gc()

  }else{  ## If not big matrix read it all at once

    ## Read incidence matrix
    incM <- data.table::fread(file = incM_filename
                              , sep = "\t"
                              , header = FALSE
                              , data.table = FALSE)
    incM <- as.matrix(incM)
    incM <- apply(incM, 2, as.logical)
  }

  ## Assign row and col names
  rownames(incM) <- peptideIDs
  colnames(incM) <- proteinIDs

  ## Return output
  return(incM)

}
