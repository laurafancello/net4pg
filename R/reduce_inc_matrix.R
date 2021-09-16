#' Reduce size of incidence matrix for downstream analyses
#'
#' Reduce the size of the incidence matrix describing peptide-to-protein
#' mappings to ease downstream analyses. The original incidence matrix is
#' reduced to only contain proteins with at least one shared peptide and all
#' peptides mapping on such proteins. This means that only proteins ambiguously
#' identified are retained, which is the most interesting ones when studying
#' ambiguity of protein identifications.
#' @param incM a \code{logical} \code{matrix} containing the incidence matrix
#' with its column and row names (respectively, protein and peptide identifiers)
#' and 0 or 1 values indicating whether or not the peptide maps on the
#' corresponding protein.
#' @return a \code{logical} \code{matrix} containing a smaller incidence matrix
#' (with column and row names respectively reporting protein and peptide
#' identifiers) and 0 or 1 values indicating whether or not the peptide maps on
#' the corresponding protein. Only proteins with at least one shared peptide and
#' all peptides mapping on such protein are reported in such reduced incidence
#' matrix.
#' @examples
#' # Read the tab-delimited file containing he proteome incidence matrix
#' incM_filename <- system.file("extdata"
#'                              , "incM_example"
#'                              , package = "net4pg"
#'                              , mustWork = TRUE)
#' rownames_filename <- system.file("extdata"
#'                                   , "peptideIDs_incM_example"
#'                                   , package = "net4pg"
#'                                   , mustWork = TRUE)
#' colnames_filename <- system.file("extdata"
#'                                  , "proteinIDs_incM_example"
#'                                  , package = "net4pg"
#'                                  , mustWork = TRUE)
#' incM <- read_inc_matrix(incM_filename = incM_filename
#'                  , colnames_filename = colnames_filename
#'                  , rownames_filename = rownames_filename)
#' # Only retain proteins with at least one shared peptide and all peptides
#' # mapping on such proteins.
#' incM_reduced <- reduce_inc_matrix(incM)
#'
#' @author Laura Fancello
#'
#' @export

reduce_inc_matrix <- function(incM) {

  # Sanity Checks  ----------------------------------------------------------
  ## Check input arguments
  if (is.null(incM)) {
    stop("argument 'incM' is missing, with no default")
  }
  if (!(methods::is(incM)[1] == "matrix")) {
    stop("argument 'incM' is not a matrix")
  }

  # Reduce incidence matrix  --------------------------------------------------
  ## First remove peptides only pointing to one protein (specific peptides)
  incM_RowFilter <- incM[-which(rowSums(incM) == 1), ]

  ## Probably useful only for smaller toy datasets where it can occur that only
  ## one peptide is left
  if (methods::is(incM_RowFilter)[2] == "vector"){
    incM_RowFilter <- t(as.matrix(incM_RowFilter))
    rownames(incM_RowFilter) <- rownames(incM)[which(rowSums(incM) > 1)]
  }

  ## Then remove proteins with 0 peptides (which is proteins only identified by
  ## specific peptides, removed on previous step)
  incM_RowColFilter <- incM_RowFilter[, -which(colSums(incM_RowFilter) == 0)]

  ## Probably useful only for smalle toy datasets where it can occur that only
  ## one peptide is left
  if (methods::is(incM_RowColFilter)[2] == "vector") {
    incM_RowColFilter <- t(as.matrix(incM_RowColFilter))
    rownames(incM_RowColFilter) <- rownames(incM_RowFilter)
  }

  ## Clean memory
  rm(incM_RowFilter)
  gc()

  ## Return output
  return(incM_RowColFilter)

}
