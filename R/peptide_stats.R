#' Calculate percentage of shared vs specific peptides
#'
#' Read in input the incidence matrix of peptide-to-protein mappings generated
#' from valid proteomic identifications
#' @param incM a \code{logical} \code{matrix} containing the incidence matrix
#' with its column and row names (respectively, protein and peptide identifiers)
#' and 0 or 1 values indicating whether or not the peptide maps on the
#' corresponding protein.
#' @return a \code{list} of three elements: i. number of shared peptides;
#' ii. number of specific peptides; iii. percentage of specific peptides
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
#' # Calculate percentage of shared vs specific peptides
#' peptideStatsOut <- peptide_stats(incM = incM)
#'
#' @author Laura Fancello
#'
#' @export


peptide_stats <- function(incM) {

  # Sanity Checks  ----------------------------------------------------------
  ## Check input arguments
  if (is.null(incM)) {
    stop("argument 'incM' is missing, with no default")
  }
  if (!(methods::is(incM)[1] == "matrix")) {
    stop("argument 'incM' is not a matrix")
  }

  ## Percentage of shared vs specific peptides  ------------------------------
  totPeptides <- dim(incM)[1]
  nbSpecific <- length(which(rowSums(incM) == 1))
  nbShared <- totPeptides - nbSpecific
  percSpecific <- paste0(round(nbSpecific / totPeptides * 100, digits = 2), "%")

  return(list(nbShared = nbShared
              , nbSpecific = nbSpecific
              , percSpecific = percSpecific))

}
