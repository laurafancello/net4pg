#' Get peptides and peptide-to-protein mappings for each connected component
#'
#' Get peptides and peptide-to-protein mappings for each connected component.
#' For each connected component, first extract its protein members; then,
#' extract all specific and shared peptides mapping on those protein; finally,
#' extract the subset of incidence matrix representing peptide-to-protein
#' mappings.
#' @param cc.proteins a \code{list} of \code{vectors} (one for each connected
#' component), each enumerating the proteins members of a connected component.
#' @param incM a \code{logical} \code{matrix} containing the incidence matrix
#' with its column and row names (respectively, protein and peptide identifiers)
#' names and 0 or 1 values indicating whether or not the peptide maps on the
#' corresponding protein.
#' @return a \code{list} of two elements: i. a \code{list} of \code{vectors}
#' (one for each connected component) enumerating peptides mapping on protein
#' members of each connected component; ii. a \code{list} of \code{matrices} or
#' \code{vectors} (one for each connected component) representing
#' peptide-to-protein mappings for each connected component; matrices are used
#' if multiple peptides identify protein members of that connected component,
#' vectors if only a single peptide.
#' @examples
#' # Read the tab-delimited file containing the proteome incidence matrix
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
#' # Generate adjacency matrix describing protein-to-protein mappings
#' adjM <- get_adj_matrix(incM_reduced)
#' # Generate graph of protein-to-protein connections and calculate its
#' # connected component
#' multProteinCC <- get_cc(adjM)
#' # For each connected component, extract peptides mapping on its protein
#' # members and the subset of the incidence matrix describing peptide-to-protein
#' # mapping
#' cc.peptides.incM <- cc_composition(cc.proteins = multProteinCC$cc
#'                                   , incM = incM)
#'
#' @author Laura Fancello
#'
#' @export

cc_composition <- function(cc.proteins, incM) {

  # Sanity Checks  ----------------------------------------------------------
  ## Check input arguments
  if (is.null(cc.proteins)) {
    stop("argument 'cc.proteins' is missing, with no default")
  }
  if (!(methods::is(cc.proteins)[1] == "list")) {
    stop("argument 'cc.proteins' is not a list")
  }
  if (is.null(incM)) {
    stop("argument 'incM' is missing, with no default")
  }
  if (!(methods::is(incM)[1] == "matrix")) {
    stop("argument 'incM' is not a matrix")
  }

  # Peptides and peptide-to-protein mappings per connected component  --------
  cc.peptides <- list()
  cc.subincM <- list()
  for (i in seq_along(cc.proteins)) {
    proteinlist <- cc.proteins[[i]]
    subX <- incM[, which(colnames(incM) %in% proteinlist)]
    cc.peptides[[i]] <- names(which(rowSums(subX) != 0))
    cc.subincM[[i]] <- subX[which(rowSums(subX) != 0), ]
    if (methods::is(cc.subincM[[i]])[2] == "vector") {
      cc.subincM[[i]] <- t(as.matrix(cc.subincM[[i]]))
      rownames(cc.subincM[[i]]) <- cc.peptides[[i]]
    }
  }

  ## Clean memory
  rm(proteinlist, subX, i)
  gc()

  ## Return output
  return(list(cc.peptides = cc.peptides, cc.subincM = cc.subincM))

}
