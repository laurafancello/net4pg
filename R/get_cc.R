#' Generate graph and calculate its connected components
#'
#' Build a graph of protein-to-protein connections from adjacency matrix and
#' calculate its connected components.
#' @param adjM a \code{numerical} \code{matrix} containing the adjacency matrix,
#' with value >0 or 0 indicating whether or not two proteins are identified by
#' shared peptide(s)
#' @return a \code{list} of two elements: i. a \code{graph} representing
#' protein-to-protein connections encoded by the adjacency matrix; ii. a
#' \code{list} of \code{vectors} (one for each connected component) enumerating
#' protein members of each connected component.
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
#' incM <- read_inc_matrix(incM_filename=incM_filename
#'                  , colnames_filename=colnames_filename
#'                  , rownames_filename=rownames_filename)
#' # Only retain proteins with at least one shared peptide and all peptides
#' # mapping on such proteins.
#' incM_reduced <- reduce_inc_matrix(incM)
#' # Generate adjacency matrix describing protein-to-protein mappings
#' adjM <- get_adj_matrix(incM_reduced)
#' # Generate graph of protein-to-protein connections and calculate its
#' # connected components
#' multProteinCC <- get_cc(adjM)
#'
#' @author Laura Fancello
#'
#' @export
#'
get_cc <- function(adjM) {

  # Sanity Checks  ----------------------------------------------------------
  ## Check input arguments
  if (is.null(adjM)) {
    stop("argument 'incM' is missing, with no default")
  }
  if (!(methods::is(adjM)[1] == "matrix")) {
    stop("argument 'adjM' is not a matrix")
  }

  # Calculate graph and its connected components  --------------------------
  ## Generate graph
  g <- graph::graphAM(adjM, edgemode = "undirected", values = NA)

  ## Calculate connected components
  ccs <- graph::connComp(methods::as(g, "graphNEL"))

  ## Return output
  result <- list(g = g, ccs = ccs)
  return(result)

}
