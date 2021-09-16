#' Plot peptide-to-protein mapping graph
#'
#' Plot a bipartite subgraph representing the connected component to which the
#' user-selected protein belongs. Peptide-to-protein mappings of that protein
#' and of all other proteins belonging to the same CC are represented. The
#' function takes in input a single Ensembl protein identifier (i.e. ENSPXXX
#' for human, ENSMUSPXXX for mouse), it identifies the connected component it
#' belongs to and plots all peptide-to-protein mappings of that connected
#' component.
#'
#' @param prot a \code{character} \code{vector} containing a single Ensembl
#' identifier (i.e. ENSPXXXXXXXXXXX for human, ENSMUSPXXXXXXXXXXX for mouse) of
#' the protein of interest.
#' @param cc.proteins a \code{list} of \code{vectors} (one for each connected
#' component) containing protein members of each connected component.
#' @param cc.subincM a \code{list} of \code{matrices} or \code{vectors} (one for
#' each connected component) representing the incidence matrix of
#' peptide-to-protein mappings for each connected component; matrices are used
#' if multiple peptides identify protein members of that connected component,
#' vectors if only a single peptide.
#' @param tagProt a \code{character} \code{vector} reporting the prefix of
#' protein identifiers (for non contaminant proteins)
#' @param tagContam a \code{character} \code{vector} reporting the tag which
#' identifies contaminant proteins
#' @param incM a \code{logical} \code{matrix} containing the incidence matrix
#' with its column and row names (respectively, protein and peptide identifiers)
#' names and 0 or 1 values indicating whether or not the peptide maps on the
#' corresponding protein.
#' @return a list of four elements: i. CC identifier; ii. protein members of the
#' CC; iii. peptide members of the CC; iv. bipartite subgraph representing
#' peptide-to-protein mapping of that CC (if multi-protein CC)
#' @examples
#' library(igraph)
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
#' # Only retain proteins with at least one shared peptide and all peptides
#' # mapping on such proteins.
#' incM_reduced <- reduce_inc_matrix(incM)
#' # Generate adjacency matrix describing protein-to-protein mappings
#' adjM <- get_adj_matrix(incM_reduced)
#' # Generate graph of protein-to-protein connections and calculate its
#' # connected components
#' multProteinCC <- get_cc(adjM)
#' # For each connected component, extract peptides mapping on its protein
#' # members and the subset of the incidence matrix describing
#' # peptide-to-protein mappings
#' cc.peptides.incM <- cc_composition(cc.proteins = multProteinCC$cc
#'                                    , incM = incM)
#' # Plot bipartite graph representing peptide-to-protein mappings for the
#' # connected component of the protein of interest (in this toy example protein
#' # "ENSP261"; note that identifiers are not authentic but made up for the
#' # example)
#' subgraphCC <- plot_cc(prot="ENSP261"
#'         , cc.proteins=multProteinCC$ccs
#'         , cc.subincM=cc.peptides.incM$cc.subincM
#'         , tagProt = "ENSP"
#'         , tagContam="Contam"
#'         , incM=incM)
#' plot.igraph(subgraphCC$g
#'       , layout = layout_as_bipartite
#'       , edge.width = 1
#'       , edge.arrow.width = 0.3
#'       , vertex.size = 10
#'       , edge.arrow.size = 0.5
#'       , vertex.size2 = 3
#'       , vertex.label.cex = 0.8
#'       , asp = 0.35
#'       , margin = -0.1) +
#' title(paste0("Protein ENSP261 in CC #", subgraphCC$cc_id), line = -1)
#'
#' @author Laura Fancello
#'
#' @export

plot_cc <- function(prot, cc.proteins, cc.subincM, tagProt, tagContam, incM) {

  # Sanity Checks  ----------------------------------------------------------
  ## Check input arguments
  if (is.null(prot)) {
    stop("argument 'prot' is missing, with no default")
  }
  if (!((methods::is(prot)[1] == "character") & (methods::is(prot)[2] == "vector"))) {
    stop("argument 'prot' is not a character vector")
  }
  if (length(prot) > 1) {
    stop("argument 'prot' longer than 1: please give in input a single protein
         identifier")
  }
  if ((length(grep(tagProt, prot)) == 0) & (length(grep(tagContam, prot)) == 0)) {
    stop("argument 'prot' is not a valid contaminant protein identifier: please
         give in input a correct Ensembl or contaminant protein identifier
         (i.e. ENSPXXX or protein with tagContam)")
  }
  if (length(grep(prot, colnames(incM))) == 0) {
    stop("argument 'prot' is not a protein identifier from the incidence matrix
         provided in input")
  }

  if (is.null(cc.proteins)) {
    stop("argument 'cc.proteins' is missing, with no default")
  }
  if (!(methods::is(cc.proteins)[1] == "list")) {
    stop("argument 'cc.proteins' is not a list")
  }

  if (is.null(cc.subincM)) {
    stop("argument 'cc.subincM' is missing, with no default")
  }
  if (!(methods::is(cc.subincM)[1] == "list")) {
    stop("argument 'cc.subincM' is not a list")
  }

  if (is.null(tagProt)) {
    stop("argument 'tagProt' is missing, with no default: please provide the
        prefix of protein identifiers (for non contaminant proteins")
  }

  if (is.null(tagContam)) {
    stop("argument 'tagContam' is missing, with no default: please provide the
         tag identifying contaminant proteins")
  }
  if (!((methods::is(tagContam)[1] == "character") & (methods::is(tagContam)[2] == "vector"))) {
    stop("argument 'tagContam' is not a character vector")
  }
  if (length(tagContam) > 1) {
    stop("argument 'tagContam' longer than 1: only one single tag allowed for
         protein contaminants")
  }

  # Plot connected components  -----------------------------------------------
  ## Find CCs containing input proteins
  res <- lapply(cc.proteins, function(cc.proteins) grep(prot, cc.proteins))
  cc_id <- which(unlist(lapply(res, function(x) length(x) > 0)))
  result <- list()

  if (length(cc_id) == 0) {
    print(paste0("Protein ", prot, " is member of a single-protein CC"))
    index_prot <- grep(prot, colnames(incM))
    peptides <- rownames(incM)[which(incM[, index_prot] == TRUE)]
    result$cc_id <- "single-prot"
    result$proteins <- prot
    result$peptides <- peptides
    edges <- data.frame(from = peptides,
                           to = rep(prot, length(peptides)))
    g <- igraph::graph_from_data_frame(edges
                                       , directed = FALSE
                                       , vertices = c(peptides, prot))
    igraph::V(g)$type <- rep(FALSE, length(names(as.list(igraph::V(g)))))
    igraph::V(g)$type[grep(tagProt, names(as.list(igraph::V(g))))] <- TRUE
    igraph::V(g)$type[grep(tagContam, names(as.list(igraph::V(g))))] <- TRUE
  }else{
    print(paste0("Protein ", prot, " is member of a multi-protein CC"))
    ## Generate bipartite graph of peptide-to-protein mappings
    g <- igraph::graph_from_incidence_matrix(cc.subincM[cc_id][[1]])
    result$cc_id <- cc_id
    proteins <- as.character(as.vector(c(names(igraph::V(g)[grep(tagProt, names(as.list(igraph::V(g))))])))
                             ,  names(igraph::V(g)[grep(tagContam, names(as.list(igraph::V(g))))]))
    result$proteins <- proteins
    result$peptides <- setdiff(names(as.list(igraph::V(g))), proteins)
  }
  ## Plot bipartite graph
  igraph::V(g)$label.cex <- 1
  igraph::V(g)$label.color <- "black"
  # Blue peptide vertices
  igraph::V(g)$color <- rep("#0072B2", length(names(as.list(igraph::V(g)))))
  # Orange protein vertices
  igraph::V(g)$color[grep(tagProt, names(as.list(igraph::V(g))))] <- "#D55E00"
  # Orange contaminant protein vertices
  igraph::V(g)$color[grep(tagContam, names(as.list(igraph::V(g))))] <- "#D55E00"
  # Output
  result$g <- g
  return(result)

}
