#' Perform transcriptome-informed post-hoc filtering
#'
#' Implement the transcriptome-informed post-hoc filtering strategy. This
#' strategy aims to reduce the ambiguity of protein identifications by
#' exploiting sample-matched transcriptome information, when available. First,
#' it takes in input the set of transcripts expressed in the sample-matched
#' transcriptome (reported using the transcript identifier in Ensembl format,
#' i.e., ENSTXXXX for human) and removes from proteomic identifications:
#' i. all proteins with no expressed transcripts and peptides exclusively
#' mapping on removed proteins ("all"); or
#' ii. only those exclusively identified
#' by shared peptides and peptides exclusively mapping on removed proteins
#' ("sharedOnly"); or
#' iii. only those exclusively identified by shared peptides,
#' whose peptides are shared with at least one protein with expressed
#' transcript, so they are not to be removed ("sharedNoRemove")
#' @param incM a \code{logical} \code{matrix} containing the incidence matrix
#' with its column and row names (respectively, protein and peptide identifiers)
#' and 0 or 1 values indicating whether or not the peptide maps on the
#' corresponding protein.
#' @param exprTranscriptsFile the name of the file containing the set of
#' transcripts expressed in the sample-matched transcriptome (one per line).
#' Transcript identifiers must be in the Ensembl format (i.e., ENSTXXXXXXXXXXX
#'  for human)
#' @param proteinToTranscriptFile the name of a tab-delimited file with protein
#' identifiers in the first column and the corresponding transcript identifiers
#' in the second column. Protein and transcript identifiers must be in the
#' Ensembl format (i.e. ENSPXXXXXXXXXXX and ENSTXXXXXXXXXXX for human)
#' @param tagContam a \code{character} \code{vector} reporting the tag which
#' identifies contaminant protein
#' @param remove \code{character} \code{vector} indicating whether to remove:
#' i. all proteins with no expressed transcripts and peptides exclusively
#' mapping on removed proteins ("all"); ii. only those exclusively identified
#' by shared peptides and peptides exclusively mapping on removed proteins
#' ("sharedOnly"); iii. only those exclusively identified by shared peptides,
#' whose peptides are shared with at least one protein with expressed
#' transcript, so they are not to be removed ("sharedNoRemove")
#' @importFrom magrittr %>%
#' @return a \code{matrix} representing a filtered incidence matrix of
#' peptide-to-protein mapping obtained by transcriptome-informed filtering.
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
#' # Perform transcriptome-informed post-hoc filtering
#' exprTranscriptsFile <- system.file("extdata"
#'                                    , "expressed_transcripts.txt"
#'                                    , package = "net4pg"
#'                                    , mustWork = TRUE)
#' protein2transcriptFile <- system.file("extdata"
#'                                         , "protein_to_transcript"
#'                                         , package = "net4pg"
#'                                         , mustWork = TRUE)
#' incM_filtered <- transcriptome_filter(incM
#'                          , exprTranscriptsFile = exprTranscriptsFile
#'                          , proteinToTranscriptFile = protein2transcriptFile
#'                          , tagContam = "Contam"
#'                          , remove = "all")
#'
#' @author Laura Fancello
#'
#' @export
#'

transcriptome_filter <- function(incM
                          , exprTranscriptsFile
                          , proteinToTranscriptFile
                          , tagContam
                          , remove) {

  # Sanity Checks  ----------------------------------------------------------
  ## Check input arguments
  if (is.null(incM)) {
    stop("argument 'incM' is missing, with no default")
  }
  if (!(methods::is(incM)[1] == "matrix")) {
    stop("argument 'incM' is not a matrix")
  }
  if (is.null(exprTranscriptsFile)) {
    stop("argument 'exprTranscriptsFile' is missing, with no default")
  }
  if (!((methods::is(exprTranscriptsFile)[1] == "character") | (methods::is(exprTranscriptsFile)[2] == "vector"))) {
    stop("argument 'exprTranscriptsFile' is not a character vector")
  }
  if (is.null(proteinToTranscriptFile)) {
    stop("argument 'proteinToTranscriptFile' is missing, with no default")
  }
  if (!((methods::is(proteinToTranscriptFile)[1] == "character") | (methods::is(proteinToTranscriptFile)[2] == "vector"))) {
    stop("argument 'proteinToTranscriptFile' is not a character vector")
  }
  if (is.null(tagContam)) {
    stop("argument 'tagContam' is missing, with no default")
  }
  if ((!(methods::is(tagContam)[1] == "character") | (!(methods::is(tagContam)[2] == "vector")))) {
    stop("argument 'tagContam' is not a character vector")
  }
  if (is.null(remove)) {
    stop("argument 'remove' is missing, with no default")
  }
  if ((!(methods::is(remove)[1] == "character") | (!(methods::is(remove)[2] == "vector")))) {
    stop("argument 'remove' is not a character vector")
  }
  if ((remove != "all") & (remove != "sharedOnly") & (remove != "sharedNoRemove")) {
    stop("argument 'remove' is not valid: please choose one between 'all',
         'sharedOnly' and 'sharedNoRemove'")
  }

  # Post-hoc filter  --------------------------------------------------------
  ## Read list of transcripts found to be expressed in the sample-matched transcriptome
  exprRNA <- scan(file = exprTranscriptsFile, what = character())

  ## Read tab-delimited file containing Ensembl transcript ID to Ensembl protein ID
  ## conversion
  trans2Prot <- utils::read.table(file = proteinToTranscriptFile
                                           , sep = "\t"
                                           , header = FALSE)
  colnames(trans2Prot) <- c("Prot", "RNA")

  ## Convert IDs of expressed transcript into the corresponding protein IDs
  exprProts <- as.character(as.vector(trans2Prot[trans2Prot$RNA %in% exprRNA, ]$Prot))

  ## Identify contaminant proteins
  proteinContam <- colnames(incM)[grep("Contam", colnames(incM))]

  if (remove == "sharedOnly") {
    ## Extract specific peptides
    specificPep <- which(rowSums(incM) == 1)
    ## Extract proteins with specific peptides
    subIncM <- incM[specificPep, ]
    specificProt <- colnames(subIncM[, which(colSums(subIncM) > 0)])
    onlySharedProt <- setdiff(colnames(incM), specificProt)

    ## Extract proteins with not expressed transcript (excluding contaminant proteins
    ## for which we do not have transcriptome information)
    noExprProts <- setdiff(colnames(incM), c(exprProts, proteinContam))

    ## Remove proteins with both the following features: 1. with not expressed
    ## transcript (excluding contaminant proteins for which we do not have
    ## transcriptome information) AND 2. with no specific peptide
    noKeepProt <- intersect(noExprProts, onlySharedProt)
    incM_filtered <- incM[, -which(colnames(incM) %in% noKeepProt)]

    ## Remove peptides only mapping on removed proteins
    filterPeptides_index <- which(rowSums(incM_filtered) == 0)
    if (length(filterPeptides_index) > 0){
      incM_filtered <- incM_filtered[-filterPeptides_index, ]
    }

    ## Clean memory
    rm(specificPep, subIncM, specificProt, onlySharedProt, noExprProts
       , noKeepProt, filterPeptides_index, exprRNA, exprProts)
    gc()

  }else{
    if (remove == "all") {
      ## Keep only proteins with expressed transcript in sample-matched
      ## transcriptome and contaminant proteins (which are not observed in
      ## transcriptomics when mapping against the reference genome)
      keepProt <- c(intersect(exprProts, colnames(incM)), proteinContam)
      incM_filtered <- incM[, which(colnames(incM) %in% keepProt)]

      ## Remove peptides only mapping on proteins whose transcript is NOT
      ## expressed in the sample-matched transcriptome
      filterPeptides_index <- which(rowSums(incM_filtered) == 0)
      incM_filtered <- incM_filtered[-filterPeptides_index, ]

      ## Clean memory
      rm(keepProt, filterPeptides_index, exprRNA, exprProts)
      gc()

    }else{
      if (remove == "sharedNoRemove") {

        ## Remove proteins fulfilling the following criteria: 1. their
        ## corresponding transcript is not expressed according to the
        ## sample-matched transcriptome; 2. they are exclusively mapped by
        ## shared and not specific peptides; 3. their shared peptides are also
        ##  mapped on at least one protein with expressed transcript so that
        ## they are not to be removed

        ## PROTEINS WITH NO TRANSCRIPT
        noRNA <- setdiff(colnames(incM), c(exprProts, proteinContam))

        ## PROTEINS WITH NO SPECIFIC PEPTIDE
        ## Extract specific peptides
        specificPep <- which(rowSums(incM) == 1)
        ## Extract proteins with specific peptides
        subIncM <- incM[specificPep, ]
        specificProt <- colnames(subIncM[, which(colSums(subIncM) > 0)])
        onlySharedProt <- setdiff(colnames(incM), specificProt)

        ## PROTEINS WITH NO TRANSCRIPT AND NO SPECIFIC PEPTIDE
        noRNAnoSpecific <- intersect(noRNA, onlySharedProt)

        ## Find peptides potentially lost removing proteins with no transcript
        ## and no specific peptide detected
        incM_noRNAnoSpecific <- incM[, -which(colnames(incM) %in% noRNAnoSpecific)]
        peptidesRemoved <- which(rowSums(incM_noRNAnoSpecific) == 0)
        ## Find proteins which are mapped by potentially lost peptides
        if (length(peptidesRemoved) == 1) {
          prots_peptidesRemoved <- which(incM[peptidesRemoved, ] > 0)
        } else{
          prots_peptidesRemoved <- which(colSums(incM[peptidesRemoved, ]) > 0)
        }
        noRNAnoSpecific_keep <- colnames(incM)[prots_peptidesRemoved]
        ## Remove proteins with no transcript except for those above identified
        noKeep <- setdiff(noRNAnoSpecific, noRNAnoSpecific_keep)
        incM_filtered <- incM[, -which(colnames(incM) %in% noKeep)]
        }
      }
    }

  ## Return output
  return(incM_filtered)

}
