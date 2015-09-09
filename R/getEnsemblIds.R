#' Get Ensembl IDs by SeqID
#'
#' Gets the Ensembl IDs associated with SomaLogic SeIDs.
#' @param seqIds A character vector of SomaLogic SeqIDs, or \code{NULL} to use
#' all 1129 SeqIDs.
#' @return A list of character vectors.  The names of the list are the input
#' SeqIDs, and the character vector associated with that element contains the
#' Ensembl IDs.
#' @examples
#' # Each SeqID may have one, many, or zero associated Ensembl IDs
#' getEnsemblIds(c("2278-61_4", "4703-87_2", "4916-2_1"))
#'
#' # Get everything in the 1129 panel.
#' getEnsemblIds()
#' @export
getEnsemblIds <- function(seqIds = NULL)
{
  e <- new.env()
  if(is.null(seqIds))
  {
    data(list = "ids1129", package = "somalogic", envir = e)
    seqIds = e$ids$SeqId
  }
  data(list = "ensembl1129", package = "somalogic", envir = e)
  e$ensemblIds[names(e$ensemblIds) %in% seqIds]
}

#' Get UniProt Keywords by SeqID
#'
#' Gets the UniProt Keywords associated with SomaLogic SeIDs.
#' @param seqIds A character vector of SomaLogic Sequence IDs, or \code{NULL} to
#' use all 1129 Sequence IDs.
#' @return A list of data frames.  The names of the list are the input
#' SeqIds, and the data frame associated with that element contains:
#' \describe{
#' \item{UniProtId}{Character. The UniProt ID that the Keyword is associated
#' with.}
#' \item{Keyword}{Character. A UniProt Keyword associated with the SeqID and
#' UniProt ID.}
#' }
#' @examples
#' # Each SeqID may have one, many, or zero associated Ensembl IDs
#' getEnsemblIds(c("2278-61_4", "4703-87_2", "4916-2_1"))
#'
#' # Get everything in the 1129 panel.
#' getEnsemblIds()
#' @export
getUniProtKeywords <- function(seqIds = NULL)
{
  e <- new.env()
  if(is.null(seqIds))
  {
    data(list = "ids1129", package = "somalogic", envir = e)
    seqIds = e$ids$SeqId
  }
  data(list = "uniprotKeywords1129", package = "somalogic", envir = e)
  e$uniprotKeywords[names(e$uniprotKeywords) %in% seqIds]
}



