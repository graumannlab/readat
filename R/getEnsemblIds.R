#' Get Ensembl IDs by SeqID
#'
#' Gets the Ensembl IDs associated with SomaLogic Sequence IDs.
#' @param seqIds A character vector of SomaLogic Sequence IDs, or \code{NULL} to use
#' all 1129 Sequence IDs.
#' @param simplify Logical.  Should the output be collapsed into a single
#' data.frame?
#' @return A list of character vectors.  The names of the list are the input
#' Sequence IDs, and the character vector associated with that element contains the
#' Ensembl IDs.
#' @examples
#' # Each SeqID may have one, many, or zero associated Ensembl IDs
#' getEnsemblIds(c("2278-61_4", "4703-87_2", "4916-2_1"))
#'
#' # Get everything in the 1129 panel.
#' \dontrun{
#' getEnsemblIds()
#' }
#' @importFrom listless list_to_data.frame
#' @export
getEnsemblIds <- function(seqIds = NULL, simplify = FALSE)
{
  e <- new.env()
  if(is.null(seqIds))
  {
    data(list = "ids1129", package = "somalogic", envir = e)
    seqIds = e$ids$SeqId
  }
  data(list = "ensembl1129", package = "somalogic", envir = e)
  ensemblIds <- e$ensemblIds[names(e$ensemblIds) %in% seqIds]
  if(simplify)
  {
    list_to_data.frame(ensemblIds, "SeqId", "EnsemblId")
  } else
  {
    ensemblIds
  }
}

#' Get UniProt Keywords by SeqID
#'
#' Gets the UniProt Keywords associated with SomaLogic Sequence IDs.
#' @param seqIds A character vector of SomaLogic Sequence IDs, or \code{NULL} to
#' use all 1129 Sequence IDs.
#' @param simplify Logical.  Should the output be collapsed into a single
#' data.frame?
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
#' getUniProtKeywords(c("2278-61_4", "4703-87_2", "4916-2_1"))
#'
#' # Get everything in the 1129 panel.
#' \dontrun{
#' getUniProtKeywords()
#' }
#' @importFrom dplyr bind_rows
#' @export
getUniProtKeywords <- function(seqIds = NULL, simplify = FALSE)
{
  e <- new.env()
  if(is.null(seqIds))
  {
    data(list = "ids1129", package = "somalogic", envir = e)
    seqIds = e$ids$SeqId
  }
  data(list = "uniprotKeywords1129", package = "somalogic", envir = e)
  uniprotKeywords <- e$uniprotKeywords[names(e$uniprotKeywords) %in% seqIds]
  if(simplify)
  {
    bind_rows(uniprotKeywords, .id = "SeqId")
  } else
  {
    uniprotKeywords
  }
}

#' Get Chromosomal Positions by SeqID
#'
#' Gets the chromosomal positions associated with SomaLogic Sequence IDs.
#' @param seqIds A character vector of SomaLogic Sequence IDs, or \code{NULL} to
#' use all 1129 Sequence IDs.
#' @param simplify Logical.  Should the output be collapsed into a single
#' data.frame?
#' @return A list of data frames.  The names of the list are the input
#' SeqIds, and the data frame associated with that element contains:
#' \describe{
#' \item{UniProtId}{Character. The UniProt ID that the Keyword is associated
#' with.}
#' \item{Keyword}{Character. A UniProt Keyword associated with the SeqID and
#' UniProt ID.}
#' }
#' @examples
#' # Each SeqID may have one, many, or zero associated chromosomal positions
#' getChromosomalPositions(c("2278-61_4", "4703-87_2", "4916-2_1"))
#'
#' # Get everything in the 1129 panel.
#' \dontrun{
#' getChromosomalPositions()
#' }
#' @importFrom dplyr bind_rows
#' @export
getChromosomalPositions <- function(seqIds = NULL, simplify = FALSE)
{
  e <- new.env()
  if(is.null(seqIds))
  {
    data(list = "ids1129", package = "somalogic", envir = e)
    seqIds = e$ids$SeqId
  }
  data(list = "chromosome1129", package = "somalogic", envir = e)
  chromosomalPositions <- e$chromosomalPositions[names(e$chromosomalPositions) %in% seqIds]
  if(simplify)
  {
    bind_rows(chromosomalPositions, .id = "SeqId")
  } else
  {
    chromosomalPositions
  }
}


