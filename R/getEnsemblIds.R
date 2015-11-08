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
  y <- getData(seqIds, FALSE, "ensembl", "ensemblIds")
  if(simplify)
  {
    list_to_data.frame(y, "SeqId", "EnsemblId")
  } else
  {
    y
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
#' # Get everything in the 1310 and 1129 panels.
#' \dontrun{
#' getUniProtKeywords()
#' }
#' @importFrom dplyr bind_rows
#' @export
getUniProtKeywords <- function(seqIds = NULL, simplify = FALSE)
{
  getData(seqIds, simplify, "uniprotKeywords", "uniprotKeywords")
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
#'  \item{UniProtId}{Character.  UniProt ID for the protein target.}
#' \item{Chromsome}{Character.  Either '1' to '22' or 'X' . (Currently no 'Y' values.)}
#' \item{StartPosition}{Integer. Distance in base pairs from the 5' end of the gene to the start of the protein.}
#' \item{EndPosition}{Integer. Distance in base pairs from the 5' end of the gene to the end of the protein.}
#' }
#' @examples
#' # Each SeqID may have one, many, or zero associated chromosomal positions
#' getChromosomalPositions(c("2278-61_4", "4703-87_2", "4916-2_1"))
#'
#' # Get everything in the 1310 and 1129 panels.
#' \dontrun{
#' getChromosomalPositions()
#' }
#' @importFrom dplyr bind_rows
#' @export
getChromosomalPositions <- function(seqIds = NULL, simplify = FALSE)
{
  getData(seqIds, simplify, "chromosome", "chromosomalPositions")
}

#' Get PFAM IDs by SeqID
#'
#' Gets the PFAM Ids and descriptions associated with SomaLogic Sequence IDs.
#' @param seqIds A character vector of SomaLogic Sequence IDs, or \code{NULL} to
#' use all 1129 Sequence IDs.
#' @param simplify Logical.  Should the output be collapsed into a single
#' data.frame?
#' @return A list of data frames.  The names of the list are the input
#' SeqIds, and the data frame associated with that element contains:
#' \describe{
#' \item{EntrezGeneId}{Character.  EntrezGene IDs for the gene that produces
#' the target protein.}
#' \item{PfamId}{Character.  PFAM ID for a property of the target protein.}
#' \item{PfamDescription}{Character.  Description of a PFAM protein property.}
#' }
#' @examples
#' # Each SeqID may have one, many, or zero associated PFAM descriptions
#' getPfam(c("2278-61_4", "4703-87_2", "4916-2_1"))
#'
#' # Get everything in the 1310 and 1129 panels.
#' \dontrun{
#' getPfam()
#' }
#' @importFrom dplyr bind_rows
#' @export
getPfam <- function(seqIds = NULL, simplify = FALSE)
{
  getData(seqIds, simplify, "pfam", "pfam")
}

#' KEGG definitions, modules, and pathways by SeqID
#'
#' Gets the KEGG definitions, modules, and pathways associated with SomaLogic
#' Sequence IDs.
#' @param seqIds A character vector of SomaLogic Sequence IDs, or \code{NULL} to
#' use all Sequence IDs.
#' @param simplify Logical.  Should the output be collapsed into a single
#' data.frame?
#' @return A list of data frames.  The names of the list are the input
#' SeqIds, and the data frame associated with that element contains:
#' \describe{
#' \item{UniProtId}{Character.  UniProt ID for the protein target.}
#' \item{KeggId}{Character.  KEGG ID for the gene that produces the
#' target protein.}
#' \item{KeggDefinition}{Character.  Description corresponding to the KEGG ID.}
#' \item{KeggCytogenicLocation}{Character.  KEGG determination of the gene's
#' locus.}
#' }
#' @examples
#' # Each SeqID may have one, many, or zero associated KEGG descriptions
#' getKeggDefinitions(c("2278-61_4", "3505-6_2", "4916-2_1"))
#' getKeggModules(c("2278-61_4", "3505-6_2", "4916-2_1"))
#' getKeggPathways(c("2278-61_4", "3505-6_2", "4916-2_1"))
#'
#' # Get everything in the 1310 and 1129 panels.
#' \dontrun{
#' getKeggDefinitions()
#' getKeggModules()
#' getKeggPathways()
#' }
#' @importFrom dplyr bind_rows
#' @export
getKeggDefinitions <- function(seqIds = NULL, simplify = FALSE)
{
  getData(seqIds, simplify, "keggDefinitions", "keggDefinitions")
}

#' @rdname getKeggDefinitions
#' @export
getKeggModules <- function(seqIds = NULL, simplify = FALSE)
{
  getData(seqIds, simplify, "keggModules", "keggModules")
}

#' @rdname getKeggDefinitions
#' @export
getKeggPathways <- function(seqIds = NULL, simplify = FALSE)
{
  getData(seqIds, simplify, "keggPathways", "keggPathways")
}


#' GO definitions by SeqID
#'
#' Gets the GO definitions associated with SomaLogic Sequence IDs.  There are
#' three datasets, one for each of these domains: molecular function, biological
#' process, and cellular compartment.
#' @param seqIds A character vector of SomaLogic Sequence IDs, or \code{NULL} to
#' use all 1129 Sequence IDs.
#' @param simplify Logical.  Should the output be collapsed into a single
#' data.frame?
#' @return A list of data frames.  The names of the list are the input
#' SeqIds, and the data frame associated with that element contains:
#' \describe{
#' \item{UniProtId}{Character.  UniProt ID for the protein target.}
#' \item{GoId}{Character.  GO ID for property of the target protein.}
#' \item{GoName}{Character.  Name corresponding to the GO ID.}
#' \item{GoDefinition}{Character.  Description corresponding to the GO ID.}
#' }
#' @examples
#' # Each SeqID may have one, many, or zero associated GO descriptions
#' getGoMolecularFunctions(c("2278-61_4", "3505-6_2", "4916-2_1"))
#' getGoBiologicalProcesses(c("2278-61_4", "3505-6_2", "4916-2_1"))
#' getGoCellularComponents(c("2278-61_4", "3505-6_2", "4916-2_1"))
#'
#' # Get everything in the 1310 and 1129 panels.
#' \dontrun{
#' getGoMolecularFunctions()
#' getGoBiologicalProcesses()
#' getGoCellularComponents()
#' }
#' @importFrom dplyr bind_rows
#' @export
getGoMolecularFunctions <- function(seqIds = NULL, simplify = FALSE)
{
  getData(seqIds, simplify, "goMolecularFunction", "goMolecularFunction")
}

#' @rdname getGoMolecularFunctions
#' @export
getGoBiologicalProcesses <- function(seqIds = NULL, simplify = FALSE)
{
  getData(seqIds, simplify, "goBiologicalProcess", "goBiologicalProcess")
}

#' @rdname getGoMolecularFunctions
#' @export
getGoCellularComponents <- function(seqIds = NULL, simplify = FALSE)
{
  getData(seqIds, simplify, "goCellularComponent", "goCellularComponent")
}

#' @importFrom utils data
getData <- function(seqIds = NULL, simplify = FALSE, dataName, dataElement)
{
  e <- new.env()
  if(is.null(seqIds))
  {
    data(list = "aptamers", package = "readat", envir = e)
    seqIds = e$aptamers$SeqId
  }
  data(list = dataName, package = "readat", envir = e)
  y <- e[[dataElement]][names(e[[dataElement]]) %in% seqIds]
  if(simplify)
  {
    bind_rows(y, .id = "SeqId")
  } else
  {
    y
  }
}
