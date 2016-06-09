#' Get Ensembl IDs by AptamerId
#'
#' Gets the Ensembl IDs associated with SomaLogic aptamer IDs.
#' @param aptamerIds A character vector of SomaLogic aptamer IDs, or \code{NULL}
#' to use all aptamer IDs.
#' @param simplify Logical.  Should the output be collapsed into a single
#' data.frame?
#' @return A list of character vectors.  The names of the list are the input
#' Sequence IDs, and the character vector associated with that element contains the
#' Ensembl IDs.
#' @examples
#' # Each AptamerId may have one, many, or zero associated Ensembl IDs
#' getEnsemblIds(c("2278-61", "4703-87", "4916-2"))
#'
#' # Get everything in the 1129 panel.
#' \dontrun{
#' getEnsemblIds()
#' }
#' @export
getEnsemblIds <- function(aptamerIds = NULL, simplify = FALSE)
{
  y <- getData(aptamerIds, FALSE, "ensembl", "ensemblIds")
  if(simplify)
  {
    list_to_data.frame(y, "SeqId", "EnsemblId")
  } else
  {
    y
  }
}

#' Get UniProt Keywords by AptamerId
#'
#' Gets the UniProt Keywords associated with SomaLogic aptamer IDs.
#' @param aptamerIds A character vector of SomaLogic aptamer IDs, or \code{NULL}
#' to use all aptamer IDs.
#' @param simplify Logical.  Should the output be collapsed into a single
#' data.frame?
#' @return A list of data frames.  The names of the list are the input
#' SeqIds, and the data frame associated with that element contains:
#' \describe{
#' \item{UniProtId}{Character. The UniProt ID that the Keyword is associated
#' with.}
#' \item{Keyword}{Character. A UniProt Keyword associated with the AptamerId and
#' UniProt ID.}
#' }
#' @examples
#' # Each AptamerId may have one, many, or zero associated Ensembl IDs
#' getUniProtKeywords(c("2278-61", "4703-87", "4916-2"))
#'
#' # Get everything in the 1310 and 1129 panels.
#' \dontrun{
#' getUniProtKeywords()
#' }
#' @importFrom dplyr bind_rows
#' @export
getUniProtKeywords <- function(aptamerIds = NULL, simplify = FALSE)
{
  getData(aptamerIds, simplify, "uniprotKeywords", "uniprotKeywords")
}

#' Get Chromosomal Positions by AptamerId
#'
#' Gets the chromosomal positions associated with SomaLogic aptamer IDs.
#' @param aptamerIds A character vector of SomaLogic aptamer IDs, or \code{NULL}
#' to use all aptamer IDs.
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
#' # Each AptamerId may have one, many, or zero associated chromosomal positions
#' getChromosomalPositions(c("2278-61", "4703-87_2", "4916-2"))
#'
#' # Get everything in the 1310 and 1129 panels.
#' \dontrun{
#' getChromosomalPositions()
#' }
#' @importFrom dplyr bind_rows
#' @export
getChromosomalPositions <- function(aptamerIds = NULL, simplify = FALSE)
{
  getData(aptamerIds, simplify, "chromosomalPositions", "chromosomalPositions")
}

#' Get PFAM IDs by AptamerId
#'
#' Gets the PFAM Ids and descriptions associated with SomaLogic aptamer IDs.
#' @param aptamerIds A character vector of SomaLogic aptamer IDs, or \code{NULL}
#' to use all aptamer IDs.
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
#' # Each AptamerId may have one, many, or zero associated PFAM descriptions
#' getPfam(c("2278-61", "4703-87", "4916-2"))
#'
#' # Get everything in the 1310 and 1129 panels.
#' \dontrun{
#' getPfam()
#' }
#' @importFrom dplyr bind_rows
#' @export
getPfam <- function(aptamerIds = NULL, simplify = FALSE)
{
  getData(aptamerIds, simplify, "pfam", "pfam")
}

#' KEGG definitions, modules, and pathways by AptamerId
#'
#' Gets the KEGG definitions, modules, and pathways associated with SomaLogic
#' aptamer IDs.
#' @param aptamerIds A character vector of SomaLogic aptamer IDs, or \code{NULL}
#' to use all aptamer IDs.
#' @param simplify Logical.  Should the output be collapsed into a single
#' data.frame?
#' @return A list of data frames.  The names of the list are the input
#' aptamerIds, and the data frame associated with that element contains:
#'
#' \describe{
#' \item{UniProtId}{Character.  UniProt ID for the protein target.}
#' \item{KeggId}{Character.  KEGG ID for the gene that produces the
#' target protein.}
#' \item{KeggDefinition}{Character.  Description corresponding to the KEGG ID.}
#' \item{KeggCytogenicLocation}{Character.  KEGG determination of the gene's
#' locus.}
#' }
#' @examples
#' # Each AptamerId may have one, many, or zero associated KEGG descriptions
#' getKeggDefinitions(c("2278-61", "3505-6", "4916-2"))
#' getKeggModules(c("2278-61", "3505-6", "4916-2"))
#' getKeggPathways(c("2278-61", "3505-6", "4916-2"))
#'
#' # Get everything in the 1310 and 1129 panels.
#' \dontrun{
#' getKeggDefinitions()
#' getKeggModules()
#' getKeggPathways()
#' }
#' @importFrom dplyr bind_rows
#' @export
getKeggDefinitions <- function(aptamerIds = NULL, simplify = FALSE)
{
  getData(aptamerIds, simplify, "keggDefinitions", "keggDefinitions")
}

#' @rdname getKeggDefinitions
#' @export
getKeggModules <- function(aptamerIds = NULL, simplify = FALSE)
{
  getData(aptamerIds, simplify, "keggModules", "keggModules")
}

#' @rdname getKeggDefinitions
#' @export
getKeggPathways <- function(aptamerIds = NULL, simplify = FALSE)
{
  getData(aptamerIds, simplify, "keggPathways", "keggPathways")
}


#' GO definitions by AptamerId
#'
#' Gets the GO definitions associated with SomaLogic Sequence IDs.  There are
#' three datasets, one for each of these domains: molecular function, biological
#' process, and cellular compartment.
#' @param aptamerIds A character vector of SomaLogic Sequence IDs, or \code{NULL} to
#' use all 1129 Sequence IDs.
#' @param simplify Logical.  Should the output be collapsed into a single
#' data.frame?
#' @return A list of data frames.  The names of the list are the input
#' aptamerIds, and the data frame associated with that element contains:
#' \describe{
#' \item{UniProtId}{Character.  UniProt ID for the protein target.}
#' \item{GoId}{Character.  GO ID for property of the target protein.}
#' \item{GoName}{Character.  Name corresponding to the GO ID.}
#' \item{GoDefinition}{Character.  Description corresponding to the GO ID.}
#' }
#' @examples
#' # Each AptamerId may have one, many, or zero associated GO descriptions
#' getGoMolecularFunctions(c("2278-61", "3505-6", "4916-2"))
#' getGoBiologicalProcesses(c("2278-61", "3505-6", "4916-2"))
#' getGoCellularComponents(c("2278-61", "3505-6", "4916-2"))
#'
#' # Get everything in the 1310 and 1129 panels.
#' \dontrun{
#' getGoMolecularFunctions()
#' getGoBiologicalProcesses()
#' getGoCellularComponents()
#' }
#' @importFrom dplyr bind_rows
#' @export
getGoMolecularFunctions <- function(aptamerIds = NULL, simplify = FALSE)
{
  getData(aptamerIds, simplify, "goMolecularFunction", "goMolecularFunction")
}

#' @rdname getGoMolecularFunctions
#' @export
getGoBiologicalProcesses <- function(aptamerIds = NULL, simplify = FALSE)
{
  getData(aptamerIds, simplify, "goBiologicalProcess", "goBiologicalProcess")
}

#' @rdname getGoMolecularFunctions
#' @export
getGoCellularComponents <- function(aptamerIds = NULL, simplify = FALSE)
{
  getData(aptamerIds, simplify, "goCellularComponent", "goCellularComponent")
}

#' @importFrom utils data
getData <- function(aptamerIds = NULL, simplify = FALSE, dataName, dataElement)
{
  e <- new.env()
  if(is.null(aptamerIds))
  {
    data(list = "aptamers", package = "readat", envir = e)
    aptamerIds = e$aptamers$AptamerId
  }
  data(list = dataName, package = "readat", envir = e)
  y <- e[[dataElement]][names(e[[dataElement]]) %in% aptamerIds]
  if(simplify)
  {
    if(inherits(y, "GRangesList")) # For, e.g., chromosomalPositions
    {
      # See https://support.bioconductor.org/p/83599/#83602
      unlist(y)
    } else
    {
      bind_rows(y, .id = "AptamerId")
    }
  } else
  {
    y
  }
}
