#' The Sequence IDs of the Somalogic SOMAscan assay
#'
#' 1129 Sequence IDs of SOMAmers in the Somalogic SOMAscan assay.
#'
#' @docType data
#' @name ids
#' @aliases ids1129
#' @format A data frame with the following columns.
#' \describe{
#' \item{SeqId}{Character, primary key.  The identifier of the SOMAmer sequence.}
#' \item{SomaId}{Character.  The SomaLogic identifier for the protein target.  For the 1129 assay, there is a one-to-one correspondence between SeqId and SomaId, but the for larger assays, there is a many-to-one correspondence.}
#' \item{UniProtId}{List of character vectors. UniProt IDs for the protein target.}
#' \item{EntrezGeneId}{List of character vectors. EntrezGene IDs for the gene that produces the protein target.}
#' \item{IsHuman}{Logical.  \code{TRUE} when the target is a human protein, and \code{FALSE} when it is a virus protein.}
#' }
#' @references The SOMAmers in the SomaLogic SOMAscan assay can be found
#' in this PDF:
#' \url{http://www.somalogic.com/somalogic/media/Assets/PDFs/SSM-011-Rev-10-SOMAscan-Assay-\%28V3-2\%29-Content.pdf}
#' @examples
#' head(ids)
NULL

#' Ensembl IDs by SomaLogic SeqID
#'
#' A lookup of Ensembl IDs by SomaLogic Sequence ID.
#'
#' @docType data
#' @name ensemblIds
#' @aliases ensembl1129
#' @format A list of character vectors.  The names of the list are SomaLogic Seq IDs, and the character vectors contain Ensembl IDs for each Seq ID.
#' @references More information on Ensembl IDs can be found at:
#' \url{http://www.ensembl.org/index.html}
#' @examples
#' head(ensemblIds)
NULL

#' UniProt Keywords by SomaLogic SeqID
#'
#' A lookup of UniProt keywords by SomaLogic Sequence ID.
#'
#' @docType data
#' @name uniprotKeywords
#' @aliases uniprotKeywords1129
#' @format A list of data frames, each with the following columns.
#' \describe{
#' \item{UniProtId}{Character.  UniProt ID for the protein target.}
#' \item{Keyword}{Character.  A UniProt keyword.}
#' }
#' @references More information on UniProt keywords can be found at:
#' \url{http://www.uniprot.org/help/keywords}
#' @examples
#' head(uniprotKeywords)
NULL

#' Chromosomal Positions by SomaLogic SeqID
#'
#' A lookup of chromosomal positions by SomaLogic Sequence ID.
#'
#' @docType data
#' @name chromosomalPositions
#' @aliases chromosome1129
#' @format A list of data frames, each with the following columns.
#' \describe{
#'  \item{UniProtId}{Character.  UniProt ID for the protein target.}
#' \item{Chromsome}{Character.  Either '1' to '22' or 'X' . (Currently no 'Y' values.)}
#' \item{StartPosition}{Integer. Distance in base pairs from the 5' end of the gene to the start of the protein.}
#' \item{EndPosition}{Integer. Distance in base pairs from the 5' end of the gene to the end of the protein.}
#' }
#' @examples
#' head(chromosomalPositions)
NULL

#' PFAM IDs by SomaLogic SeqID
#'
#' A lookup of PFAM (Protein FAMilies) IDs and descriptions by SomaLogic
#' Sequence ID.
#'
#' @docType data
#' @name pfam
#' @aliases pfam1129
#' @format A list of data frames, each with the following columns.
#' \describe{
#' \item{EntrezGeneId}{Character.  EntrezGene IDs for the gene that produces
#' the target protein.}
#' \item{PfamId}{Character.  PFAM ID for a property of the target protein.}
#' \item{PfamDescription}{Character.  Description of a PFAM protein property.}
#' }
#' @references More information on PFAM IDs can be found at:
#' \url{http://pfam.xfam.org/}
#' @examples
#' head(pfam)
NULL

#' KEGG definitions, modules, and pathways by SomaLogic SeqID
#'
#' A lookup of KEGG (Kyoto Encyclopedia of Genes and Genomes) definitions,
#' modules, and pathways by SomaLogic Sequence ID.
#'
#' @docType data
#' @name keggDefinitions
#' @aliases keggModules keggPathways keggDefinitions1129 keggModules1129 keggPathways1129
#' @format A list of data frames, each with the following columns.
#' \describe{
#' \item{UniProtId}{Character.  UniProt ID for the protein target.}
#' \item{KeggId}{Character.  KEGG ID for the gene that produces the
#' target protein.}
#' \item{KeggDefinition}{Character.  Description corresponding to the KEGG ID.}
#' \item{KeggCytogenicLocation}{Character.  KEGG determination of the gene's
#' locus.}
#' }
#' @references More information on KEGG can be found at:
#' \url{http://www.kegg.jp/}
#' @examples
#' head(keggDefinitions)
#' head(keggModules)
#' head(keggPathways)
NULL

#' GO definitions by SomaLogic SeqID
#'
#' A lookup of GO (Gene Ontology) definitions, for molecular function,
#' biological process, and cellular component by SomaLogic Sequence ID.
#'
#' @docType data
#' @name goMolecularFunction
#' @aliases goBiologicalProcess goCellularComponent goMolecularFunction1129 goBiologicalProcess1129 goCellularComponent1129
#' @format A list of data frames, each with the following columns.
#' \describe{
#' \item{UniProtId}{Character.  UniProt ID for the protein target.}
#' \item{GoId}{Character.  GO ID for property of the target protein.}
#' \item{GoName}{Character.  Name corresponding to the GO ID.}
#' \item{GoDefinition}{Character.  Description corresponding to the GO ID.}
#' }
#' @references More information on GO can be found at:
#' \url{http://geneontology.org/}
#' @examples
#' head(goMolecularFunction)
#' head(goBiologicalProcess)
#' head(goCellularComponent)
NULL

