#' Sequence data for the Somalogic SOMAscan assay
#'
#' Sequence data for the SOMAmers in the 1129 and 1310 panels of the Somalogic
#' SOMAscan assay.
#'
#' @docType data
#' @name aptamers
#' @format A data frame with the following columns.
#' \describe{
#' \item{AptamerId}{Character, primary key.  The identifier of the aptamer sequence
#' identified using SELEX.  It consists of an SomaLogic SeqId truncated at the
#' underscore, to remove its sequence version number.}
#' \item{SomaId}{Character. The SomaLogic identifier for the protein target.
#' For the 1129 and 1310 assays, there is a one-to-one correspondence between
#' SeqId and SomaId, but in theory there is a many-to-one correspondence.}
#' \item{Target}{Character. The name of the protein target, from the protein
#' standard supplier, sometimes with additional annotation by SomaLogic.}
#' \item{TargetFullName}{Character. The name of the protein target, from UniProt.}
#' \item{UniProt}{List of character vectors. UniProt IDs for the protein target.}
#' \item{EntrezGeneID}{List of character vectors. Entrez Gene IDs for the gene
#' associated with protein target.}
#' \item{EntrezGeneSymbol}{List of character vectors. Entrez Gene symbols for
#' the gene associated with protein target.}
#' \item{Organism}{Either "Human" or the name of a virus.}
#' \item{Units}{Should always be RFU, short for Relative Fluorescence Units.}
#' \item{Type}{The protein standard source for the original SELEX performed for
#' the sequence. "Protein" refers to human protein; a few sequences used rat
#' protein standard sources in the SELEX.}
#' \item{PlasmaDilution}{What dilution factor is used in plasma?}
#' \item{SerumDilution}{What dilution factor is used in serum?}
#' \item{IsIn1129Panel}{Is the aptamer in the 1129 panel of the SOMAscan assay?}
#' \item{IsIn1310Panel}{Is the aptamer in the 1310 panel of the SOMAscan assay?}
#' }
#' @references The SOMAmers in the SomaLogic SOMAscan 1310 assay can be found
#' in this PDF:
#' \url{http://www.somalogic.com/somalogic/media/Assets/PDFs/SSM-045-REV-1-SOMAscan-Assay-1-3k-Content.pdf}
#' Those from the 1129 assay can be found here:
#' \url{http://www.somalogic.com/somalogic/media/Assets/PDFs/SSM-011-Rev-11-SOMAscan-Assay-\%28V1-1k\%29-Content.pdf}
#' @examples
#' head(aptamers)
NULL

#' Ensembl IDs by SomaLogic AptamerId
#'
#' A lookup of Ensembl IDs by SomaLogic Aptamer ID.
#'
#' @docType data
#' @name ensemblIds
#' @format A list of character vectors.  The names of the list are SomaLogic Seq IDs, and the character vectors contain Ensembl IDs for each Seq ID.
#' @references More information on Ensembl IDs can be found at:
#' \url{http://www.ensembl.org/index.html}
#' @examples
#' head(ensemblIds)
NULL

#' UniProt Keywords by SomaLogic AptamerId
#'
#' A lookup of UniProt keywords by SomaLogic Aptamer ID.
#'
#' @docType data
#' @name uniprotKeywords
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

#' Chromosomal Positions by SomaLogic AptamerId
#'
#' A lookup of chromosomal positions by SomaLogic Aptamer ID.
#'
#' @docType data
#' @name chromosomalPositions
#' @format A \code{\link[GenomicRanges]{GRangesList}} with the hg38 genome.
#' Each Apatmer ID has an element that is a \code{\link[GenomicRanges]{GRanges}}
#' object with 3 metadata columns.
#' \describe{
#' \item{AptamerId}{Character, primary key.  The identifier of the aptamer sequence
#' identified using SELEX.  It consists of an SomaLogic SeqId truncated at the
#' underscore, to remove its sequence version number.}
#'  \item{UniProtId}{Character.  UniProt ID for the protein target.}
#' \item{EntrezGeneId}{Character.  EntrezGene IDs for the gene that produces
#' the target protein.}
#' }
#' @examples
#' head(chromosomalPositions)
NULL

#' PFAM IDs by SomaLogic AptamerId
#'
#' A lookup of PFAM (Protein FAMilies) IDs and descriptions by SomaLogic
#' Aptamer ID.
#'
#' @docType data
#' @name pfam
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

#' KEGG definitions, modules, and pathways by SomaLogic AptamerId
#'
#' A lookup of KEGG (Kyoto Encyclopedia of Genes and Genomes) definitions,
#' modules, and pathways by SomaLogic Aptamer ID.
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

#' GO definitions by SomaLogic AptamerId
#'
#' A lookup of GO (Gene Ontology) definitions, for molecular function,
#' biological process, and cellular component by SomaLogic Aptamer ID.
#'
#' @docType data
#' @name goMolecularFunction
#' @aliases goBiologicalProcess goCellularComponent
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

