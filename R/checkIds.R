# See http://www.uniprot.org/help/accession_numbers
UNIPROT_ID <- "^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"

ENTREZ_GENE_ID <- "^(?:[0-9]+)$"

# From http://www.genenames.org/about/guidelines#genesymbols
## 1. The initial character of the symbol should always be a letter.
## 2. Subsequent characters may be other letters, or if necessary, Arabic numerals.
## 3. No punctuation may be used, with the exception of hyphens in specific families and locus types

# From http://www.genenames.org/about/guidelines#specificusage
# @, # punctuation character allowed (all other allowed prefixes/suffixes are
# letters)

# From http://www.genenames.org/about/guidelines#Appendix1
# When referring to splice variants, the symbol can be followed by an
# underscore and the lower case letter "v" then a consecutive number
# This last rule is fudged a little to allow arbitrary letters and numbers after
# the underscore.  There are quite a few cases of hyphen then number in the
# EBI HGNC HS dataset
# ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt

# In the broader scheme of things, Entrez Gene Symbols have no common pattern.
# ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/All_Data.gene_info.gz
ENTREZ_GENE_SYMBOL <- "^(?:[A-Z][A-Z0-9]*)(?:[_-][A-Z0-9]+)*[#@]?$"

# Predicates
#' @importFrom assertive.base call_and_name
#' @importFrom assertive.base set_cause
#' @importFrom stringi stri_detect_regex
is_uniprot_id <- function(x)
{
  x <- as.character(x)
  call_and_name(
    function(x)
    {
      ok <- stri_detect_regex(x, UNIPROT_ID)
      set_cause(ok, "bad format")
    },
    x
  )
}

is_entrez_gene_id <- function(x)
{
  x <- as.character(x)
  call_and_name(
    function(x)
    {
      ok <- stri_detect_regex(x, ENTREZ_GENE_ID)
      set_cause(ok, "bad format")
    },
    x
  )
}

is_entrez_gene_symbol <- function(x)
{
  x <- as.character(x)
  call_and_name(
    function(x)
    {
      ok <- stri_detect_regex(
        x,
        paste0(ENTREZ_GENE_SYMBOL, "|^Human-virus$"),
        case_insensitive = TRUE
      )
      set_cause(ok, "bad format")
    },
    x
  )
}

# assertions
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_all_are_uniprot_ids <- function(x, na_ignore = FALSE,
  severity = getOption("assertive.severity", "stop"))
{
  msg <- gettextf(
    "The values of %s are not all UniProt IDs.",
    get_name_in_parent(x),
    domain = "R-readat"
  )

  assert_engine(
    is_uniprot_id,
    x,
    msg = msg,
    na_ignore = na_ignore,
    severity = severity
  )
}

assert_all_are_entrez_gene_ids <- function(x, na_ignore = FALSE,
  severity = getOption("assertive.severity", "stop"))
{
  msg <- gettextf(
    "The values of %s are not all Entrez Gene IDs.",
    get_name_in_parent(x),
    domain = "R-readat"
  )

  assert_engine(
    is_entrez_gene_id,
    x,
    msg = msg,
    na_ignore = na_ignore,
    severity = severity
  )
}

assert_all_are_entrez_gene_symbols <- function(x, na_ignore = FALSE,
  severity = getOption("assertive.severity", "stop"))
{
  msg <- gettextf(
    "The values of %s are not all Entrez Gene symbols.",
    get_name_in_parent(x),
    domain = "R-readat"
  )
  assert_engine(
    is_entrez_gene_symbol,
    x,
    msg = msg,
    na_ignore = na_ignore,
    severity = severity
  )
}

# High-level checks.

#' Check that the IDs are genuine
#'
#' Checks on character vectors to see if they contain valid UniProt IDs,
#' Entrez Gene IDs, and Entrez Gene symbols.
#' @param sequenceData A data.table of sequence data.  It Should contain columns
#' named "UniProt", "EntrezGeneID", and "EntrezGeneSymbol".
#' @return The character vector of IDs/symbols is silently returned, but the
#' functions are mostly invoked for the side effect of throwing a warning if
#' there are any bad elements.
#' @noRd
#' @importFrom dplyr %>%
checkUniprotIds <- function(sequenceData)
{
  upIds <- strsplit(
    as.character(sequenceData$UniProt),
    " ",
    fixed = TRUE
  ) %>%
    unlist(use.names = FALSE)
  assert_all_are_uniprot_ids(
    ignoreHce(upIds), na_ignore = TRUE, severity = "warning"
  )
}

checkEntrezGeneIds <- function(sequenceData)
{
  egIds <- strsplit(
    as.character(sequenceData$EntrezGeneID),
    " ",
    fixed = TRUE
  ) %>%
    unlist(use.names = FALSE)
  assert_all_are_entrez_gene_ids(
    ignoreHce(egIds), na_ignore = TRUE, severity = "warning"
  )
}


checkEntrezGeneSymbols <- function(sequenceData)
{
  egSymbols <- strsplit(
    as.character(sequenceData$EntrezGeneSymbol),
    " ",
    fixed = TRUE
  ) %>%
    unlist(use.names = FALSE)
  assert_all_are_entrez_gene_symbols(
    ignoreHce(egSymbols), na_ignore = TRUE, severity = "warning"
  )
}

#' @importFrom stringi stri_detect_regex
ignoreHce <- function(x)
{
  x[stri_detect_regex(x, "^HCE[[:digit:]]{6}$")] <- NA_character_
  x
}
