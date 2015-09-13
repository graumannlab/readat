# See http://www.uniprot.org/help/accession_numbers
UNIPROT_ID <- "(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"

ENTREZ_GENE_ID <- "^(?:[0-9]+)$"

ENTREZ_GENE_SYMBOL <- "^(?:(?:OK/|b|DKFZp)?(?:[A-Z0-9.]+)(?:orf[A-Z0-9]+)?(?:-[A-Z0-9.]+)?)(?:orf[A-Z0-9]+)?(?:_[2HB]|@|-[A-Z0-9]+)?$"

is_uniprot_id <- function(x)
{
  x <- coerce_to(x, "character")
  call_and_name(
    function(x)
    {
      ok <- str_detect(x, UNIPROT_ID)
      set_cause(ok, "bad format")
    },
    x
  )
}

is_entrezgene_id <- function(x)
{
  x <- coerce_to(x, "character")
  call_and_name(
    function(x)
    {
      ok <- str_detect(x, ENTREZ_GENE_ID)
      set_cause(ok, "bad format")
    },
    x
  )
}

is_entrez_gene_symbol <- function(x)
{
  x <- coerce_to(x, "character")
  call_and_name(
    function(x)
    {
      ok <- str_detect(x, ENTREZ_GENE_SYMBOL)
      set_cause(ok, "bad format")
    },
    x
  )
}

assert_all_are_uniprot_ids <- function(x,
  severity = getOption("assertive.severity", "stop"))
{
  msg <- gettextf(
    "The values of %s are not all UniProt IDs.",
    get_name_in_parent(x),
    domain = "R-somalogic"
  )
  assert_engine(is_uniprot_id, x, msg = msg, severity = severity)
}

assert_all_are_entrez_gene_ids <- function(x,
  severity = getOption("assertive.severity", "stop"))
{
  msg <- gettextf(
    "The values of %s are not all Entrez Gene IDs.",
    get_name_in_parent(x),
    domain = "R-somalogic"
  )
  assert_engine(is_entrez_gene_id, x, msg = msg, severity = severity)
}

assert_all_are_entrez_gene_symbols <- function(x,
  severity = getOption("assertive.severity", "stop"))
{
  msg <- gettextf(
    "The values of %s are not all Entrez Gene symbols.",
    get_name_in_parent(x),
    domain = "R-somalogic"
  )
  assert_engine(is_entrez_gene_symbol, x, msg = msg, severity = severity)
}

check_uniprot_ids <- function(sequenceData)
{
  upIds <- strsplit(
    as.character(sequenceData$UniProt),
    " ",
    fixed = TRUE
  ) %>%
    unlist(use.names = FALSE)
  assert_all_are_uniprot_ids(upIds)
}

check_entrez_gene_ids <- function(sequenceData)
{
  egIds <- strsplit(
    as.character(sequenceData$EntrezGeneID),
    " ",
    fixed = TRUE
  ) %>%
    unlist(use.names = FALSE)
  assert_all_are_entrez_gene_ids(egIds)
}


check_entrez_gene_symbols <- function(sequenceData)
{
  egSymbols <- strsplit(
    as.character(sequenceData$EntrezGeneSymbol),
    " ",
    fixed = TRUE
  ) %>%
    unlist(use.names = FALSE)
  assert_all_are_entrez_gene_symbols(egSymbols)
}


