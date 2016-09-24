test_that(
  "is_entrez_gene_symbol matches all HGNC gene symbols",
  {
    symbolFile <- "hgnc_gene_symbols.rds"
    symbols <- readRDS(symbolFile)
    actual <- readat:::is_entrez_gene_symbol(symbols)
    expect_true(all(actual), info = toString(symbols[!actual]))
  }
)
