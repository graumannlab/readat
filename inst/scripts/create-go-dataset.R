library(somalogic)
library(magrittr)
library(listless)
library(dplyr)

library(biomaRt)

source("somalogic/inst/scripts/backend.R")

load("somalogic/data/ids1129.rda")


uniProtIds <- seqInfo %>%
  filter_(~ EntrezGeneSymbol != "Human-virus") %$%
  strsplit(as.character(UniProt), " ") %>%
  unlist %>%
  unique

ensemblMart <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl" # ensembl code for humans
)
