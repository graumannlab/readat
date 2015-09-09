library(somalogic)
library(magrittr)
library(listless)
library(dplyr)
library(org.Hs.eg.db)

source("somalogic/inst/scripts/backend.R")

load("somalogic/data/ids1129.rda")

entrezGeneIds <- ids$EntrezGeneId %>%
  setNames(seq_along(.)) %>%
  list_to_data.frame("index", "EntrezGeneId")

entrezGeneAndEnsemblIdLookup <- mGetData(entrezGeneIds$EntrezGeneId, org.Hs.egENSEMBL) %>%
  setNames(c("EntrezGeneId", "EnsemblId"))
ensemblIds <- inner_join(
  entrezGeneIds,
  entrezGeneAndEnsemblIdLookup,
  by = "EntrezGeneId"
) %>%
  split(.$index1) %>%
  lapply(function(d) unique(d$EnsemblId[!is.na(d$EnsemblId)])) %>%
  setNames(ids$SeqId[as.integer(names(.))])


save(
  ensemblIds,
  file = "somalogic/data/ensembl1129.rda",
  compress = "xz"
)





