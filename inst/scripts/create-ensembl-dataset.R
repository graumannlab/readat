library(readat)
library(magrittr)
library(listless)
library(org.Hs.eg.db)
library(dplyr)

source("readat/inst/scripts/backend.R")

load("readat/data/ids1129.rda")

entrezGeneIds <- ids$EntrezGeneId %>%
  setNames(seq_along(.)) %>%
  list_to_data.frame("index", "EntrezGeneId", stringsAsFactors = FALSE)

entrezGeneAndEnsemblIdLookup <- mGetData(entrezGeneIds$EntrezGeneId, org.Hs.egENSEMBL) %>%
  setNames(c("EntrezGeneId", "EnsemblId"))
ensemblIds <- inner_join(
  entrezGeneIds,
  entrezGeneAndEnsemblIdLookup,
  by = "EntrezGeneId"
) %>%
  split(.$index) %>%
  lapply(function(d) unique(d$EnsemblId[!is.na(d$EnsemblId)])) %>%
  setNames(ids$SeqId[as.integer(names(.))])


save(
  ensemblIds,
  file = "readat/data/ensembl.rda",
  compress = "xz"
)





