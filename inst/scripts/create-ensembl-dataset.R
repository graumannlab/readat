library(readat)
library(magrittr)
library(listless)
library(org.Hs.eg.db)
library(dplyr)

source("inst/scripts/backend.R")

load("data/aptamers.rda")

entrezGeneIds <- aptamers$EntrezGeneId %>%
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
  file = "data/ensembl.rda",
  compress = "xz"
)





