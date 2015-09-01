library(somalogic)
library(magrittr)
library(listless)
library(dplyr)
library(org.Hs.eg.db)

source("somalogic/inst/scripts/backend.R")

load("somalogic/data/ids1129.rda")

entrezGeneIds <- seqInfo %$%
  strsplit(as.character(EntrezGeneID), " ") %>%
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
  #lapply(function(d) paste(d$EnsemblId, collapse = " ")) %>%
  lapply(function(d) unique(d$EnsemblId[!is.na(d$EnsemblId)])) %>%
  setNames(seqInfo$SeqId[as.integer(names(.))])


saveRDS(ensemblIds, "somalogic/data/ensembl1129.rda")





