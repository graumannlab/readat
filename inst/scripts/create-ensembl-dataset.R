library(readat)
library(magrittr)
library(org.Hs.eg.db)
library(dplyr)

source("inst/scripts/backend.R")

load("data/aptamers.rda")

entrezGeneIds <- with(
  aptamers,
  strsplit(EntrezGeneID, " +") %>%
  setNames(AptamerId) %>%
  readat:::list_to_data.frame("AptamerId", "EntrezGeneID", stringsAsFactors = FALSE)
)

entrezGeneAndEnsemblIdLookup <- mGetData(entrezGeneIds$EntrezGeneID, org.Hs.egENSEMBL) %>%
  setNames(c("EntrezGeneID", "EnsemblId"))
ensemblIds <- inner_join(
  entrezGeneIds,
  entrezGeneAndEnsemblIdLookup,
  by = "EntrezGeneID"
) %>%
  split(.$AptamerId) %>%
  lapply(function(d) unique(d$EnsemblId[!is.na(d$EnsemblId)]))


save(
  ensemblIds,
  file = "data/ensembl.rda",
  compress = "xz"
)





