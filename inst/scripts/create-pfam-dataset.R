library(somalogic)
library(magrittr)
library(listless)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(PFAM.db)

source("somalogic/inst/scripts/backend.R")

load("somalogic/data/ids1129.rda")

entrezGeneIds <- ids$EntrezGeneId %>%
  setNames(ids$SeqId) %>%
  list_to_data.frame("SeqId", "EntrezGeneId", stringsAsFactors = FALSE)

pfam <- entrezGeneIds$EntrezGeneId %>% lapply(
  function(egId)
  {
    tryCatch(
        {
          ids <- AnnotationDbi::select(org.Hs.eg.db, egId, "PFAM")$PFAM
          ids[!is.na(ids)]
        },
        error = function(e) character()
      )
  }
) %>% setNames(entrezGeneIds$EntrezGeneId) %>%
  list_to_data.frame("EntrezGeneId", "PfamId", stringsAsFactors = FALSE)

pfam$PfamDescription <- pfam$PfamId %>%
  BiocGenerics::mget(envir = PFAMDE, ifnotfound = NA_character_) %>%
  unlist %>%
  unname

joined <- entrezGeneIds %>%
  inner_join(
    pfam,
    by = c(EntrezGeneId = "EntrezGeneId1")
  )

pfam <- joined %>%
  split(.$SeqId1) %>%
  lapply(select_, ~ EntrezGeneId, ~ PfamId, ~ PfamDescription)

save(pfam, file = "somalogic/data/pfam1129.rda")



