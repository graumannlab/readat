library(readat)
library(magrittr)
library(listless)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(PFAM.db)

source("inst/scripts/backend.R")

load("data/aptamers.rda")


entrezGeneIds <- aptamers$EntrezGeneId %>%
  setNames(aptamers$SeqId) %>%
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
    by = c("EntrezGeneId")
  )

pfam <- joined %>%
  split(.$SeqId) %>%
  lapply(select_, ~ EntrezGeneId, ~ PfamId, ~ PfamDescription) %>%
  lapply(distinct_)

# Merge in cases where PFAM ids were not found
notFound <- setdiff(ids$SeqId, names(pfam))
notFoundList <- vector("list", length(notFound)) %>%
  setNames(notFound)
pfam <- pfam %>% c(notFoundList)

save(
  pfam,
  file = "data/pfam.rda",
  compress = "xz"
)



