library(readat)
library(magrittr)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(PFAM.db)

source("inst/scripts/backend.R")

load("data/aptamers.rda")

a <- aptamers %>%
  filter_(~ Type != "Hybridization Control Elution")

entrezGeneIds <- a %$%
  strsplit(EntrezGeneID, " ") %>%
  setNames(a$AptamerId) %>%
  readat:::list_to_data.frame("AptamerId", "EntrezGeneID", stringsAsFactors = FALSE)

pfamData <- entrezGeneIds$EntrezGeneID %>% lapply(
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
) %>% setNames(entrezGeneIds$EntrezGeneID) %>%
  readat:::list_to_data.frame("EntrezGeneID", "PfamId", stringsAsFactors = FALSE)

pfamData$PfamDescription <- pfamData$PfamId %>%
  BiocGenerics::mget(envir = PFAMDE, ifnotfound = NA_character_) %>%
  unlist %>%
  unname

joined <- entrezGeneIds %>%
  inner_join(
    pfamData,
    by = c("EntrezGeneID")
  )

pfam <- joined %>%
  split(.$AptamerId) %>%
  lapply(select_, ~ - AptamerId) %>%
  lapply(distinct_)

# Merge in cases where PFAM ids were not found
notFound <- setdiff(a$AptamerId, names(pfam))
notFoundList <- vector("list", length(notFound)) %>%
  setNames(notFound)
pfam %<>% c(notFoundList)

save(
  pfam,
  file = "data/pfam.rda",
  compress = "xz"
)



