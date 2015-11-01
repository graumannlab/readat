library(readat)
library(magrittr)
library(listless)
library(dplyr)
library(stringr)
library(biomaRt)
library(rebus)
library(tidyr)
library(data.table)

source("inst/scripts/backend.R")

load("data/aptamers.rda")

uniProtIds <- aptamers %$%
  unlist(UniProtId) %>%
  unique

keywordFiles <- downloadUniprotKeywords(uniProtIds) # dir(choose.dir(getwd()), full.names = TRUE)

keywordData <- lapply(keywordFiles, readRDS)

keywordData <- keywordData %>%
  bind_rows() %>%
  as.data.table

flatIds <- aptamers %>%
  unnest_("UniProtId")

joined <- flatIds %>%
  inner_join(
    keywordData,
    by = "UniProtId"
  ) %>%
  select_(~ SeqId, ~ UniProtId, ~ Keyword)

uniprotKeywords <- joined %>%
  as.data.frame %$%
  split(., SeqId) %>%
  lapply(select_, ~ UniProtId, ~ Keyword) %>%
  lapply(distinct_)

save(
  uniprotKeywords,
  file = "data/uniprotKeywords.rda",
  compress = "xz"
)
