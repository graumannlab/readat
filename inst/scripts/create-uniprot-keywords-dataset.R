library(somalogic)
library(magrittr)
library(listless)
library(dplyr)
library(stringr)
library(biomaRt)
library(rebus)
library(tidyr)
library(data.table)


source("somalogic/inst/scripts/backend.R")

load("somalogic/data/ids1129.rda")


uniProtIds <- ids %>%
  filter_(~ IsHuman) %$%
  unlist(UniProtId) %>%
  unique

keywordFiles <- downloadUniprotKeywords(uniProtIds)

keywordData <- lapply(keywordFiles, readRDS)

keywordData <- keywordData %>%
  bind_rows() %>%
  as.data.table


flatIds <- ids %>%
  unnest_("UniProtId") %>%
  unnest_("EntrezGeneId")

joined <- flatIds %>%
  inner_join(
    keywordData,
    by = "UniProtId"
  ) %>%
  select_(~ SeqId, ~ UniProtId, ~ Keyword)

uniprotKeywords <- joined %>%
  as.data.frame %$%
  split(., SeqId) %>%
  lapply(select_, ~ UniProtId, ~ Keyword)

save(uniprotKeywords, file = "somalogic/data/uniprotKeywords1129.rda")
