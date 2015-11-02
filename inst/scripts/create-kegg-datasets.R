library(assertive)
library(readat)
library(magrittr)
library(listless)
library(dplyr)
library(data.table)
library(stringr)
library(biomaRt)
library(KEGGREST)
library(tidyr)

source("inst/scripts/backend.R")

load("data/aptamers.rda")

uniProtIds <- aptamers %$%
  unlist(UniProtId) %>%
  unique


keggFiles <- downloadKeggData(uniProtIds) # dir(choose.dir(getwd()), full.names = TRUE)
keggFiles <- keggFiles[file.exists(keggFiles)]

keggData <- lapply(keggFiles, readRDS)


matches <- str_match(
    basename(keggFiles),
    "([[:digit:]]+)_KEGG_([[:alnum:]]+)\\.rds"
  )
foundUniProtIds  <- matches %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  setNames(c("FullMatch", "Number", "UniProtId")) %>%
  mutate_(Number = ~ as.numeric(Number)) %>%
  arrange_(~ Number) %>%
  extract2("UniProtId")

keggDefinitions <- combineKeggDefinitions(keggData, foundUniProtIds)

keggModules <- combineKeggModules(keggData, foundUniProtIds)

keggPathways <- combineKeggPathways(keggData, foundUniProtIds)

flatIds <- aptamers %>%
  unnest_("UniProtId") %>%
  unnest_("EntrezGeneId")

joinedDefinitions <- flatIds %>%
  inner_join(
    keggDefinitions,
    by = "UniProtId",
    copy = TRUE
  ) %>%
  select_(~ SeqId, ~ UniProtId, ~ KeggId, ~ KeggDefinition, ~ KeggCytogenicLocation) %>%
  distinct_

keggDefinitions <- joinedDefinitions %>%
  as.data.frame %$%
  split(., SeqId) %>%
  lapply(select_, ~ - SeqId) %>%
  lapply(distinct_)

save(
  keggDefinitions,
  file = "data/keggDefinitions.rda",
  compress = "xz"
)


joinedModules <- flatIds %>%
  inner_join(
    keggModules,
    by = "UniProtId"
  ) %>%
  select_(~ SeqId, ~ UniProtId, ~ KeggModuleId, ~ KeggModule) %>%
  distinct_

keggModules <- joinedModules %>%
  as.data.frame %$%
  split(., SeqId) %>%
  lapply(select_, ~ - SeqId) %>%
  lapply(distinct_)

save(
  keggModules,
  file = "data/keggModules1129.rda",
  compress = "xz"
)


joinedPathways <- flatIds %>%
  inner_join(
    keggPathways,
    by = "UniProtId"
  ) %>%
  select_(~ SeqId, ~ UniProtId, ~ KeggPathwayId, ~ KeggPathway) %>%
  distinct_

keggPathways <- joinedPathways %>%
  as.data.frame %$%
  split(., SeqId) %>%
  lapply(select_, ~ - SeqId) %>%
  lapply(distinct_)

save(
  keggPathways,
  file = "data/keggPathways.rda",
  compress = "xz"
)

