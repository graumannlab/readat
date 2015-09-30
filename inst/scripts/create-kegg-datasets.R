library(readat)
library(magrittr)
library(listless)
library(dplyr)
library(stringr)
library(biomaRt)
library(assertive)
library(KEGGREST)

source("readat/inst/scripts/backend.R")

load("readat/data/ids1129.rda")

uniProtIds <- ids %>%
  filter_(~ IsHuman) %$%
  unlist(UniProtId) %>%
  unique

keggFiles <- dir(choose.dir(getwd()), full.names = TRUE) # downloadKeggData(uniProtIds)

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

flatIds <- ids %>%
  unnest_("UniProtId") %>%
  unnest_("EntrezGeneId")

joinedDefinitions <- flatIds %>%
  inner_join(
    keggDefinitions,
    by = "UniProtId"
  ) %>%
  select_(~ SeqId, ~ UniProtId, ~ KeggId, ~ KeggDefinition, ~ KeggCytogenicLocation)

keggDefinitions <- joinedDefinitions %>%
  as.data.frame %$%
  split(., SeqId) %>%
  lapply(select_, ~ - SeqId) %>%
  lapply(distinct_)

save(
  keggDefinitions,
  file = "readat/data/keggDefinitions1129.rda",
  compress = "xz"
)


joinedModules <- flatIds %>%
  inner_join(
    keggModules,
    by = "UniProtId"
  ) %>%
  select_(~ SeqId, ~ UniProtId, ~ KeggModuleId, ~ KeggModule) %>%
  lapply(distinct_)

keggModules <- joinedModules %>%
  as.data.frame %$%
  split(., SeqId) %>%
  lapply(select_, ~ - SeqId) %>%
  lapply(distinct_)

save(
  keggModules,
  file = "readat/data/keggModules1129.rda",
  compress = "xz"
)


joinedPathways <- flatIds %>%
  inner_join(
    keggPathways,
    by = "UniProtId"
  ) %>%
  select_(~ SeqId, ~ UniProtId, ~ KeggPathwayId, ~ KeggPathway)

keggPathways <- joinedPathways %>%
  as.data.frame %$%
  split(., SeqId) %>%
  lapply(select_, ~ - SeqId) %>%
  lapply(distinct_)

save(
  keggPathways,
  file = "readat/data/keggPathways1129.rda",
  compress = "xz"
)






