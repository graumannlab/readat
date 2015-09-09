library(somalogic)
library(magrittr)
library(listless)
library(dplyr)
library(stringr)
library(biomaRt)
library(assertive)
library(KEGGREST)

source("somalogic/inst/scripts/backend.R")

load("somalogic/data/ids1129.rda")

uniProtIds <- ids %>%
  filter_(~ IsHuman) %$%
  unlist(UniProtId) %>%
  unique

keggFiles <- downloadKeggData(uniProtIds)

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
  lapply(select_, ~ - SeqId)

save(
  keggDefinitions,
  file = "somalogic/data/keggDefinitions1129.rda",
  compress = "xz"
)


joinedModules <- flatIds %>%
  inner_join(
    keggModules,
    by = "UniProtId"
  ) %>%
  select_(~ SeqId, ~ UniProtId, ~ KeggModuleId, ~ KeggModule)

keggModules <- joinedModules %>%
  as.data.frame %$%
  split(., SeqId) %>%
  lapply(select_, ~ - SeqId)

save(
  keggModules,
  file = "somalogic/data/keggModules1129.rda",
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
  lapply(select_, ~ - SeqId)

save(
  keggPathways,
  file = "somalogic/data/keggPathways1129.rda",
  compress = "xz"
)






