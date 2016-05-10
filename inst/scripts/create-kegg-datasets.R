library(assertive)
library(readat)
library(magrittr)
library(dplyr)
library(data.table)
library(stringr)
library(biomaRt)
library(KEGGREST)
library(tidyr)

source("inst/scripts/backend.R")

load("data/aptamers.rda")

uniProtIds <- aptamers %>%
  filter_(~ Type != "Hybridization Control Elution") %$%
  strsplit(UniProt, " ") %>%
  unlist %>%
  unique

# KEGG not working?  https://support.bioconductor.org/p/77847
options(KEGGREST_DEBUG=TRUE)
keggData <- downloadKeggData(uniProtIds)
saveRDS(keggData, "keggData.rds")



keggDefinitions <- combineKeggDefinitions(keggData, uniProtIds)

keggModules <- combineKeggModules(keggData, uniProtIds)

keggPathways <- combineKeggPathways(keggData, uniProtIds)

flatIds <- aptamers %>%
  mutate_(UniProt = ~ strsplit(UniProt, " ")) %>%
  unnest_("UniProt") %>%
  mutate_(EntrezGeneID = ~ strsplit(EntrezGeneID, " ")) %>%
  unnest_("EntrezGeneID")

joinedDefinitions <- flatIds %>%
  inner_join(
    keggDefinitions,
    by = "UniProt",
    copy = TRUE
  ) %>%
  select_(~ AptamerId, ~ UniProt, ~ KeggId, ~ KeggDefinition, ~ KeggCytogenicLocation) %>%
  distinct_

keggDefinitions <- joinedDefinitions %>%
  as.data.frame %$%
  split(., AptamerId) %>%
  lapply(select_, ~ - AptamerId) %>%
  lapply(distinct_)

save(
  keggDefinitions,
  file = "data/keggDefinitions.rda",
  compress = "xz"
)


joinedModules <- flatIds %>%
  inner_join(
    keggModules,
    by = "UniProt"
  ) %>%
  select_(~ AptamerId, ~ UniProt, ~ KeggModuleId, ~ KeggModule) %>%
  distinct_

keggModules <- joinedModules %>%
  as.data.frame %$%
  split(., AptamerId) %>%
  lapply(select_, ~ - AptamerId) %>%
  lapply(distinct_)

save(
  keggModules,
  file = "data/keggModules.rda",
  compress = "xz"
)


joinedPathways <- flatIds %>%
  inner_join(
    keggPathways,
    by = "UniProt"
  ) %>%
  select_(~ AptamerId, ~ UniProt, ~ KeggPathwayId, ~ KeggPathway) %>%
  distinct_

keggPathways <- joinedPathways %>%
  as.data.frame %$%
  split(., AptamerId) %>%
  lapply(select_, ~ - AptamerId) %>%
  lapply(distinct_)

save(
  keggPathways,
  file = "data/keggPathways.rda",
  compress = "xz"
)

