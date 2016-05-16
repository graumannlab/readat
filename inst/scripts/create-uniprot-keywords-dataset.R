library(readat)
library(magrittr)
library(dplyr)
library(stringr)
library(UniProt.ws)
library(rebus)
library(tidyr)
library(data.table)

source("inst/scripts/backend.R")

load("data/aptamers.rda")

uniProtIds <- aptamers %>%
  filter_(~ Type != "Hybridization Control Elution") %$%
  strsplit(UniProt, " ") %>%
  unlist %>%
  unique

homoSapiens <- as.integer(
  availableUniprotSpecies(pattern = "^Homo sapiens$")$`taxon ID`
)
uniprotWebService <- UniProt.ws(taxId = homoSapiens)

keywordData <- select(uniprotWebService, uniProtIds, "KEYWORDS", "UNIPROTKB")

keywordData %<>%
  setNames(c("UniProt", "Keyword")) %>%
  as.data.table %>%
  mutate_(Keyword = ~ strsplit(Keyword, "; ", fixed = TRUE)) %>%
  unnest_("Keyword")

flatIds <- aptamers %>%
  unnest_("UniProt")

joined <- flatIds %>%
  inner_join(
    keywordData,
    by = "UniProt"
  ) %>%
  select_(~ AptamerId, ~ UniProt, ~ Keyword)

uniprotKeywords <- joined %>%
  as.data.frame %$%
  split(., AptamerId) %>%
  lapply(select_, ~ UniProt, ~ Keyword) %>%
  lapply(distinct_)

save(
  uniprotKeywords,
  file = "data/uniprotKeywords.rda",
  compress = "xz"
)
