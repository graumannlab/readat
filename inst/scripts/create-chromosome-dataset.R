library(readat)
library(magrittr)
library(dplyr)
library(stringr)
library(biomaRt)
library(rebus)
library(tidyr)

source("inst/scripts/backend.R")

load("data/aptamers.rda")

uniProtIds <- aptamers %>%
  filter_(~ Type != "Hybridization Control Elution") %$%
  strsplit(UniProt, " ") %>%
  unlist %>%
  unique

chromosomalData <- downloadChromosomalData(uniProtIds)

chromosomalData %<>%
  mutate_(chromosome_name = ~ str_extract(chromosome_name, "(\\d{1,2}|X)")) %>%
  rename_(
    UniProt       = ~ uniprot_swissprot,
    EntrezGeneID  = ~ entrezgene,
    Chromosome    = ~ chromosome_name,
    StartPosition = ~ start_position,
    EndPosition   = ~ end_position
  )


# Some values not found using UniProt. Try again using EntrezGene.
notFound <- setdiff(uniProtIds, chromosomalData$UniProtId)

entrezGeneIds <- aptamers %>%
  filter_(~ UniProt %in% notFound, ~ !is.na(EntrezGeneID)) %$%
  unlist(EntrezGeneID) %>%
  unique()



chromosomalData2 <- downloadChromosomalData(entrezGeneIds, idType = "EntrezGene")

chromosomalData2 %<>%
  mutate_(
    chromosome_name = ~ str_extract(chromosome_name, "(\\d{1,2}|X)"),
    entrezgene = ~ as.character(entrezgene)
  ) %>%
  rename_(
    UniProt       = ~ uniprot_swissprot,
    EntrezGeneID  = ~ entrezgene,
    Chromosome    = ~ chromosome_name,
    StartPosition = ~ start_position,
    EndPosition   = ~ end_position
  )

flatIds <- aptamers %>%
  mutate_(UniProt = ~ strsplit(UniProt, " ")) %>%
  unnest_("UniProt") %>%
  mutate_(EntrezGeneID = ~ strsplit(EntrezGeneID, " ")) %>%
  unnest_("EntrezGeneID")




joined <- flatIds %>%
  inner_join(
    chromosomalData %>% select_(~ -EntrezGeneID),
    by = "UniProt"
  )

joined2 <- flatIds %>%
  inner_join(
    chromosomalData2 %>% select_(~ -UniProt),
    by = "EntrezGeneID"
  )

chromosomalPositions <- bind_rows(joined, joined2) %$%
  split(., AptamerId) %>%
  lapply(
    function(x)
    {
      x %>%
        select_(~ UniProt, ~ Chromosome, ~ StartPosition, ~ EndPosition) %>%
        distinct_() %>%
        as.data.frame
    }
  )

save(
  chromosomalPositions,
  file = "data/chromosome1129.rda",
  compress = "xz"
)
