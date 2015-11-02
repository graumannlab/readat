library(readat)
library(magrittr)
library(listless)
library(dplyr)
library(stringr)
library(biomaRt)
library(rebus)
library(tidyr)

source("inst/scripts/backend.R")

load("data/aptamers.rda")

uniProtIds <- aptamers %$%
  unlist(UniProtId) %>%
  unique

chromosomalFiles <- downloadChromosomalData(uniProtIds)

chromosomalData <- lapply(chromosomalFiles, readRDS)

chromosomalData <- combineChromosomalData(chromosomalData)

# Some values not found using UniProt. Try again using EntrezGene.
notFound <- setdiff(uniProtIds, chromosomalData$UniProtId)

entrezGeneIds <- ids %>%
  filter_(~ UniProtId %in% notFound) %$%
  unlist(EntrezGeneId) %>%
  unique()



chromosomalFiles2 <- downloadChromosomalData(entrezGeneIds, idType = "EntrezGene")


chromosomalData2 <- lapply(chromosomalFiles2, readRDS)

chromosomalData2 <- combineChromosomalData(chromosomalData2)

flatIds <- ids %>%
  unnest_("UniProtId") %>%
  unnest_("EntrezGeneId")




joined <- flatIds %>%
  inner_join(
    chromosomalData %>% select_(~ -EntrezGeneId),
    by = "UniProtId"
  )

joined2 <- flatIds %>%
  inner_join(
    chromosomalData2 %>% select_(~ -UniProtId),
    by = "EntrezGeneId"
  )

chromosomalPositions <- bind_rows(joined, joined2) %$%
  split(., SeqId) %>%
  lapply(
    function(x)
    {
      x %>%
        select_(~ UniProtId, ~ Chromosome, ~ StartPosition, ~ EndPosition) %>%
        distinct_() %>%
        as.data.frame
    }
  )

save(
  chromosomalPositions,
  file = "data/chromosome1129.rda",
  compress = "xz"
)
