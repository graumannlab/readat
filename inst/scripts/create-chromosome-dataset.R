library(somalogic)
library(magrittr)
library(listless)
library(dplyr)
library(stringr)
library(biomaRt)
library(rebus)
library(tidyr)


source("somalogic/inst/scripts/backend.R")

load("somalogic/data/ids1129.rda")


uniProtIds <- ids %>%
  filter_(~ IsHuman) %$%
  unlist(UniProtId) %>%
  unique

# chromosomalFiles <- downloadChromosomalData(uniProtIds)
chromosomalFiles <- dir("D:\\workspace\\chromosomeec867892e40", full.names = TRUE)

chromosomalData <- lapply(chromosomalFiles, readRDS)

chromosalData <- combineChromosomalData(chromosomalData)

# Some values not found using UniProt. Try again using EntrezGene.
notFound <- setdiff(uniProtIds, chromosalData$UniProtId)

entrezGeneIds <- ids %>%
  filter_(~ UniProtId %in% notFound) %$%
  unlist(EntrezGeneId) %>%
  unique()





chromosomalFiles2 <- downloadChromosomalData(entrezGeneIds, idType = "EntrezGene")


chromosomalData2 <- lapply(chromosomalFiles2, readRDS)

chromosalData2 <- combineChromosomalData(chromosomalData2)

flatIds <- ids %>%
  unnest_("UniProtId") %>%
  unnest_("EntrezGeneId")




joined <- flatIds %>%
  inner_join(
    chromosalData %>% select_(~ -EntrezGeneId),
    by = "UniProtId"
  )

joined2 <- flatIds %>%
  inner_join(
    chromosalData2 %>% select_(~ -UniProtId),
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

save(chromosomalPositions, file = "somalogic/data/chromosome1129.rda")
