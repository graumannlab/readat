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

goFiles <- downloadGoData(uniProtIds)

goData <- lapply(goFiles, readRDS)

goData <- combineGoData(goData)

# Some values not found using UniProt. Try again using EntrezGene.
notFound <- setdiff(
  uniProtIds,
  lapply(goData, extract2, "UniProtId") %>% unlist(use.names = FALSE) %>% unique
)

entrezGeneIds <- ids %>%
  filter_(~ UniProtId %in% notFound) %$%
  unlist(EntrezGeneId) %>%
  unique()





goFiles2 <- downloadGoData(entrezGeneIds, idType = "EntrezGene")


goData2 <- lapply(goFiles2, readRDS)

goData2 <- combineGoData(goData2)

flatIds <- ids %>%
  unnest_("UniProtId") %>%
  unnest_("EntrezGeneId")




joined <- goData %>%
  lapply(
    function(ns)
    {
      flatIds %>%
        inner_join(
          ns %>% select_(~ -EntrezGeneId),
          by = "UniProtId"
        )
    }
  )

joined2 <- goData %>%
  lapply(
    function(ns)
    {
      flatIds %>%
        inner_join(
          ns %>% select_(~ -UniProtId),
          by = "EntrezGeneId"
        )
    }
  )

go <- Map(
  function(j1, j2)
  {
    bind_rows(j1, j2) %$%
      split(., SeqId) %>%
      lapply(
        function(x)
        {
          x %>%
            select_(~ UniProtId, ~ GoId, ~ GoName, ~ GoDefinition) %>%
            distinct_() %>%
            as.data.frame
        }
      )
  },
  joined,
  joined2
)

save(
  go,
  file = "somalogic/data/go1129.rda",
  compress = "xz"
)






