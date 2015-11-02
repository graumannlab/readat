library(readat)
library(magrittr)
library(listless)
library(dplyr)
library(data.table)
library(stringr)
library(biomaRt)
library(rebus)
library(tidyr)


source("inst/scripts/backend.R")

load("data/aptamers.rda")

uniProtIds <- aptamers %$%
  unlist(UniProtId) %>%
  unique

goFiles <- downloadGoData(uniProtIds) # dir(choose.dir(getwd()), full.names = TRUE)

goData <- lapply(goFiles, readRDS)

goData <- combineGoData(goData)

# Some values not found using UniProt. Try again using EntrezGene.
notFound <- setdiff(
  uniProtIds,
  lapply(goData, extract2, "UniProtId") %>%
    unlist(use.names = FALSE) %>%
    unique
)

entrezGeneIds <- aptamers %>%
  filter_(~ UniProtId %in% notFound) %$%
  unlist(EntrezGeneId) %>%
  unique()





goFiles2 <- downloadGoData(entrezGeneIds, idType = "EntrezGene") # goFiles2 <- dir(choose.dir(getwd()), full.names = TRUE)


goData2 <- lapply(goFiles2, readRDS)

goData2 <- combineGoData(goData2)

flatIds <- aptamers %>%
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

goMolecularFunction <- go$molecular_function
goBiologicalProcess <- go$biological_process
goCellularComponent <- go$cellular_component

save(
  goMolecularFunction,
  file = "data/goMolecularFunction.rda",
  compress = "xz"
)

save(
  goBiologicalProcess,
  file = "data/goBiologicalProcess.rda",
  compress = "xz"
)

save(
  goCellularComponent,
  file = "data/goCellularComponent.rda",
  compress = "xz"
)






