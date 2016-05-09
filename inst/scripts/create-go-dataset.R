library(readat)
library(magrittr)
library(dplyr)
library(data.table)
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

idType <- "UniProt"
ensemblMart <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl" # ensembl code for humans
)

ensemblAttrs <- listAttributes(ensemblMart)
# go attrs, ignoring "go_linkage_type"
goAttrs <- c("go_id", "name_1006", "definition_1006", "namespace_1003")


goData1 <- getBM(
  attributes = c("uniprot_swissprot", goAttrs),
  filters    = "uniprot_swissprot",
  values     = uniProtIds,
  mart       = ensemblMart,
  uniqueRows = TRUE
)

goData1 %<>%
  rename_(
    UniProt = ~ uniprot_swissprot,
    GoId = ~ go_id,
    GoName = ~ name_1006,
    GoDefinition = ~ definition_1006,
    GoNamespace = ~ namespace_1003
  )




# Some values not found using UniProt. Try again using EntrezGene.
notFound <- setdiff(uniProtIds, unique(goData$UniProt))

entrezGeneIds <- aptamers %>%
  filter_(~ UniProt %in% notFound, ~!is.na(EntrezGeneID), ~ Type != "Hybridization Control Elution") %$%
  strsplit(EntrezGeneID, " ") %>%
  unlist %>%
  unique

goData2 <- getBM(
  attributes = c("entrezgene", goAttrs),
  filters    = "entrezgene",
  values     = entrezGeneIds,
  mart       = ensemblMart,
  uniqueRows = TRUE
)

goData2 %<>%
  rename_(
    EntrezGeneID = ~ entrezgene,
    GoId = ~ go_id,
    GoName = ~ name_1006,
    GoDefinition = ~ definition_1006,
    GoNamespace = ~ namespace_1003
  ) %>%
  mutate_(EntrezGeneID = ~ as.character(EntrezGeneID))



flatIds <- aptamers %>%
  mutate_(UniProt = ~ strsplit(UniProt, " ")) %>%
  unnest_("UniProt") %>%
  mutate_(EntrezGeneID = ~ strsplit(EntrezGeneID, " ")) %>%
  unnest_("EntrezGeneID")




joined <- goData1 %>%
  filter_(~ nzchar(GoNamespace)) %$%
  split(., GoNamespace) %>%
  lapply(
    function(ns)
    {
      flatIds %>%
        inner_join(ns, by = "UniProt")
    }
  )

joined2 <- goData2 %>%
  filter_(~ nzchar(GoNamespace)) %$%
  split(., GoNamespace) %>%
  lapply(
    function(ns)
    {
      flatIds %>%
        inner_join(ns, by = "EntrezGeneID")
    }
  )

go <- Map(
  function(j1, j2)
  {
    bind_rows(j1, j2) %$%
      split(., AptamerId) %>%
      lapply(
        function(x)
        {
          x %>%
            select_(~ UniProt, ~ GoId, ~ GoName, ~ GoDefinition) %>%
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






