library(data.table)
library(readat)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(magrittr)
library(stringi)
library(tidyr)

fixId <- function(x)
{
  stri_replace_all_regex(x, "[, ]+", " ")
}

# Start with the dilutions file, then add in HCE sequences + organism etc. info
# From 1.3k & 1.1k adat files.  Need to use pre-processed files from kneadat
# that still contain HCE sequences (almost a circular dependency).

dilutionData <- fread(
  system.file("extdata/SOMAscan_1.3k_1.1k_dilutions.csv", package = "readat")
) %>%
  select_(~ -SeqId, ~ -note) %>%
  mutate_(
    UniProt = ~ fixId(UniProt), # ~ strsplit(UniProt, "[, ]+"),
    IsIn1129Panel = ~ as.logical(IsIn1129Panel),
    IsIn1310Panel = ~ as.logical(IsIn1310Panel)
  ) %>%
  as.data.frame

adatSequenceData <- lapply(
  c("1.3k", "1.1k") %>% setNames(., .),
  function(panel)
  {
    kneadat::readCtrlFile(panel, keepOnlyPasses = FALSE) %>% # extractSampleData(panel)
    getSequenceData %>%
    mutate_(
      AptamerId = ~ convertSeqIdToAptamer(SeqId),
      EntrezGeneID = ~ fixId(EntrezGeneID)
    ) %>%
    select_(~ -SeqId, ~ -Dilution)
  }
) %>%
  bind_rows %>%
  distinct_

aptamers <- dilutionData %>%
  full_join(adatSequenceData, by = c("AptamerId", "SomaId", "Target", "TargetFullName", "UniProt", "Type")) %>%
  mutate_(
    IsHce = ~ Type == "Hybridization Control Elution",
    IsIn1129Panel = ~ ifelse(IsHce, TRUE, IsIn1129Panel),
    IsIn1310Panel = ~ ifelse(IsHce, TRUE, IsIn1310Panel),
    PlasmaDilution = ~ ifelse(IsHce, 0, PlasmaDilution),
    SerumDilution = ~ ifelse(IsHce, 0, SerumDilution),
    Units = ~ "RFU"
  ) %>%
  select_(~ -IsHce)

cols <- c(
  "AptamerId", "SomaId", "Target", "TargetFullName", "UniProt", "EntrezGeneID",
  "EntrezGeneSymbol", "Organism", "Units", "Type", "PlasmaDilution",
  "SerumDilution", "IsIn1310Panel", "IsIn1129Panel"
)
aptamers <- aptamers[, cols]

save(aptamers, file = "data/aptamers.rda")
