library(data.table)
library(dplyr)
# library(stringi)

aptamers <- fread(file.choose(), skip = 2)

# Fix repeated colnames
colnames(aptamers) <- colnames(aptamers) %>%
  make.names(unique = TRUE)

# Only get useful bits
aptamers %<>%
  select_(
    ~ SeqId,
    AptamerId = ~ aptamer, # SeqId, truncated at underscore
    ~ SomaId,
    ~ Target,
    ~ TargetFullName,
    UniProtId = ~ UniProt,
    IsIn1129Panel = ~ X1.1k,
    IsIn1310Panel = ~ X1.3k,
    PlasmaDilution = ~ Plasma,
    SerumDilution = ~ Serum
  ) %>%
  mutate_(
    UniProtId = ~ strsplit(UniProtId, "[, ]+"),
    IsIn1129Panel = ~ as.logical(IsIn1129Panel),
    IsIn1310Panel = ~ as.logical(IsIn1310Panel)
  )
setkeyv(aptamers, "SeqId")

save(aptamers, file = "data/aptamers.rda")
