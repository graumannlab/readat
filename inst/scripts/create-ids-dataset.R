source("somalogic/inst/scripts/backend.R")

file <- system.file("extdata", "WEI_15-046_20150330.adat", package = "koraproteomics")

sl <- readSomaLogic(file, keepOnlyPasses = FALSE)

seqInfo <- getSequenceInfo(sl)


ids <- seqInfo[
  j = list(
    SeqId,
    SomaId,
    UniProtId = strsplit(as.character(UniProt), " ", fixed = TRUE),
    EntrezGeneId = strsplit(as.character(EntrezGeneID), " ", fixed = TRUE),
    IsHuman = Organism == "Human"
  )
]

save(ids, file = "somalogic/data/ids1129.rda")
