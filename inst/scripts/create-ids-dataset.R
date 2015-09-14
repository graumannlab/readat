source("somalogic/inst/scripts/backend.R")

file <- file.choose()

sl <- readAdat(file, keepOnlyPasses = FALSE)

seqInfo <- getSequenceInfo(sl)


ids <- seqInfo[
  j = list(
    SeqId = as.character(SeqId),
    SomaId = as.character(SomaId),
    UniProtId = strsplit(as.character(UniProt), " ", fixed = TRUE),
    EntrezGeneId = strsplit(as.character(EntrezGeneID), " ", fixed = TRUE),
    EntrezGeneSymbol = strsplit(as.character(EntrezGeneSymbol), " ", fixed = TRUE),
    Target = as.character(Target),
    TargetFullName = as.character(TargetFullName),
    IsHuman = Organism == "Human"
  )
]

save(
  ids,
  file = "readat/data/ids1129.rda",
  compress = "xz"
)
