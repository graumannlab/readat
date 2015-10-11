source("readat/inst/scripts/backend.R")

file <- file.choose()

sl <- readAdat(file, keepOnlyPasses = FALSE)

seqData <- getSequenceData(sl)


ids <- seqData[
  j = list(
    SeqId = as.character(SeqId),
    SomaId = as.character(SomaId),
    UniProtId = strsplit(as.character(UniProt), " ", fixed = TRUE),
    EntrezGeneId = strsplit(as.character(EntrezGeneID), " ", fixed = TRUE),
    EntrezGeneSymbol = strsplit(as.character(EntrezGeneSymbol), " ", fixed = TRUE),
    Target = as.character(Target),
    TargetFullName = as.character(TargetFullName),
    Type = as.character(Type),
    IsHuman = Organism == "Human"
  )
]

save(
  ids,
  file = "readat/data/ids1129.rda",
  compress = "xz"
)
