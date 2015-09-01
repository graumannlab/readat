source("somalogic/inst/scripts/backend.R")

file <- "//acfs/proteomics/Projects/2014/SomaLogic/Data/Original data/WCQ-14-130_20140925/WCQ-14-130_Set_A_RPT.HybMedNormCal_20140925.adat"

sl <- readSomaLogic(file, keepOnlyPasses = FALSE)

seqInfo <- getSequenceInfo(sl)


ids <- seqInfo[
  j = list(
    SeqId,
    SomaId,
    UniProtId = strsplit(as.character(UniProt), " ", fixed = TRUE),
    EntrezGeneId = strsplit(as.character(EntrezGeneID), " ", fixed = TRUE)
  )
]


saveRDS(ids, "somalogic/data/ids1129.rda")
