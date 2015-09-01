library(somalogic)
library(magrittr)
library(listless)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)

source("somalogic/inst/scripts/backend.R")

file <- "//acfs/proteomics/Projects/2014/SomaLogic/Data/Original data/WCQ-14-130_20140925/WCQ-14-130_Set_A_RPT.HybMedNormCal_20140925.adat"

sl <- readSomaLogic(file, keepOnlyPasses = FALSE)

seqInfo <- getSequenceInfo(sl)


# Ensembl IDs

entrezGeneIds <- seqInfo %$%
  strsplit(as.character(EntrezGeneID), " ") %>%
  setNames(seq_along(.)) %>%
  list_to_data.frame("index", "EntrezGeneId")

entrezGeneAndEnsemblIdLookup <- mGetData(entrezGeneIds$EntrezGeneId, org.Hs.egENSEMBL) %>%
  setNames(c("EntrezGeneId", "EnsemblId"))
ensemblIds <- inner_join(
  entrezGeneIds,
  entrezGeneAndEnsemblIdLookup,
  by = "EntrezGeneId"
) %>%
  split(.$index1) %>%
  #lapply(function(d) paste(d$EnsemblId, collapse = " ")) %>%
  lapply(function(d) unique(d$EnsemblId[!is.na(d$EnsemblId)])) %>%
  setNames(seqInfo$SeqId[as.integer(names(.))])


saveRDS(ensemblIds, "somalogic/data/ensembl1129.rda")




# GO terms

uniProtIds <- seqInfo %>%
  filter_(~ EntrezGeneSymbol != "Human-virus") %$%
  strsplit(as.character(UniProt), " ") %>%
  unlist %>%
  unique

ensemblMart <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl" # ensembl code for humans
)

