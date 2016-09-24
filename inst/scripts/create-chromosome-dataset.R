library(readat)
library(magrittr)
library(dplyr)
library(stringi)
library(biomaRt)
library(rebus)
library(tidyr)
library(GenomicRanges)

source("inst/scripts/backend.R")

load("data/aptamers.rda")

clean_raw <- function(x)
{
  # Single digit chromo has to come after multi-digit chromos in regex
  chromosome_rx <- "(1[0-9]|2[0-2]|[1-9]|X|Y)"
  x %>%
    mutate_(
      Chromosome = ~ stri_extract_first_regex(chromosome_name, paste0("^", chromosome_rx, "$")),
      Chromosome = ~ ifelse(
        is.na(Chromosome),
        stri_match_first_regex(chromosome_name, paste0("CHR", chromosome_rx))[, 2],
        Chromosome
      ),
      Strand = ~ ifelse(strand == 1, "+", "-"),
      entrezgene = ~ as.character(entrezgene)
    ) %>%
    select_(~ - chromosome_name, ~ - strand) %>%
    filter_(
      # 1 record in CHR_HG2030_PATCH doesn't have a chromosome, I think
      ~ !is.na(Chromosome)
    ) %>%
    mutate_(
      # Chromosome numbers prefixed with "chr" to match the values in Seqinfo(genome="hg38")
      Chromosome = ~ paste0("chr", Chromosome)
    ) %>%
    rename_(
      UniProt       = ~ uniprot_swissprot,
      EntrezGeneID  = ~ entrezgene,
      StartPosition = ~ start_position,
      EndPosition   = ~ end_position
    ) %>%
    magrittr::extract(c(1, 2, 5, 3, 4, 6))
}

uniProtIds <- aptamers %>%
  filter_(~ Type != "Hybridization Control Elution") %$%
  strsplit(UniProt, " ") %>%
  unlist %>%
  unique

chromosomalDataRaw <- downloadChromosomalData(uniProtIds)
chromosomalData <- clean_raw(chromosomalDataRaw)



# Some values not found using UniProt. Try again using EntrezGene.
notFound <- setdiff(uniProtIds, chromosomalData$UniProt)

entrezGeneIds <- aptamers %>%
  filter_(~ UniProt %in% notFound, ~ !is.na(EntrezGeneID)) %$%
  unlist(EntrezGeneID) %>%
  unique()



chromosomalDataRaw2 <- downloadChromosomalData(entrezGeneIds, idType = "EntrezGene")(uniProtIds)
chromosomalData2 <- clean_raw(chromosomalDataRaw2)

flatIds <- aptamers %>%
  mutate_(UniProt = ~ strsplit(UniProt, " ")) %>%
  unnest_("UniProt") %>%
  mutate_(EntrezGeneID = ~ strsplit(EntrezGeneID, " ")) %>%
  unnest_("EntrezGeneID")

joined <- flatIds %>%
  inner_join(
    chromosomalData %>% select_(~ -EntrezGeneID),
    by = "UniProt"
  )

joined2 <- flatIds %>%
  inner_join(
    chromosomalData2 %>% select_(~ -UniProt),
    by = "EntrezGeneID"
  )

# chromosomalPositions <- bind_rows(joined, joined2) %$%
#   split(., AptamerId) %>%
#   lapply(
#     function(x)
#     {
#       x %>%
#         select_(~ UniProt, ~ Chromosome, ~ StartPosition, ~ EndPosition) %>%
#         distinct_() %>%
#         as.data.frame
#     }
#   )

# TODO
# 1.  Either rename chromosome values as chr1, chr2, etc. to match rownames of
# Seqinfo(genome = "hg38") or
# Maybe don't process them after download.  Perhaps the extra junk is useful

chromosomalPositions <- bind_rows(joined, joined2) %>%
  select_(~ AptamerId, ~ UniProt, ~ EntrezGeneID, ~ Chromosome, ~ StartPosition, ~ EndPosition, ~ Strand) %>%
  distinct_ %>%
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqinfo = Seqinfo(genome = "hg38"),
    start.field = "StartPosition",
    end.field = "EndPosition",
    strand.field = "Strand"
  )  %$%
  split(., AptamerId)

save(
  chromosomalPositions,
  file = "data/chromosomalPositions.rda",
  compress = "xz"
)
