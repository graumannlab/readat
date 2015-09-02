mGetData <- function(x, envir)
{
  x <- as.character(x)
  mget(x[!is.na(x) & nzchar(x)], envir, ifnotfound = NA_character_) %>%
    listless::list_to_data.frame(stringsAsFactors = FALSE)
}


downloadChromosomalData <- function(ids, outdir = tempfile("chromosome"),
  idType = c("UniProt", "EntrezGene"))
{
  idType <- match.arg(idType)
  # Create the 'mart' (ensembl, people)
  ensemblMart <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl" # ensembl code for humans
  )

  # Choose the columns to fetch
  attrs <- c(
    "uniprot_swissprot", "entrezgene",
    "chromosome_name", "start_position", "end_position"
  )

  # Create a place to put them
  dir.create(outdir, recursive = TRUE)

  message("Saving the chromosome files in ", normalizePath(outdir))

  # Since connecting to databases is dangerous, download values one at a time,
  # and save to file
  oneToN <- seq_along(ids)
  outfiles <- file.path(outdir, paste0(oneToN, "_chromosome_", ids, ".rds"))
  tryCatch(
    for(i in oneToN)
    {
      result <- getBM(
        attributes = attrs,
        filters    = switch(
          idType,
          UniProt    = "uniprot_swissprot",
          EntrezGene = "entrezgene"
        ),
        values     = ids[i],
        mart       = ensemblMart,
        uniqueRows = TRUE
      )
      saveRDS(result, outfiles[i])
    },
    error = function(e)
    {
      message(
        sprintf(
          "Failed to retrieve data from ensembl on iteration %d (%s ID = %s).",
          i,
          idType,
          ids[i]
        )
      )
      print(e)
    }
  )

  # Return location of downloaded files
  invisible(outfiles)
}

combineChromosomalData <- function(chromosomalData)
{
  # Chromosome is 1-22 or X
  rxChr <- capture(repeated(char_class(ASCII_DIGIT %R% "X"), 1, 2))

  chromosomalData %>%
    lapply(
      function(x)
      {
        if(nrow(x) == 0)
        {
          return(
            data.frame(
              UniProtId = character(),
              Chromosome = character(),
              StartPosition = integer(),
              EndPosition = integer(),
              stringsAsFactors = FALSE
            )
          )
        }
        x %>%
          mutate_(
            Chromosome   = ~ str_match(chromosome_name, rxChr)[, 2],
            EntrezGeneId = ~ as.character(entrezgene)) %>%
          select_(
            UniProtId     = ~ uniprot_swissprot,
            EntrezGeneId  = ~ EntrezGeneId,
            ~ Chromosome,
            StartPosition = ~ start_position,
            EndPosition   = ~ end_position
          ) %>%
          distinct_()
      }
    ) %>%
    bind_rows() %>%
    as.data.table
}
