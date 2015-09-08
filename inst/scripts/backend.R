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
              EntrezGeneId = character(),
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

downloadGoData <- function(ids, outdir = tempfile("GO"),
  idType = c("UniProt", "EntrezGene"))
{
  idType <- match.arg(idType)
  # Create the 'mart' (ensembl, people)
  ensemblMart <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl" # ensembl code for humans
  )

  # Choose the columns to fetch
  ensemblAttrs <- listAttributes(ensemblMart)
  goAttrs <- with(ensemblAttrs, name[str_detect(description, fixed("GO "))]) #space needed to exclude GOSlim terms
  attrs <- c(
    "uniprot_swissprot", "entrezgene",
    goAttrs
  )

  # Create a place to put them
  dir.create(outdir, recursive = TRUE)

  message("Saving the GO files in ", normalizePath(outdir))

  # Since connecting to databases is dangerous, download values one at a time,
  # and save to file
  oneToN <- seq_along(ids)
  outfiles <- file.path(outdir, paste0(oneToN, "_GO_", ids, ".rds"))
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

combineGoData <- function(goData)
{
  by_namespace <- goData %>%
    lapply(
      function(x)
      {
        if(nrow(x) == 0)
        {
          return(
            data.frame(
              UniProtId        = character(),
              EntrezGeneId     = character(),
              GoId             = character(),
              GoName           = character(),
              GoDefinition     = character(),
              GoNamespace      = character(),
              stringsAsFactors = FALSE
            )
          )
        }
        x %>%
          mutate_(
            EntrezGeneId = ~ as.character(entrezgene)) %>%
          select_(
            UniProtId     = ~ uniprot_swissprot,
            EntrezGeneId  = ~ EntrezGeneId,
            GoId          = ~ go_id,
            GoName        = ~ name_1006,
            GoDefinition  = ~ definition_1006,
            GoNamespace   = ~ namespace_1003
          ) %>%
          distinct_()
      }
    ) %>%
    bind_rows() %>%
    as.data.table %$%
    split(., GoNamespace)
  by_namespace <- by_namespace[nzchar(names(by_namespace))]
  lapply(by_namespace, select_, ~ - GoNamespace)
}

downloadUniprotKeywords <- function(ids, outdir = tempfile("uniprot_keywords"))
{
  mart <- useMart("unimart", "uniprot")

  dir.create(outdir, recursive = TRUE)

  oneToN <- seq_along(ids)
  outfiles <- file.path(outdir, paste0(oneToN, "_uniprot_keywords_", ids, ".rds"))
  for(i in oneToN)
  {
    keywords <- getBM(
      c("accession", "keyword"),
      filters = "accession",
      values  = ids[i],
      mart
    ) %>%
      rename_(UniProtId = ~ accession, Keyword = ~ keyword) %>%
      select_(~ UniProtId, ~ Keyword)
    saveRDS(keywords, outfiles[i])
  }

  saveRDS(keywordData, "uniprotKeywords.rds")

  # Return location of downloaded files
  invisible(outfiles)
}

downloadKeggData <- function(ids, outdir = tempfile("KEGG"))
{
  # Create a place to put them
  dir.create(outdir, recursive = TRUE)

  message("Saving the KEGG files in ", normalizePath(outdir))

  oneToN <- seq_along(uniProtIds)
  outfiles <- file.path(outdir, paste0(oneToN, "_KEGG_", uniProtIds, ".rds"))

  for(i in oneToN)
  {
    message("UniProt ID = ", uniProtIds[i])
    hsaIds <- keggConv("hsa", paste0("up:", uniProtIds[i]))
    if(is_empty(hsaIds))
    {
      message("No KEGG ID corresponding to ", uniProtIds[i])
      next
    }
    result <- keggGet(hsaIds) %>% setNames(hsaIds)
    saveRDS(result, outfiles[i])
  }

  # Return location of downloaded files
  invisible(outfiles)
}

combineKeggDefinitions <- function(keggData, uniProtIds)
{
  Map(
    function(kegg, uniProtId) # loop on a UniProt ID level
    {
      lapply(
        kegg,
        function(keggid)  # loop on a KEGG ID level
        {
          with(
            keggid,
            {
              d <- data.frame(
                UniProtId             = uniProtId,
                KeggId                = unname(ENTRY),
                KeggDefinition        = DEFINITION,
                KeggCytogenicLocation = POSITION,
                stringsAsFactors      = FALSE
              )
            }
          )
        }
      ) %>%
        bind_rows()
    },
    keggData,
    uniProtIds
  ) %>%
    bind_rows() %>%
    distinct_() %>%
    as.data.table
}

combineKeggModules <- function(keggData, uniProtIds)
{
  Map(
    function(kegg, uniProtId) # loop on a UniProt ID level
    {
      lapply(
        kegg,
        function(keggid)  # loop on a KEGG ID level
        {
          if(is.null(keggid$MODULE))
          {
            return(NULL)
          }
          with(
            keggid,
            {
              d <- data.frame(
                UniProtId             = uniProtId,
                KeggModuleId         = names(MODULE),
                KeggModule           = MODULE,
                stringsAsFactors      = FALSE
              )
            }
          )
        }
      ) %>%
        bind_rows()
    },
    keggData,
    uniProtIds
  ) %>%
    bind_rows() %>%
    distinct_() %>%
    as.data.table
}

combineKeggPathways <- function(keggData, uniProtIds)
{
  Map(
    function(kegg, uniProtId) # loop on a UniProt ID level
    {
      lapply(
        kegg,
        function(keggid)  # loop on a KEGG ID level
        {
          if(is.null(keggid$PATHWAY))
          {
            return(NULL)
          }
          with(
            keggid,
            {
              d <- data.frame(
                UniProtId             = uniProtId,
                KeggPathwayId         = names(PATHWAY),
                KeggPathway           = PATHWAY,
                stringsAsFactors      = FALSE
              )
            }
          )
        }
      ) %>%
        bind_rows()
    },
    keggData,
    uniProtIds
  ) %>%
    bind_rows() %>%
    distinct_() %>%
    as.data.table
}
