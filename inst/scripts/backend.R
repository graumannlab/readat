mGetData <- function(x, envir)
{
  x <- as.character(x)
  mget(x[!is.na(x) & nzchar(x)], envir, ifnotfound = NA_character_) %>%
    readat:::list_to_data.frame(stringsAsFactors = FALSE)
}


downloadChromosomalData <- function(ids, idType = c("UniProt", "EntrezGene"))
{
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

  getBM(
    attributes = attrs,
    filters    = switch(
      idType,
      UniProt    = "uniprot_swissprot",
      EntrezGene = "entrezgene"
    ),
    values     = ids,
    mart       = ensemblMart,
    uniqueRows = TRUE
  )
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
  uniProtIds <- aptamers %>%
    filter_(~ Type != "Hybridization Control Elution") %$%
    strsplit(UniProt, " ") %>%
    unlist %>%
    unique

  idType <- "UniProt"
  ensemblMart <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl" # ensembl code for humans
  )

  ensemblAttrs <- listAttributes(ensemblMart)
  # go attrs, ignoring "go_linkage_type"
  goAttrs <- c("go_id", "name_1006", "definition_1006", "namespace_1003")


  getBM(
    attributes = c(idType, goAttrs),
    filters    = switch(
      idType,
      UniProt    = "uniprot_swissprot",
      EntrezGene = "entrezgene"
    ),
    values     = ids,
    mart       = ensemblMart,
    uniqueRows = TRUE
  )
}


downloadUniprotKeywords <- function(ids, outdir = tempfile("uniprot_keywords"))
{
  mart <- useMart("unimart", "uniprot")

  message("Saving files to ", outdir)
  dir.create(outdir, recursive = TRUE)

  oneToN <- seq_along(ids)
  outfiles <- file.path(outdir, paste0(oneToN, "_uniprot_keywords_", ids, ".rds"))

  for(i in oneToN)
  {
    tryCatch(
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
      },
      error = function(e)
      {
        message(
          sprintf(
            "Failed to retrieve data from UniProt on iteration %d (ID = %s).",
            i,
            ids[i]
          )
        )
        print(e)
      }
    )
  }

  # Return location of downloaded files
  invisible(outfiles)
}

downloadKeggData <- function(uniProtIds)
{
  lapply(
    seq_along(uniProtIds),
    function(i)
    {
      message("UniProt ID = ", uniProtIds[i])
      hsaIds <- keggConv("hsa", paste0("uniprot:", uniProtIds[i]))
      if(is_empty(hsaIds))
      {
        message("No KEGG ID corresponding to ", uniProtIds[i])
        return(character())
      }
      keggGet(hsaIds) %>% setNames(hsaIds)
    }
  )
}

combineKeggDefinitions <- function(keggData, uniProtIds)
{
  Map(
    function(kegg, uniProtId) # loop on a UniProt ID level
    {
      lapply(
        kegg,
        function(keggDataI)  # loop on a KEGG ID level
        {
          with(
            keggDataI,
            {
              d <- data.frame(
                UniProt               = uniProtId,
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
        function(keggDataI)  # loop on a KEGG ID level
        {
          if(is.null(keggDataI$MODULE))
          {
            return(NULL)
          }
          with(
            keggDataI,
            {
              d <- data.frame(
                UniProt             = uniProtId,
                KeggModuleId        = names(MODULE),
                KeggModule          = MODULE,
                stringsAsFactors    = FALSE
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
        function(keggDataI)  # loop on a KEGG ID level
        {
          if(is.null(keggDataI$PATHWAY))
          {
            return(NULL)
          }
          with(
            keggDataI,
            {
              d <- data.frame(
                UniProt               = uniProtId,
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
