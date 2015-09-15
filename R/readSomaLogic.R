utils::globalVariables("SeqId")
utils::globalVariables("SomaId")
utils::globalVariables("Target")
utils::globalVariables("TargetFullName")
utils::globalVariables("UniProt")
utils::globalVariables("EntrezGeneID")
utils::globalVariables("EntrezGeneSymbol")
utils::globalVariables("Organism")
utils::globalVariables("Units")
# utils::globalVariables("CalReference")
# utils::globalVariables("Cal_Set_A_RPT")
utils::globalVariables("ColCheck")
# utils::globalVariables("Dilution")
utils::globalVariables("PlateId")
utils::globalVariables("SlideId")
utils::globalVariables("SampleId")
utils::globalVariables("SampleType")
utils::globalVariables("SampleMatrix")
utils::globalVariables("Barcode")
utils::globalVariables("Barcode2d")
utils::globalVariables("SampleNotes")
utils::globalVariables("SampleDescription")
utils::globalVariables("TimePoint")
utils::globalVariables("ExtIdentifier")
utils::globalVariables("SampleGroup")
utils::globalVariables("SiteId")
utils::globalVariables("SampleUniqueID")
utils::globalVariables("Subject_ID")
utils::globalVariables("RowCheck")
utils::globalVariables("Intensity")

#' Read a SomaLogic data file
#'
#' Reads a SomaLogic ADAT data file.
#' @param file A string containing the path to the file to be read.
#' @param keepOnlyPasses A logical value indicating whether or not to keep
#' only the rows and columns where the data quality was considered to be
#' passable.
#' @param dateFormat A string describing the format of the dates contained in
#' the file's metadata.  See \code{\link[base]{strptime}} for how to specify
#' these.
#' @return An object of class \code{WideSomaLogicData}, which inherits from
#' \code{data.table}.
#' The return value consists of a data frame where each row represents a
#' sample.  Initial columns contain sample metadata and later columns contain
#' intensities of proteins.  The specific metadata columns are not fixed, but
#' the most useful ones described below should always be present.
#' \describe{
#' \item{SampleUniqueID}{A unique identifier for the sample.}
#' \item{Subject_ID}{A unique identifier for the person being sampled.}
#' \item{RowCheck}{This should have the value "PASS" if the data quality is
#' acceptable.}
#' }
#' Columns of proteins intensities have a name beginning \code{SeqId.}.
#' Return value also has three attributes.
#' \describe{
#' \item{SequenceInfo}{A data frame }
#' \item{Metadata}{A list of experimental metadata values.}
#' \item{Checksum}{A SHA1 checksum to ensure file integrity.}
#' }
#' @examples
#' \donttest{
#' unzip(
#'   system.file("extdata", "soma_atkin_diabetes.zip", package = "koraproteomics"),
#'   exdir = tempdir()
#' )
#' soma_file <- file.path(tempdir(), "soma_atkin_diabetes.adat")
#' wide_soma_data <- readSomaLogic(soma_file)
#' str(wide_soma_data, list.len = 35)
#' unlink(soma_file)
#' }
#' @importFrom assertive assert_any_are_true
#' @importFrom assertive is_a_string
#' @importFrom assertive is_readable_connection
#' @importFrom data.table ":="
#' @importFrom data.table as.data.table
#' @importFrom data.table fread
#' @importFrom data.table setattr
#' @importFrom data.table setkey
#' @importFrom data.table setnames
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace_all
#' @export
#' @author Richard Cotton
readAdat <- function(file, keepOnlyPasses = TRUE, dateFormat = "%d/%m/%Y")
{
  assert_any_are_true(c(is_a_string(file), is_readable_connection(file)))

  # The file is split into several groups of data:
  # A checksum, header data, column data and row data.  Read each separately.
  nFields <- count.fields(file, sep = "\t", quote = "")
  dataGroupRow <- which(nFields == 1)

  if(length(dataGroupRow) < 4L)
  {
    stop("The input file is malformed; there should be four rows with only one column.")
  }

  # Read SHA1 checksum
  # For a single line, read.table is faster than data.table::fread. Compare
  # microbenchmark(
  #   read.table = read.table(file, sep = "\t", nrows = 1, stringsAsFactors = FALSE),
  #   fread = fread(file, sep = "\t", nrows = 1, header = FALSE, stringsAsFactors = FALSE, skip = 0)
  # )
  checksum <- read.table(
    file,
    sep              = "\t",
    nrows            = 1,
    stringsAsFactors = FALSE
  )$V2[1]


  # Read header
  # (Ab)using fread for this tends to results in unnecessary warnings.
  headerData <- suppressWarnings(fread(
    file,
    sep        = "\t",
    nrows      = dataGroupRow[2] - dataGroupRow[1] - 1,
    header     = FALSE,
    skip       = dataGroupRow[1],
    integer64  = 'numeric',
    na.strings = c("", "NA")
  ))
  header <- with(headerData, setNames(as.list(V2), substring(V1, 2)))
  header <- within(
    header,
    {
      Version <- as.package_version(Version)
      CreatedDate <- as.Date(CreatedDate, format = dateFormat)
      ExpDate <- as.Date(ExpDate, format = dateFormat)
      ProteinEffectiveDate <- as.Date(ProteinEffectiveDate, format = dateFormat)
    }
  )


  # Read column data
  # This causes an erroneous false-positive warning about not reading to the end
  # of the file.  See https://github.com/Rdatatable/data.table/issues/1330
  sequenceData <- fread(
    file,
    sep              = "\t",
    nrows            = nFields[dataGroupRow[2] + 1] - 1,
    colClasses       = "character",
    skip             = dataGroupRow[4],
    header           = FALSE,
    stringsAsFactors = FALSE,
    integer64        = 'numeric',
    na.strings      = c("", "NA")
  )
  # Get the column that contains the headers
  sequenceHeaderColumnNumber <- nFields[dataGroupRow[3] + 1]
  sequenceHeaderNames <- sequenceData[[sequenceHeaderColumnNumber]]

  # Remove leading blank columns
  sequenceData <- sequenceData[
    j = -seq_len(sequenceHeaderColumnNumber),
    with = FALSE
  ]

  #Transpose, and undo the conversion to matrix
  sequenceData <- as.data.table(t(sequenceData))
  setnames(sequenceData, sequenceHeaderNames)

  # Update column types. Columns change with different SOMA versions.
  # Therefore, only update columns which we currently use, leave others unchanged
  sequenceData <- sequenceData[
    j = `:=`(
      SeqId            = factor(SeqId),
      SomaId           = factor(SomaId),
      Target           = factor(Target),
      TargetFullName   = factor(TargetFullName),
      UniProt          = factor(str_replace_all(UniProt, "[, ]+", " ")),
      EntrezGeneID     = factor(str_replace_all(EntrezGeneID, "[, ]+", " ")),
      EntrezGeneSymbol = factor(str_replace_all(EntrezGeneSymbol, "[, ]+", " ")),
      Organism         = factor(Organism),
      Units            = factor(Units),
      ColCheck         = factor(ColCheck),
      CalReference     = as.numeric(CalReference),
      Dilution         = as.numeric(Dilution)
    )
  ]

  # Check for bad UniProt and EntrezGene Ids.
  check_uniprot_ids(sequenceData)
  check_entrez_gene_ids(sequenceData)
  check_entrez_gene_symbols(sequenceData)

  # There are some more columns that need fixing, which should have the names
  # sprintf("Cal_%s", str_replace(levels(intensityData$PlateId), " ", "_"))
  calCols <- colnames(sequenceData)[str_detect(colnames(sequenceData), "^Cal_")]
  for(i in seq_along(calCols))
  {
    sequenceData[[calCols[i]]] <- as.numeric(sequenceData[[calCols[i]]])
  }

  # Read row data
  intensityData <- fread(
    file,
    sep              = "\t",
    nrows            = length(nFields) - dataGroupRow[4] - ncol(sequenceData),
    header           = TRUE,
    skip             = dataGroupRow[4] + ncol(sequenceData),
    integer64        = 'numeric',
    na.strings       = c("", "NA")
  )
  # Remove blank column between sample data and intensity data
  intensityData <- intensityData[, -sequenceHeaderColumnNumber, with = FALSE]

  # Give intensity data columns a name
  setnames(
    intensityData,
    seq.int(sequenceHeaderColumnNumber, ncol(intensityData)),
    paste(
      "SeqId",
      sequenceData$SeqId,
      sep = "."
    )
  )

  # As with sequence data, ensure correct datatypes
  # TODO: Waiting to hear from SomaLogic about which columns are compulsory,
  # and which are optional.  Update this next code chunk when we know.
  # Warnings are suppressed due to 'adding' non-existent columns as NULL
  # which does nothing.  (This is intentional.)
  suppressWarnings(intensityData[
    j = `:=`(
      PlateId           = factor(PlateId),
      SlideId           = factor(SlideId),
      SampleId          = factor(SampleId),
      SampleType        = if(exists("SampleType")) factor(SampleType) else NULL,
      SampleMatrix      = if(exists("SampleMatrix")) factor(SampleMatrix) else NULL,
      Barcode           = if(exists("Barcode")) factor(Barcode) else NULL,
      Barcode2d         = if(exists("Barcode2d")) factor(Barcode2d) else NULL,
      SampleNotes       = if(exists("SampleNotes")) factor(SampleNotes) else NULL,
      SampleDescription = if(exists("SampleDescription")) factor(SampleDescription) else NULL,
      TimePoint         = if(exists("TimePoint")) as.numeric(TimePoint) else NULL,
      ExtIdentifier     = factor(ExtIdentifier),
      SampleGroup       = if(exists("SampleGroup")) factor(SampleGroup) else NULL,
      SiteId            = if(exists("SiteId")) factor(SiteId) else NULL,
      #  SampleUniqueID    = factor(SampleUniqueID), # no longer present in soma file version 1.2
      Subject_ID        = if(exists("Subject_ID")) factor(Subject_ID) else NULL,
      RowCheck          = factor(RowCheck)
    )
  ])


  # Remove failures
  if(keepOnlyPasses)
  {
    okSeqColumns <- isPass(sequenceData$ColCheck)
    sequenceData <- sequenceData[okSeqColumns, ]
    okColumns <- !str_detect(colnames(intensityData), "^SeqId\\.") #metadata
    okColumns[!okColumns] <- okSeqColumns
    intensityData <- intensityData[
      isPass(intensityData$RowCheck),
      okColumns,
      with = FALSE
    ]
  }

  setkey(sequenceData, SeqId)
  setkey(intensityData, SampleId)

  # Return everything
#   structure(
#     intensityData,
#     SequenceInfo = sequenceData,
#     Metadata     = header,
#     Checksum     = checksum,
#     class        = c("WideSomaLogicData", "data.table", "data.frame")
#   )
  setattr(intensityData, "SequenceInfo", sequenceData)
  setattr(intensityData, "Metadata", header)
  setattr(intensityData, "Checksum", checksum)
  setattr(intensityData, "class", c("WideSomaLogicData", "data.table", "data.frame"))
  intensityData
}

#' @rdname readAdat
#' @export
readSomaLogic <- function(file, keepOnlyPasses = TRUE, dateFormat = "%d/%m/%Y")
{
  .Deprecated("readAdat")
  readAdat(file, keepOnlyPasses, dateFormat)
}

#' Is the value a pass
#'
#' Checks if a string is "PASS".
#' @param x A character vector or factor.
#' @return A logical vector, the same length as the input, which is \code{TRUE}
#' whenever the input is the string "PASS".  Missing values return \code{FALSE}.
#' @author Richard Cotton
isPass <- function(x)
{
  !is.na(x) & x == "PASS"
}

#' Melt a WideSomaLogicData object
#'
#' Convert a \code{WideSomaLogicData} object from wide format to long format.
#' @param x An object of class \code{WideSomaLogicData}.
#' @param ... Currently unused.
#' @return An object of class \code{LongSomaLogicData} that inherits from
#' \code{data.frame}.
#' This function melts the sample data contained in a \code{WideSomaLogicData}
#' object so the sequence IDs are contained in a single column \code{SeqID},
#' with the corresponding intensities in a single column named \code{Intensity}.
#' the \code{SequenceInfo} attribute of the input is then merged into this.
#' the \code{Metadata} and \code{Checksum} attributes are preserved.
#' @importFrom data.table setkey
#' @importFrom reshape2 melt
#' @importFrom stringr str_detect
#' @export
#' @author Richard Cotton
melt.WideSomaLogicData <- function(x, ...)
{
  isSeqColumn <- str_detect(colnames(x), "^SeqId\\.")
  class(x) <- c("data.table", "data.frame")
  long <- melt(
    x,
    id.vars       = colnames(x)[!isSeqColumn],
    measure.vars  = colnames(x)[isSeqColumn],
    variable.name = "SeqId",
    value.name    = "Intensity"
  )
  long$SeqId <- substring(long$SeqId, 7)
  setkeyv(long, "SeqId")
  long <- long[attr(x, "SequenceInfo")]
  structure(
    long,
    Metadata     = attr(x, "Metadata"),
    Checksum     = attr(x, "Checksum"),
    class        = c("LongSomaLogicData", "data.table", "data.frame")
  )
}

#' Convert an object into a molten data frame
#'
#' See \code{\link[reshape2]{melt}}.
#' @name melt
#' @export
melt <- reshape2::melt

#' Get the intensities from a WideSomaLogicData or LongSomaLogicData object
#'
#' Gets the intensities from an object of class \code{WideSomaLogicData} or
#' \code{LongSomaLogicData}.
#' @param x An object of class \code{WideSomaLogicData} or
#' \code{LongSomaLogicData}.
#' @param ... Variables passed from other methods. Currently
#' ignored.
#' @return A numeric matrix of intensities for each protein. Row names are taken
#' from the \code{SampleId} of the input.  Column names are the protein
#' sequence IDs.
#' @examples
#' \donttest{
#' unzip(
#'   system.file("extdata", "soma_atkin_diabetes.zip", package = "pdapmain"),
#'   exdir = tempdir()
#' )
#' soma_file <- file.path(tempdir(), "soma_atkin_diabetes.adat")
#' wide_soma_data <- readSomaLogic(soma_file)
#' unlink(soma_file)
#' getIntensities(wide_soma_data)         # A matrix
#' long_soma_data <- melt(wide_soma_data)
#' getIntensities(long_soma_data)         # A data.table
#' }
#' @importFrom stringr str_detect
#' @export
#' @author Richard Cotton
getIntensities <- function(x, ...)
{
  UseMethod("getIntensities")
}

#' @export
getIntensities.WideSomaLogicData <- function(x, ...)
{
  isSeqColumn <- str_detect(colnames(x), "^SeqId\\.")
  class(x) <- c("data.table", "data.frame")
  m <- as.matrix(x[, isSeqColumn, with = FALSE])
  rownames(m) <- x$SampleId
  colnames(m) <- substring(colnames(m), 7)
  m
}

#' @export
getIntensities.LongSomaLogicData <- function(x, ...)
{
  class(x) <- c("data.table", "data.frame")
  x[, list(SeqId, SampleId, Intensity)]
}

#' @rdname getIntensities
as.matrix.WideSomaLogicData <- function(x, ...)
{
  .Deprecated("getIntensities")
  getIntensities(x, ...)
}

#' Get WideSomaLogicData attributes
#'
#' Shortcut functions to get attributes of \code{WideSomaLogicData} objects.
#' @param x An object of class \code{WideSomaLogicData}.
#' @param value Value to set the attribute to.
#' @return An attribute of the input.
#' For inputs that are not \code{WideSomaLogicData} objects, the return value
#' may be \code{NULL}.
#' @seealso \code{\link{readSomaLogic}}
#' @name WideSomaLogicDataAttributes
NULL

#' @rdname WideSomaLogicDataAttributes
#' @export
getSequenceInfo <- function(x)
{
  attr(x, "SequenceInfo", exact = TRUE)
}

#' @rdname WideSomaLogicDataAttributes
#' @export
getMetadata <- function(x)
{
  attr(x, "Metadata", exact = TRUE)
}

#' @rdname WideSomaLogicDataAttributes
#' @export
getChecksum <- function(x)
{
  attr(x, "Checksum", exact = TRUE)
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom assertive assert_is_inherited_from
#' @export
setSequenceInfo <- function(x, value)
{
  assert_is_inherited_from(value, "data.table")
  attr(x, "SequenceInfo") <- value
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom assertive assert_is_list
#' @export
setMetadata <- function(x, value)
{
  assert_is_list(value)
  attr(x, "Metadata") <- value
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom assertive assert_is_character
#' @export
setChecksum <- function(x, value)
{
  assert_is_character(value)
  attr(x, "Checksum") <- value
}

#' Indexing for WideSomaLogicData objects
#'
#' Wrapper to \code{[.data.table}, ensuring that the \code{SequenceInfo},
#' \code{Metadata} and \code{Checksum} attributes are preserved.
#' @param x A \code{WideSomaLogicData} object.
#' @param ... Passed to \code{[.data.table}.
#' @return A \code{WideSomaLogicData} object.
#' @seealso \code{\link[data.table]{data.table}}
#' @export
`[.WideSomaLogicData` <- function(x, ...)
{
  sequenceInfo <- getSequenceInfo(x)
  metadata     <- getMetadata(x)
  checksum     <- getChecksum(x)
  class(x) <- c("data.table", "data.frame")
  structure(
    x[...],
    SequenceInfo = sequenceInfo,
    Metadata     = metadata,
    Checksum     = checksum,
    class        = c("WideSomaLogicData", "data.table", "data.frame")
  )
}
