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
#' \item{SequenceData}{A data frame.}
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
#' @importFrom assertive is_connection
#' @importFrom assertive is_scalar
#' @importFrom assertive assert_is_a_string
#' @importFrom assertive assert_all_are_existing_files
#' @importFrom data.table ":="
#' @importFrom data.table as.data.table
#' @importFrom data.table fread
#' @importFrom data.table setattr
#' @importFrom data.table setkeyv
#' @importFrom data.table setnames
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_replace
#' @export
#' @author Richard Cotton
readAdat <- function(file, keepOnlyPasses = TRUE, dateFormat = "%Y-%m-%d",
  verbose = getOption("verbose"))
{
  # stri_read_lines and fead don't behave well with file connections
  # Also, the logic gets complicated because the position in the file
  # keeps moving.
  if(is_connection(file))
  {
    file <- summary(file)$description
  }

  assert_is_a_string(file)
  assert_all_are_existing_files(file)

 # The file is split into several groups of data:
  # A checksum, header data, column data and row data.  Read each separately.

  # Opening and resaving in Excel appends TAB characters to make the data
  # rectangular.  This means that the number of fields is not reliable.
  # (Can't use count.fields.)
  firstChar <- readFirstChar(file)

  dataGroupRow <- which(firstChar == "^")

  if(length(dataGroupRow) < 4L)
  {
    stop(
      "The input file is malformed; there should be four rows begining with the ^ character but there are ",
      length(dataGroupRow)
    )
  }

  # Read SHA1 checksum
  # First line, without "!Checksum\t"
  checksum <- substring(readLines(file, 1), 11)


  # Read header
  metadata <- readMetadata(file, dataGroupRow[1], dataGroupRow[2], dateFormat)

  nSequenceFields <- getNFields(file, dataGroupRow[2], verbose = verbose)
  nSampleFields <- getNFields(file, dataGroupRow[3], verbose = verbose)

  # Read column data
  sequenceData <- readSequenceData(
    file,
    nSequenceFields,
    nSampleFields,
    dataGroupRow[4],
    verbose = verbose
  )

  # Check for bad UniProt and EntrezGene Ids/symbols.
  checkUniprotIds(sequenceData)
  checkEntrezGeneIds(sequenceData)
  checkEntrezGeneSymbols(sequenceData)

  sampleAndIntensityData <- readSampleAndIntensityData(
    file,
    nSequenceFields,
    nSampleFields,
    dataGroupRow[4] + nSequenceFields,
    sequenceData$SeqId,
    verbose = verbose
  )


  # Remove failures
  if(keepOnlyPasses)
  {
    okSeqColumns <- isPass(sequenceData$ColCheck)
    sequenceData <- sequenceData[okSeqColumns, ]
    okColumns <- !str_detect(colnames(sampleAndIntensityData), "^SeqId\\.") #metadata
    okColumns[!okColumns] <- okSeqColumns
    sampleAndIntensityData <- sampleAndIntensityData[
      isPass(sampleAndIntensityData$RowCheck),
      okColumns,
      with = FALSE
    ]
  }

  setkeyv(sequenceData, "SeqId")
  setkeyv(sampleAndIntensityData, "SampleId")

  # Return everything
  setattr(sampleAndIntensityData, "SequenceData", sequenceData)
  setattr(sampleAndIntensityData, "Metadata", metadata)
  setattr(sampleAndIntensityData, "Checksum", checksum)
  setattr(sampleAndIntensityData, "class", c("WideSomaLogicData", "data.table", "data.frame"))
  sampleAndIntensityData
}

#' Read the first character of each line
#'
#' Reads the first character of each line of a text file.
#' @param file A string or readable connection.
#' @return A character vector of single characters.
#' @references See \url{http://stackoverflow.com/q/27747426/134830}
#' @importFrom stringi stri_read_lines
#' @importFrom stringi stri_sub
readFirstChar <- function(file)
{
  # substring(readLines(file), 1, 1)
  stri_sub(stri_read_lines(file), 1, 1)
}

#' @importFrom data.table fread
readMetadata <- function(file, headerRow, colDataRow, dateFormat)
{
  # (Ab)using fread for this tends to results in unnecessary warnings.
  metadata <- suppressWarnings(fread(
    file,
    sep        = "\t",
    nrows      = colDataRow - headerRow - 1,
    header     = FALSE,
    skip       = headerRow,
    integer64  = "numeric",
    na.strings = c("", "NA", "null")
  ))
  metadata <- with(metadata, setNames(as.list(V2), substring(V1, 2)))
  within(
    metadata,
    {
      Version <- as.package_version(Version)
      CreatedDate <- as.Date(CreatedDate, format = dateFormat)
      ExpDate <- as.Date(ExpDate, format = dateFormat)
      ProteinEffectiveDate <- as.Date(ProteinEffectiveDate, format = dateFormat)
    }
  )
}

#' @importFrom data.table fread
#' @importFrom data.table as.data.table
#' @importFrom stringi stri_replace_all_regex
#' @importFrom stringi stri_detect_regex
readSequenceData <- function(file, nSequenceFields, nSampleFields, skip,
  verbose = getOption("verbose"))
{
  sequenceData <- fread(
    file,
    sep              = "\t",
    nrows            = nSequenceFields,
    colClasses       = "character",
    skip             = skip,
    header           = FALSE,
    stringsAsFactors = FALSE,
    integer64        = 'numeric',
    na.strings       = c("", "NA", "null"),
    verbose          = verbose
  )
  # Get the column that contains the headers
  sequenceHeaderColumnNumber <- nSampleFields + 1
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
      UniProt          = factor(stri_replace_all_regex(UniProt, "[, ]+", " ")),
      EntrezGeneID     = factor(stri_replace_all_regex(EntrezGeneID, "[, ]+", " ")),
      EntrezGeneSymbol = factor(stri_replace_all_regex(EntrezGeneSymbol, "[, ]+", " ")),
      Organism         = factor(Organism),
      Units            = factor(Units),
      ColCheck         = factor(ColCheck),
      CalReference     = as.numeric(CalReference),
      Dilution         = as.numeric(Dilution)
    )
  ]

  # There are some more columns that need fixing, which should have the names
  # sprintf("Cal_%s", str_replace(levels(intensityData$PlateId), " ", "_"))
  calCols <- colnames(sequenceData)[stri_detect_regex(colnames(sequenceData), "^Cal_")]
  for(i in seq_along(calCols))
  {
    sequenceData[[calCols[i]]] <- as.numeric(sequenceData[[calCols[i]]])
  }

  sequenceData
}


readSampleAndIntensityData <- function(file, nSequenceFields, nSampleFields, skip,
  seqIds, verbose = getOption("verbose"))
{
  # Read row data
  sampleAndIntensityData <- fread(
    file,
    sep              = "\t",
    nrows            = nSampleFields,
    header           = TRUE,
    skip             = skip,
    integer64        = 'numeric',
    na.strings       = c("", "NA", "null"),
    verbose          = verbose
  )
  # Remove blank column between sample data and intensity data
  sequenceHeaderColumnNumber <- nSampleFields + 1
  sampleAndIntensityData <- sampleAndIntensityData[
    j = -sequenceHeaderColumnNumber,
    with = FALSE
  ]

  # Give intensity data columns a name
  setnames(
    sampleAndIntensityData,
    seq.int(sequenceHeaderColumnNumber, ncol(sampleAndIntensityData)),
    paste0("SeqId.", seqIds)
  )

  # As with sequence data, ensure correct datatypes
  # TODO: Waiting to hear from SomaLogic about which columns are compulsory,
  # and which are optional.  Update this next code chunk when we know.
  # Warnings are suppressed due to 'adding' non-existent columns as NULL
  # which does nothing.  (This is intentional.)
  suppressWarnings(sampleAndIntensityData[
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
      # TimePoint is sometimes numeric, sometimes character data.
      # TimePoint         = if(exists("TimePoint")) as.numeric(TimePoint) else NULL,
      ExtIdentifier     = factor(ExtIdentifier),
      SampleGroup       = if(exists("SampleGroup")) factor(SampleGroup) else NULL,
      SiteId            = if(exists("SiteId")) factor(SiteId) else NULL,
      #  SampleUniqueID    = factor(SampleUniqueID), # no longer present in soma file version 1.2
      Subject_ID        = if(exists("Subject_ID")) factor(Subject_ID) else NULL,
      RowCheck          = factor(RowCheck)
    )
  ])
  sampleAndIntensityData
}

getNFields <- function(file, skip, verbose = getOption("verbose"))
{
  y <- scan(
    file,
    character(),
    sep = "\t",
    skip = skip,
    nlines = 1L,
    quiet = !verbose
  )
  as.integer(sum(nzchar(y))) - 1L
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

