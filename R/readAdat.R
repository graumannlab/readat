utils::globalVariables("SeqId")
utils::globalVariables("SomaId")
utils::globalVariables("Target")
utils::globalVariables("TargetFullName")
utils::globalVariables("UniProt")
utils::globalVariables("EntrezGeneID")
utils::globalVariables("EntrezGeneSymbol")
utils::globalVariables("Organism")
utils::globalVariables("Units")
utils::globalVariables("CalReference")
utils::globalVariables("AssayNotes")
utils::globalVariables("ScannerID")
utils::globalVariables("TubeUniqueID")
utils::globalVariables("SsfExtId")
utils::globalVariables("SampleUniqueId")
utils::globalVariables("Subarray")
utils::globalVariables("PlatePosition")
utils::globalVariables("ColCheck")
utils::globalVariables("Dilution")

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
utils::globalVariables("HybControlNormScale")
utils::globalVariables("NormScale_40")
utils::globalVariables("NormScale_1")
utils::globalVariables("NormScale_0_005")

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
#' @param verbose Logical value indicating whether (lots of) diagnostic messages
#' should be shown.
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
#'   system.file("extdata", "PLASMA.1.3k.20151030.adat.zip", package = "readat"),
#'   exdir = tempdir()
#' )
#' soma_file <- file.path(tempdir(), "PLASMA.1.3k.HybNorm.MedNorm.Cal.20151030.adat")
#' wide_soma_data <- readAdat(soma_file)
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
#' @importFrom stringi stri_detect_regex
#' @importFrom stringi stri_replace_all_regex
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
    l <- removeFailures(sequenceData, sampleAndIntensityData)
    sequenceData <- l$sequenceData
    sampleAndIntensityData <- l$sampleAndIntensityData
  }

  setkeyv(sequenceData, "SeqId")
  setkeyv(sampleAndIntensityData, "ExtIdentifier")

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
  metadata <- fread(
    file,
    sep        = "\t",
    nrows      = colDataRow - headerRow - 1,
    header     = FALSE,
    skip       = headerRow,
    integer64  = "numeric",
    na.strings = c("", "NA", "null")
  )
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

  # Update column types.
  # SeqId and Target are compulsory.
  # Other common fields are updated if they exist.
  # Therefore, only update columns which we currently use, leave others unchanged
  suppressWarnings(sequenceData[
    j = `:=`(
      SeqId            = factor(SeqId),
      Target           = factor(Target),
      SomaId           = if(exists("SomaId")) factor(SomaId),
      TargetFullName   = if(exists("TargetFullName")) factor(TargetFullName) else NULL,
      UniProt          = if(exists("UniProt")) factor(stri_replace_all_regex(UniProt, "[, ]+", " ")) else NULL,
      EntrezGeneID     = if(exists("EntrezGeneID")) factor(stri_replace_all_regex(EntrezGeneID, "[, ]+", " ")) else NULL,
      EntrezGeneSymbol = if(exists("EntrezGeneSymbol")) factor(stri_replace_all_regex(EntrezGeneSymbol, "[, ]+", " ")) else NULL,
      Organism         = if(exists("Organism")) factor(Organism) else NULL,
      Units            = if(exists("Units")) factor(Units) else NULL,
      ColCheck         = if(exists("ColCheck")) factor(ColCheck) else NULL,
      CalReference     = if(exists("CalReference")) as.numeric(CalReference) else NULL,
      Dilution         = if(exists("Dilution")) as.numeric(Dilution) else NULL
    )
  ])

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
  # Don't set nrows arg; just read to the end of the file.
  sampleAndIntensityData <- fread(
    file,
    sep              = "\t",
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
  # ExtIdentifier is the only compulsory field, though SampleId, TimePoint,
  # SampleGroup, SampleNotes and AssayNotes should also be included in all
  # files from SomaLogic.
  # TimePoint is sometimes numeric, sometimes character data, so don't try to
  # coerce that field.
  # Warnings are suppressed due to 'adding' non-existent columns as NULL
  # which does nothing.  (This is intentional.)
  suppressWarnings(sampleAndIntensityData[
    j = `:=`(
      # Compulsory, and SomaLogic-compulsory
      ExtIdentifier     = factor(ExtIdentifier),
      SampleId          = if(exists("SampleDescription")) factor(SampleId) else NULL,
      SampleGroup       = if(exists("SampleGroup")) factor(SampleGroup) else NULL,
      SampleNotes       = if(exists("SampleNotes")) factor(SampleNotes) else NULL,
      AssayNotes        = if(exists("AssayNotes")) factor(AssayNotes) else NULL,
      # Optional IDs
      PlateId           = if(exists("PlateId")) factor(PlateId) else NULL,
      SlideId           = if(exists("SlideId")) factor(SlideId) else NULL,
      ScannerID         = if(exists("ScannerID")) factor(ScannerID) else NULL,
      Subject_ID        = if(exists("Subject_ID")) factor(Subject_ID) else NULL,
      SiteId            = if(exists("SiteId")) factor(SiteId) else NULL,
      TubeUniqueID      = if(exists("TubeUniqueID")) factor(TubeUniqueID) else NULL,
      SsfExtId          = if(exists("SsfExtId")) factor(SsfExtId) else NULL,
      Barcode           = if(exists("Barcode")) factor(Barcode) else NULL,
      Barcode2d         = if(exists("Barcode2d")) factor(Barcode2d) else NULL,
      SampleUniqueId    = if(exists("SampleUniqueId")) factor(SampleUniqueId) else NULL,
      # Optional Experimental conditions
      SampleType        = if(exists("SampleType")) factor(SampleType) else NULL,
      SampleMatrix      = if(exists("SampleMatrix")) factor(SampleMatrix) else NULL,
      SampleDescription = if(exists("SampleDescription")) factor(SampleDescription) else NULL,
      Subarray          = if(exists("Subarray")) as.integer(Subarray) else NULL,
      PlatePosition     = if(exists("PlatePosition")) factor(PlatePosition) else NULL,
      # Normalization, calibration, QC
      HybControlNormScale  = if(exists("HybControlNormScale")) as.numeric(HybControlNormScale) else NULL,
      NormScale_40      = if(exists("NormScale_40")) as.numeric(NormScale_40) else NULL,
      NormScale_1       = if(exists("NormScale_1")) as.numeric(NormScale_1) else NULL,
      NormScale_0_005   = if(exists("NormScale_0_005")) as.numeric(NormScale_0_005) else NULL,
      RowCheck          = if(exists("RowCheck")) factor(RowCheck) else NULL
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

removeFailures <- function(sequenceData, sampleAndIntensityData)
{
  if(is.null(sequenceData$ColCheck))
  {
    warning("There is no 'ColCheck field', so failing sequences cannot be removed.")
  } else
  {
    okSeqColumns <- isPass(sequenceData$ColCheck)
    if(!all(okSeqColumns))
    {
      nBadSeqs <- sum(!okSeqColumns)
      message(
        sprintf(
          ngettext(
            nBadSeqs,
            "Removing %d sequence that failed QC.",
            "Removing %d sequences that failed QC.",
            domain = NA # don't translate, at least for now
          ),
          nBadSeqs
        )
      )
    }
    sequenceData <- sequenceData[okSeqColumns, ]
    okColumns <- !stri_detect_regex(
      colnames(sampleAndIntensityData),
      "^SeqId\\."
    )
    okColumns[!okColumns] <- okSeqColumns
    sampleAndIntensityData <- sampleAndIntensityData[
      j = okColumns,
      with = FALSE
    ]
  }
  if(is.null(sampleAndIntensityData$RowCheck))
  {
    warning("There is no 'RowCheck field', so failing samples cannot be removed.")
  } else
  {
    okSampleRows <- isPass(sampleAndIntensityData$RowCheck)
    if(!all(okSampleRows))
    {
      nBadSamples <- sum(!okSampleRows)
      message(
        sprintf(
          ngettext(
            nBadSamples,
            "Removing %d sample that failed QC.",
            "Removing %d samples that failed QC.",
            domain = NA # don't translate, at least for now
          ),
          nBadSamples
        )
      )
    }
    sampleAndIntensityData <- sampleAndIntensityData[okSampleRows]
  }
  list(
    sequenceData = sequenceData,
    sampleAndIntensityData = sampleAndIntensityData
  )
}
