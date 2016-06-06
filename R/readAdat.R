#' @include sfread.R
NULL

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
#' @param keepOnlySamples A logical value indicating whether or not to keep
#' only the rows containing actual samples (as opposed to QC, buffer, and
#' calibrator samples).
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
#' somaFile <- extractSampleData()
#' wideSomaData <- readAdat(somaFile)
#' str(wideSomaData, list.len = 35)
#' unlink(somaFile)
#' @include checkIds.R
#' @importFrom assertive.files is_connection
#' @importFrom assertive.properties is_scalar
#' @importFrom assertive.types assert_is_a_string
#' @importFrom assertive.files assert_all_are_existing_files
#' @importFrom data.table ":="
#' @importFrom data.table as.data.table
#' @importFrom data.table fread
#' @importFrom data.table setattr
#' @importFrom data.table setkeyv
#' @importFrom data.table setnames
#' @export
#' @author Richard Cotton
readAdat <- function(file, keepOnlyPasses = TRUE, keepOnlySamples = TRUE,
  dateFormat = "%Y-%m-%d",  verbose = getOption("verbose"))
{
  # stri_read_lines and fread don't behave well with file connections
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

  # Remove QC, buffer, calibrator samples
  if(keepOnlySamples)
  {
    sampleAndIntensityData <- removeNonSamples(sampleAndIntensityData)
  }

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
  WideSomaLogicData(sampleAndIntensityData, sequenceData, metadata, checksum)
}

#' Create a WideSomaLogicData object
#'
#' Creates and object of class \code{WideSomaLogicData}.
#' @param sampleAndIntensityData A data.table of sample and intensity data.
#' @param sequenceData A data.table of sequence data.
#' @param metadata A list of metadata.
#' @param checksum A string containing a SHA1 checksum.
#' @return An object of class \code{WideSomaLogicData}.
#' @importFrom assertive.types assert_is_data.frame
#' @importFrom assertive.types assert_is_list
#' @importFrom assertive.types assert_is_a_string
#' @importFrom data.table is.data.table
#' @importFrom data.table copy
#' @importFrom data.table setattr
#' @examples
#' wsld <- WideSomaLogicData(
#'   data.frame(SampleStuff = letters, IntensityStuff = rnorm(26)),
#'   data.frame(SequenceStuff = LETTERS),
#'   list(MetadataStuff = Sys.Date())
#' )
#' str(wsld)
#' @export
WideSomaLogicData <- function(sampleAndIntensityData, sequenceData, metadata,
  checksum = paste0(rep.int("f", 40), collapse = ""))
{
  assert_is_data.frame(sampleAndIntensityData)
  sampleAndIntensityData <- if(is.data.table(sampleAndIntensityData))
  {
    copy(sampleAndIntensityData)
  } else # is.data.frame(sampleAndIntensityData)
  {
    as.data.table(sampleAndIntensityData)
  }

  assert_is_data.frame(sequenceData)
  sequenceData <- if(is.data.table(sequenceData))
  {
    copy(sequenceData)
  } else # is.data.frame(sequenceData)
  {
    as.data.table(sequenceData)
  }

  assert_is_list(metadata)
  assert_is_a_string(checksum)

  setattr(sampleAndIntensityData, "SequenceData", sequenceData)
  setattr(sampleAndIntensityData, "Metadata", metadata)
  setattr(sampleAndIntensityData, "Checksum", checksum)
  setattr(
    sampleAndIntensityData,
    "class",
    c("WideSomaLogicData", "data.table", "data.frame")
  )
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
#' @importFrom stringi stri_replace_first_regex
#' @importFrom stringi stri_detect_regex
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
  # V1 and V2 are column names autogenerated by fread
  # (since there are no headers)
  metadata <- with(
    metadata,
    setNames(as.list(V2), stri_replace_first_regex(V1, "^!", ""))
  )
  plateTestFields <- names(metadata)[
    stri_detect_regex(names(metadata), "(?:PlateMedianCal|PlateTailPercent)")
  ]
  metadata %>%
    updateFields("Version", as.package_version) %>%
    updateFields(
      c("CreatedDate", "ExpDate", "ProteinEffectiveDate"),
      as.Date,
      dateFormat
    ) %>%
    updateFields(plateTestFields, as.numeric)
}

#' @importFrom data.table fread
#' @importFrom data.table as.data.table
#' @importFrom magrittr %<>%
readSequenceData <- function(file, nSequenceFields, nSampleFields, skip,
  verbose = getOption("verbose"))
{
  # Can't use sfread here because data needs transposing before type checking
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
  cnames <- colnames(sequenceData)
  calCols <- cnames[stri_detect_regex(cnames, "^Cal_")]
  sequenceData %<>%
    updateFields(
      c("SeqId", "Target", "SomaId", "TargetFullName", "Organism", "Units"),
      factor
    ) %>%
    updateFields(
      c("UniProt", "EntrezGeneID", "EntrezGeneSymbol"),
      fixMultiValueSeparators
    ) %>%
    updateFields("ColCheck", fixFailFlag) %>%
    updateFields(c("CalReference", "Dilution", calCols), as.numeric)

  sequenceData
}


readSampleAndIntensityData <- function(file, nSequenceFields, nSampleFields,
  skip, seqIds, verbose = getOption("verbose"))
{
  # Read row data
  # Don't set nrows arg; just read to the end of the file.
  sampleAndIntensityData <- sfread(
    file,
    sep              = "\t",
    header           = TRUE,
    na.strings       = c("", "NA", "null"),
    stringsAsFactors = TRUE,
    skip             = skip,
    colClasses       = list(
      character = c(
        "ExtIdentifier", "SampleId", "SampleGroup", "SampleGroup",
        "SampleNotes", "AssayNotes", "PlateId", "SlideId",
        "ScannerID", "Subject_ID", "SiteId", "TubeUniqueID",
        "SsfExtId", "Barcode", "Barcode2d", "SampleUniqueId",
        "SampleType", "SampleMatrix", "SampleDescription", "PlatePosition",
        "RowCheck"
      ),
      numeric = c("HybControlNormScale", "NormScale_40", "NormScale_1", "NormScale_0_005"),
      integer = "Subarray"
    ),
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

  sampleAndIntensityData %>%
    updateFields("RowCheck", fixFailFlag)

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

removeNonSamples <- function(sampleAndIntensityData)
{
  if(is.null(sampleAndIntensityData$SampleType))
  {
    warning("There is no 'SampleType' field, so QC/buffer/calibrator samples cannot be removed.")
  } else
  {
    okSampleRows <- sampleAndIntensityData$SampleType == "Sample"
    if(!all(okSampleRows))
    {
      if(!all(okSampleRows))
      {
        nBadSamples <- sum(!okSampleRows)
        message(
          sprintf(
            ngettext(
              nBadSamples,
              "Removing %d QC/buffer/calibrator sample.",
              "Removing %d QC/buffer/calibrator samples.",
              domain = NA # don't translate, at least for now
            ),
            nBadSamples
          )
        )
      }
      sampleAndIntensityData <- sampleAndIntensityData[okSampleRows]
    }
  }
  sampleAndIntensityData
}

removeFailures <- function(sequenceData, sampleAndIntensityData)
{
  if(is.null(sequenceData$ColCheck))
  {
    warning("There is no 'ColCheck' field, so failing sequences cannot be removed.")
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
      sequenceData <- sequenceData[okSeqColumns, ]
      okColumns <- !colnamesStartWithSeqId(sampleAndIntensityData)
      okColumns[!okColumns] <- okSeqColumns
      sampleAndIntensityData <- sampleAndIntensityData[
        j = okColumns,
        with = FALSE
      ]
    }
  }
  if(is.null(sampleAndIntensityData$RowCheck))
  {
    warning("There is no 'RowCheck' field, so failing samples cannot be removed.")
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
      sampleAndIntensityData <- sampleAndIntensityData[okSampleRows]
    }
  }
  list(
    sequenceData = sequenceData,
    sampleAndIntensityData = sampleAndIntensityData
  )
}

#' Transmute a field, if it exists
#'
#' If a field in a data frame of data table exists, transform it and replace it.
#' @param data A data frame (or derivative) or a list.
#' @param fields A character vector of names of fields that may exist in
#' \code{data}.
#' @param transform A function or string naming a function.
#' @param ... Passed to \code{transform}.
#' @note data.table's := and dplyr's filter are a bit painful to use with
#' columns that may no exist in the data frame.
#' @noRd
updateFields <- function(data, fields, transform, ...)
{
  transform <- match.fun(transform)
  for(field in fields)
  {
    if(!is.null(data[[field]]))
    {
      data[[field]] <- transform(data[[field]], ...)
    }
  }
  data
}

#' @importFrom stringi stri_replace_all_regex
fixMultiValueSeparators <- function(x)
{
  x %>%
    as.character %>%
    stri_replace_all_regex("[, ]+", " ") %>%
    factor
}

#' Replace "FAIL" with "FLAG"
#'
#' SomaLogic are inconsistent about using "FAIL" and "FLAG" for row and column
#' checks.  This enforces the use of "FLAG".
#' @param x A character vector or factor of a row or column check.
#' @return A factor with two levels: "PASS" and "FLAG".
#' @examples
#' readat:::fixFailFlag(c("PASS", "FLAG", "FAIL", "something else", NA))
#' @noRd
fixFailFlag <- function(x)
{
  x <- as.character(x)
  if(any(x == "FAIL"))
  {
    warning("Changing 'FAIL' to 'FLAG'.")
    x[x == "FAIL"] <- "FLAG"
  }
  factor(x, levels = c("PASS", "FLAG"))
}

#' Does the input contain Sequence ID column names?
#'
#' Checks if the input begins with "SeqID.".
#' @param x A character vector.  Typically the column names of sequence data.
#' @return A logical vector with the same length as \code{x}.
#' @examples
#' somaFile <- extractSampleData()
#' wideSomaData <- readAdat(somaFile)
#' colnamesStartWithSeqId(wideSomaData)
#' unlink(somaFile)
#' @importFrom stringi stri_detect_regex
#' @export
colnamesStartWithSeqId <- function(x)
{
  stri_detect_regex(colnames(x), "^SeqId\\.")
}

#' Convert SeqIds to Aptamers
#'
#' Converts SomaLogic sequence IDs to aptamer IDs.
#' @param seqId A character vector or factor of sequence IDs.
#' @return A character vector of factors.
#' @importFrom stringi stri_split_fixed
#' @examples
#' convertSeqIdToAptamer("2717-12_5")
#' @export
convertSeqIdToAptamer <- function(seqId)
{
  stri_split_fixed(seqId, "_", n = 2, simplify = TRUE)[, 1L]
}
