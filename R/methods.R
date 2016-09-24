#' Melt a WideSomaLogicData object
#'
#' Convert a \code{WideSomaLogicData} object from wide format to long format.
#' @param data An object of class \code{WideSomaLogicData}.
#' @param ... Currently unused.
#' @param na.rm TODO
#' @param value.name TODO
#' @return An object of class \code{LongSomaLogicData} that inherits from
#' \code{data.frame}.
#' This function melts the sample data contained in a \code{WideSomaLogicData}
#' object so the sequence IDs are contained in a single column \code{SeqID},
#' with the corresponding intensities in a single column named \code{Intensity}.
#' the \code{SequenceData} attribute of the input is then merged into this.
#' the \code{Metadata} and \code{Checksum} attributes are preserved.
#' @importFrom data.table copy
#' @importFrom data.table setkeyv
#' @importFrom data.table melt.data.table
#' @importFrom reshape2 melt
#' @importFrom stringi stri_detect_regex
#' @export
#' @include readAdat.R
melt.WideSomaLogicData <- function(data, ..., na.rm = FALSE,
  value.name = "Intensity")
{
  isSeqColumn <- colnamesStartWithSeqId(data)
  data <- copy(data)
  class(data) <- c("data.table", "data.frame")

  long <- suppressMessages(data.table::melt.data.table(
    data,
    id.vars       = colnames(data)[!isSeqColumn],
    measure.vars  = colnames(data)[isSeqColumn],
    variable.name = "SeqId",
    value.name    = value.name
  ))
  long$SeqId <- substring(long$SeqId, 7)
  setkeyv(long, "SeqId")
  setattr(long, "Metadata", attr(data, "Metadata"))
  setattr(long, "Checksum", attr(data, "Checksum"))
  setattr(long, "SequenceData", attr(data, "SequenceData"))
  setattr(long, "class", c("LongSomaLogicData", "data.table", "data.frame"))
  long
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom stringi stri_detect_regex
#' @export
getIntensities <- function(x, ...)
{
  UseMethod("getIntensities")
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom data.table copy
#' @export
getIntensities.WideSomaLogicData <- function(x,
  rowsContain = c("samples", "sequences"), reorder = FALSE, ...)
{
  x <- copy(x)
  class(x) <- c("data.table", "data.frame")
  rowsContain <- match.arg(rowsContain)
  isSeqColumn <- colnamesStartWithSeqId(x)
  m <- as.matrix(x[, isSeqColumn, with = FALSE])
  rownames(m) <- x$ExtIdentifier
  colnames(m) <- substring(colnames(m), 7)
  if(reorder)
  {
    rowOrder <- order(x$ExtIdentifier)
    colOrder <- order(colnames(m))
    m <- m[rowOrder, colOrder]
  }
  if(rowsContain == "samples") m else t(m)
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom data.table copy
#' @export
getIntensities.LongSomaLogicData <- function(x, ...)
{
  x <- copy(x)
  class(x) <- c("data.table", "data.frame")
  x[, c("SeqId", "ExtIdentifier", "Intensity"), with = FALSE]
}

#' @rdname WideSomaLogicDataAttributes
as.matrix.WideSomaLogicData <- function(x, ...)
{
  .Deprecated("getIntensities")
  getIntensities(x, ...)
}

#' @rdname WideSomaLogicDataAttributes
#' @export
getSampleData <- function(x, ...)
{
  UseMethod("getSampleData")
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom data.table copy
#' @export
getSampleData.WideSomaLogicData <- function(x, ...)
{
  x <- copy(x)
  class(x) <- c("data.table", "data.frame")
  isSampleColumn <- !colnamesStartWithSeqId(x)
  x[, isSampleColumn, with = FALSE]
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom data.table copy
#' @export
getSampleData.LongSomaLogicData <- function(x, ...)
{
  x <- copy(x)
  class(x) <- c("data.table", "data.frame")
  x[, -"Intensity", with = FALSE]
}

#' Get WideSomaLogicData attributes
#'
#' Accessors and mutators (getters and setters) for objects of class
#' \code{WideSomaLogicData} or \code{LongSomaLogicData}.
#' @param x An object of class \code{WideSomaLogicData} or
#' \code{LongSomaLogicData}.
#' @param value Value to set the attribute to.
#' @param rowsContain Either samples or sequences.
#' @param reorder If \code{TRUE}, rows are reordered by \code{ExtIdentifier} and
#' columns are reordered by \code{SeqId}, alphabetically.
#' @param prependSeqIdToColNames Logical.  Should "SeqId." be prepended to the
#' column names of the intensities? If \code{NA}, auto-guess whether they
#' should be.
#' @param ... Variables passed to and from from other methods.
#' @return \code{getIntensities} returns a numeric matrix of intensities for
#' each protein. Row names are taken from the \code{ExtIdentifier} of the input.
#' Column names are the protein sequence IDs.  (If you set
#' \code{rowsContain = "sequences"}, then rows and columns and their names are
#' swapped.)
#'
#' \code{getSampleData} returns a data table containing the sample data.  That
#' is, the input without the intensities or attributes.  There is one compulsory
#' column:
#' \describe{
#' \item{ExternalId}{Character. A unique identifer for the sample.}
#' }
#' These columns are not compulsory for the file spec, but are always included by SomaLogic:
#' \describe{
#' \item{SampleId}{Character. The customer's original sample identifier. This is either the same as the subject identifier, or the calibrator and buffer separated by a hyphen. For example "1234" or "PPT-09".}
#' \item{TimePoint}{Character or numeric.  The time point identified by the customer.}
#' \item{SampleGroup}{Character. Cohort information provided by the customer.}
#' \item{SampleNotes}{Character. Lab comments about the sample condition.}
#' \item{AssayNotes}{Character. Lab comments about assay adverse events.}
#' }
#' Common optional columns:
#' \describe{
#' \item{PlateId}{Character. Identifier for the plate. For assay experiments,
#' this also functions as an experiment identifier.}
#' \item{SlideId}{Character. Agilent slide barcode.}
#' \item{Subarray}{Integer. Agilent subarray number from 1 to 8.}
#' \item{SampleType}{Character. Either "Sample, "QC", "Buffer" or "Calibrator".}
#' \item{BarCode}{Character. 1D Barcode of aliquot.}
#' \item{Barcode2d}{Character. 2D Barcode of aliquot.}
#' \item{SiteId}{Character. The laboratory that ran the experiment, applicable
#' if subcontracted.}
#' \item{SampleUniqueId}{Character.  Where multiple experiments are combined,
#' this column can be used to ensure a unique sample ID.}
#' \item{Subject_ID}{Character. Identifier for the person/animal/thing being
#' sampled.}
#' \item{HybControlNormScale}{Numeric. Hybridization control normalization scale
#' factor. For example, 1.23.}
#' \item{NormScale_40}{Numeric. Median normalization scale factor, for only the
#' sequences diluted to 40\% concentration. For example, 1.23.}
#' \item{NormScale_0_005}{Numeric. Median normalization scale factor, for only
#' the sequences diluted to 0.005\% concentration. For example, 1.23.}
#' \item{NormScale_1}{Numeric. Median normalization scale factor, for only the
#' sequences diluted to 1\% concentration. For example, 1.23.}
#' \item{PercentDilution}{Numeric. How much was the sample diluted? Either 40,
#' 1, or 0.005.}
#' \item{SampleMatrix}{Character. Either "EDTA-Plasma", "Sodium Citrate Plasma"
#' or "Serum". Applicable if different matrices are used with different samples
#' in the file, otherwise use the StudyMatrix metadata value.}
#' \item{SampleDescription}{Character. Free text describing the sample.}
#' \item{TubeUniqueID}{Character. A unique identifier for every sample tube,
#' assigned by the customer.}
#' \item{RowCheck}{Character. Either "PASS" or "FAIL". Did the sample pass
#' quality control checks?}
#' }
#'
#' \code{getSequenceData} returns a data table containing the sequence data.
#' "Sequence" can be considered synonymous with "feature" or "SOMAmer reagent".
#' The following columns are compulsory.
#' \describe{
#' \item{SeqId}{A unique identifer for the SOMAmer reagent.}
#' \item{Target}{ The unique name for the targeted proteins.}
#' }
#' The following columns are optional but should always be provided by
#' SomaLogic.
#' \describe{
#' \item{SomaId}{Character.  The SomaLogic identifier for the protein target.
#' For the 1129 and 1310 assays, there is a one-to-one correspondence between
#' SeqId and SomaId, but in theory there is a many-to-one correspondence.}
#' \item{TargetFullName}{Character.  A description of the proteins targeted by
#' the SOMAmer reagent, taken from UniProt.}
#' \item{UniProt}{Character. The UniProt identifiers for the targeted proteins.
#' Older versions of the file format may also contain the value 'Family',
#' referring to a whole family of proteins, though this behaviour is considered
#' deprecated.}
#' \item{EntrezGeneID}{Character. The Entrez Gene identifiers for the genes
#' corresponding to the targeted proteins.}
#' \item{EntrezGeneSymbol}{Character. The Entrez Gene symbols for the genes
#' corresponding to the targeted proteins.}
#' \item{Organism}{Character. What animal does the target protein refer to?}
#' \item{ColCheck}{Character. Either PASS or FAIL. Did the sequence pass quality
#' control checks?}
#' \item{Units}{Character.  Unit of measurement for the SOMAmer reagent
#' abundances. Should always be RFU.}
#' \item{Type}{Character. One of "Spuriomer", "Protein",
#' "Hybridization Control", "Non-human Protein".}
#' \item{Cal_*PlateId*}{Numeric.  Calibration scale factor for each plate. For
#' example, 1.23.}
#' \item{CalReference}{Numeric. Calibration reference intensity, in RFU. For
#' example, 543.21.}
#' \item{Dilution}{Numeric. Concentration of the diluted sample compared to the
#' neat sample, as a percentage. Either 40, 1, or 0.005.}
#' }
#'
#' \code{getCheckSum} and \code{getMetadata} return the checksum (a string) and
#' metadata (a list) attributes of the input.
#' @seealso \code{\link{readAdat}}
#' @examples
#' # Get the sample dataset
#' soma_file <- extractSampleData()
#' wide_soma_data <- readAdat(soma_file)
#'
#' # Access its components
#' checksum <- getChecksum(wide_soma_data)
#' metadata <- getMetadata(wide_soma_data)
#' sequenceData <- getSequenceData(wide_soma_data)
#' sampleData <- getSampleData(wide_soma_data)
#'
#' # Intensities of a WideSomaLogicData object are a matrix
#' intWideSamp <- getIntensities(wide_soma_data)
#' \donttest{
#' # Not run due to side effects of View
#'   View(intWideSamp, "Wide intensities, samples per row")
#' }
#' intWideSeq <- getIntensities(                 # The transpose
#'   wide_soma_data,
#'   rowsContain = "sequences"
#' )
#' \donttest{
#'   View(intWideSeq, "Wide intensities, seqs per row")
#' }
#'
#' # Sample data is always a data table
#' sampWide <- getSampleData(wide_soma_data)
#' \donttest{
#'   View(sampWide, "Wide sample data")
#' }
#'
#' # For LongSomaLogicData objects, the intensities are returned
#' # as a data.table
#' long_soma_data <- reshape2::melt(wide_soma_data)
#' intLong <- getIntensities(long_soma_data)
#' \donttest{
#'   View(intLong, "Long intensities")
#' }
#'
#' # Sample data has a different shape now
#' sampLong <- getSampleData(long_soma_data)
#' \donttest{
#'   View(sampLong, "Long sample data")
#' }
#' @name WideSomaLogicDataAttributes
NULL

#' @rdname WideSomaLogicDataAttributes
#' @export
getSequenceData <- function(x)
{
  attr(x, "SequenceData", exact = TRUE)
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
#' @importFrom assertive.types assert_is_inherited_from
#' @importFrom data.table setattr
#' @export
setSequenceData <- function(x, value)
{
  assert_is_inherited_from(value, "data.table")
  setattr(x, "SequenceData", value)
  invisible(x)
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom assertive.types assert_is_list
#' @export
setMetadata <- function(x, value)
{
  assert_is_list(value)
  attr(x, "Metadata") <- value
  invisible(x)
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom assertive.types assert_is_character
#' @export
setChecksum <- function(x, value)
{
  assert_is_character(value)
  attr(x, "Checksum") <- value
  invisible(x)
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom assertive.base assert_are_identical
#' @importFrom assertive.properties assert_is_not_null
#' @importFrom assertive.types assert_is_data.frame
#' @importFrom data.table as.data.table
#' @export
setSampleData <- function(x, value)
{
  assert_is_data.frame(value)
  value <- as.data.table(value)
  assert_is_not_null(value$ExtIdentifier)
  intensities <- getIntensities(x)
  assert_are_identical(nrow(intensities), nrow(value))
  # Can't use dplyr::bind_cols with matrices
  sampleAndIntensityData <- cbind(value, intensities)
  invisible(
    WideSomaLogicData(
      sampleAndIntensityData,
      getSequenceData(x),
      getMetadata(x),
      getChecksum(x)
    )
  )
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom assertive.base are_identical
#' @importFrom assertive.base assert_are_identical
#' @importFrom assertive.types assert_is_inherited_from
#' @export
setIntensities <- function(x, value, prependSeqIdToColNames = NA)
{
  assert_is_inherited_from(value, c("data.frame", "matrix"))
  sampleData <- getSampleData(x)
  assert_are_identical(nrow(sampleData), nrow(value))
  # Should column names be prefixed with "SeqId.", or do they have it already?
  if(is.na(prependSeqIdToColNames))
  {
    seqIdPrefix <- colnamesStartWithSeqId(value)
    prependSeqIdToColNames <- if(any(seqIdPrefix))
    {
      if(!all(seqIdPrefix))
      {
        stop("Some column names are prefixed with 'SeqId.' but not others.")
      }
      FALSE
    } else
    {
      TRUE
    }
  }
  if(prependSeqIdToColNames)
  {
    colnames(value) <- paste0("SeqId.", colnames(value))
  }

  # If rows match in x and value, we can just cbind, otherwise need to merge
  sampleAndIntensityData <- if(
    are_identical(as.character(sampleData$ExtIdentifier), rownames(value))
  )
  {
    cbind(sampleData, value)
  } else
  {
    valueDF <- data.frame(ExtIdentifier = rownames(value), value)
    merge(sampleData, valueDF, by = "ExtIdentifier")
  }
  invisible(
    WideSomaLogicData(
      sampleAndIntensityData,
      getSequenceData(x),
      getMetadata(x),
      getChecksum(x)
    )
  )
}

#' Indexing for WideSomaLogicData objects
#'
#' Wrapper to \code{[.data.table}, ensuring that the \code{SequenceData},
#' \code{Metadata} and \code{Checksum} attributes are preserved.
#' @param x A \code{WideSomaLogicData} object.
#' @param ... Passed to \code{[.data.table}.
#' @return If the indexing returns a A \code{WideSomaLogicData} object.
#' @seealso \code{\link[data.table]{data.table}}
#' @importFrom data.table is.data.table
#' @examples
#' soma_file <- extractSampleData()
#' wide_soma_data <- readAdat(soma_file)
#'
#' # Indexing returns a data.table, so the WideSomaLogicClass is preserved
#' wide_soma_data[1:5, list(`SeqId.3896-5_2`)]
#' # Indexing simplifies to a numeric vector, so the class is lost
#' wide_soma_data[1:5, `SeqId.3896-5_2`]
#'
#' # Ignore the intensity columns (as per getSampleData)
#' j <- !colnamesStartWithSeqId(wide_soma_data)
#' wide_soma_data[1:5, j, with = FALSE]
#' unlink(soma_file)
#' @importFrom data.table as.data.table
#' @importFrom assertive.properties is_empty
#' @export
`[.WideSomaLogicData` <- function(x, ...)
{
  sequenceData <- getSequenceData(x)
  metadata     <- getMetadata(x)
  checksum     <- getChecksum(x)
  x <- copy(x)
  y <- NextMethod("[")

  # result may be a data.table, or have been simplified to a vector
  if(is.data.table(y))
  {
    setattr(y, "SequenceData", sequenceData)
    setattr(y, "Metadata", metadata)
    setattr(y, "Checksum", checksum)
    setattr(y, "class", c("WideSomaLogicData", "data.table", "data.frame"))
    y
  } else
  {
    structure(
      y,
      SequenceData = sequenceData,
      Metadata     = metadata,
      Checksum     = checksum
    )
  }
}
