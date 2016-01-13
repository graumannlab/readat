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
#' @importFrom data.table setDT
#' @importFrom data.table setkeyv
#' @importFrom data.table melt.data.table
#' @importFrom stringi stri_detect_regex
#' @export
#' @author Richard Cotton
melt.WideSomaLogicData <- function(data, ..., na.rm = FALSE, value.name = "Intensity")
{
  isSeqColumn <- stri_detect_regex(colnames(data), "^SeqId\\.")
  #class(data) <- c("data.table", "data.frame")
  setDT(data)
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
#' @export
getIntensities.WideSomaLogicData <- function(x, rowsContain = c("samples", "sequences"), reorder = FALSE, ...)
{
  rowsContain <- match.arg(rowsContain)
  isSeqColumn <- stri_detect_regex(colnames(x), "^SeqId\\.")
  class(x) <- c("data.table", "data.frame")
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
#' @export
getIntensities.LongSomaLogicData <- function(x, ...)
{
  class(x) <- c("data.table", "data.frame")
  x[, list(SeqId, ExtIdentifier, Intensity)]
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
#' @export
getSampleData.WideSomaLogicData <- function(x, ...)
{
  isSampleColumn <- !stri_detect_regex(colnames(x), "^SeqId\\.")
  class(x) <- c("data.table", "data.frame")
  x[, isSampleColumn, with = FALSE]
}

#' @rdname WideSomaLogicDataAttributes
#' @export
getSampleData.LongSomaLogicData <- function(x, ...)
{
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
#' @return A numeric matrix of intensities for each protein. Row names are taken
#' from the \code{ExtIdentifer} of the input.  Column names are the protein
#' sequence IDs.
#' @seealso \code{\link{readAdat}}
#' @examples
#' \donttest{
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
#' View(intWideSamp, "Wide intensities, samples per row")
#'
#' intWideSeq <- getIntensities(                 # The transpose
#'   wide_soma_data,
#'   rowsContain = "sequences"
#' )
#' View(intWideSeq, "Wide intensities, seqs per row")
#'
#' # Sample data is always a data table
#' sampWide <- getSampleData(wide_soma_data)
#' View(sampWide, "Wide sample data")
#'
#' if(requireNamespace("reshape2"))
#' {
#'   # For LongSomaLogicData objects, the intensities are returned
#'   # as a data.table
#'   long_soma_data <- reshape2::melt(wide_soma_data)
#'   intLong <- getIntensities(long_soma_data)
#'   View(intLong, "Long intensities")
#'
#'   # Sample data has a different shape now
#'   sampLong <- getSampleData(long_soma_data)
#'   View(sampLong, "Long sample data")
#' }
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
getSequenceInfo <- function(x)
{
  .Deprecated("getSequenceData")
  getSequenceData(x)
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

# calculateChecksum <- function(x)
# {
#   digest::digest(list(getMetadata(x), getSequenceData(x), getSampleData(x), getIntensities(x)), "sha1")
# }

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
#' @importFrom assertive.types assert_is_inherited_from
#' @export
setSequenceInfo <- function(x, value)
{
  .Deprecated("setSequenceData")
  setSequenceData(x, value)
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
    seqIdPrefix <- stri_detect_regex(colnames(value), "^SeqId\\.")
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
#' @param drop Should dimensions be dropped? Passed to \code{[.data.table}.
#' Note that this defaults to FALSE, unlike the data.tale method.
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
#' unlink(soma_file)
#' @importFrom data.table as.data.table
#' @export
`[.WideSomaLogicData` <- function(x, ..., drop = FALSE)
{
  sequenceData <- getSequenceData(x)
  metadata     <- getMetadata(x)
  checksum     <- getChecksum(x)
  # x <- as.data.table(x)
  setDT(x)
  y <- x[..., drop = drop]
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
