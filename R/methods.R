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
#' @importFrom data.table setkeyv
#' @importFrom data.table melt.data.table
#' @importFrom stringi stri_detect_regex
#' @export
#' @author Richard Cotton
melt.WideSomaLogicData <- function(data, ..., na.rm = FALSE, value.name = "Intensity")
{
  isSeqColumn <- stri_detect_regex(colnames(data), "^SeqId\\.")
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
  setattr(long, "class", c("LongSomaLogicData", "data.table", "data.frame"))
  long
}

#' Get the intensities or sample data
#'
#' Gets the intensities or sample data from an object of class
#' \code{WideSomaLogicData} or \code{LongSomaLogicData}.
#' @param x An object of class \code{WideSomaLogicData} or
#' \code{LongSomaLogicData}.
#' @param rowsContain Either samples or sequences.
#' @param ... Variables passed to and from from other methods.
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
#' wide_soma_data <- suppressWarnings(
#'   readSomaLogic(soma_file)
#' )
#' unlink(soma_file)
#'
#' intWideSamp <- getIntensities(wide_soma_data) # A matrix
#' View(intWideSamp, "Wide intensities, samples per row")
#'
#' intWideSeq <- getIntensities(                 # The transpose
#'   wide_soma_data,
#'   rowsContain = "sequences"
#' )
#' View(intWideSeq, "Wide intensities, seqs per row")
#'
#' sampWide <- getSampleData(wide_soma_data)     # A data.table
#' View(sampWide, "Wide sample data")
#'
#' if(requireNamespace("reshape2"))
#' {
#'   long_soma_data <- reshape2::melt(wide_soma_data)
#'   intLong <- getIntensities(long_soma_data)   # A data.table
#'   View(intLong, "Long intensities")
#'
#'   sampLong <- getSampleData(long_soma_data)   # A data.table
#'   View(sampLong, "Long sample data")
#' }
#' }
#' @importFrom stringi stri_detect_regex
#' @export
getIntensities <- function(x, ...)
{
  UseMethod("getIntensities")
}

#' @rdname getIntensities
#' @export
getIntensities.WideSomaLogicData <- function(x, rowsContain = c("samples", "sequences"), ...)
{
  rowsContain <- match.arg(rowsContain)
  isSeqColumn <- stri_detect_regex(colnames(x), "^SeqId\\.")
  class(x) <- c("data.table", "data.frame")
  m <- as.matrix(x[, isSeqColumn, with = FALSE])
  rownames(m) <- x$SampleId
  colnames(m) <- substring(colnames(m), 7)
  if(rowsContain == "samples") m else t(m)
}

#' @rdname getIntensities
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

#' @rdname getIntensities
#' @export
getSampleData <- function(x, ...)
{
  UseMethod("getSampleData")
}

#' @rdname getIntensities
#' @export
getSampleData.WideSomaLogicData <- function(x, ...)
{
  isSampleColumn <- !stri_detect_regex(colnames(x), "^SeqId\\.")
  class(x) <- c("data.table", "data.frame")
  x[, isSampleColumn, with = FALSE]
}

#' @rdname getIntensities
#' @export
getSampleData.LongSomaLogicData <- function(x, ...)
{
  class(x) <- c("data.table", "data.frame")
  x[, -"Intensity", with = FALSE]
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

#' Indexing for WideSomaLogicData objects
#'
#' Wrapper to \code{[.data.table}, ensuring that the \code{SequenceData},
#' \code{Metadata} and \code{Checksum} attributes are preserved.
#' @param x A \code{WideSomaLogicData} object.
#' @param i Row index, passed to \code{[.data.table}.
#' @param ... Passed to \code{[.data.table}.
#' @param drop Should dimensions be dropped? Passed to \code{[.data.table}.
#' @return If the indexing returns a A \code{WideSomaLogicData} object.
#' @seealso \code{\link[data.table]{data.table}}
#' @importFrom data.table is.data.table
#' @examples
#' \donttest{
#' unzip(
#'   system.file("extdata", "soma_atkin_diabetes.zip", package = "koraproteomics"),
#'   exdir = tempdir()
#' )
#' soma_file <- file.path(tempdir(), "soma_atkin_diabetes.adat")
#' wide_soma_data <- readSomaLogic(soma_file)
#'
#' # Indexing returns a data.table, so the WideSomaLogicClass is preserved
#' wide_soma_data[1:5, list(`SeqId.3896-5_2`)]
#' # Indexing simplifies to a numeric vector, so the class is lost
#' wide_soma_data[1:5, `SeqId.3896-5_2`]
#' unlink(soma_file)
#' }
#' @export
`[.WideSomaLogicData` <- function(x, i, ..., drop = FALSE)
{
  sequenceData <- getSequenceData(x)
  metadata     <- getMetadata(x)
  checksum     <- getChecksum(x)
  class(x) <- c("data.table", "data.frame")
  y <- x[i, ..., drop = drop]
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
