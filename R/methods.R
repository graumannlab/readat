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
#' the \code{SequenceData} attribute of the input is then merged into this.
#' the \code{Metadata} and \code{Checksum} attributes are preserved.
#' @importFrom data.table setkeyv
#' @importFrom data.table melt.data.table
#' @importFrom stringr str_detect
#' @export
#' @author Richard Cotton
melt.WideSomaLogicData <- function(x, ...)
{
  isSeqColumn <- str_detect(colnames(x), "^SeqId\\.")
  class(x) <- c("data.table", "data.frame")
  long <- suppressMessages(melt.data.table(
    x,
    id.vars       = colnames(x)[!isSeqColumn],
    measure.vars  = colnames(x)[isSeqColumn],
    variable.name = "SeqId",
    value.name    = "Intensity"
  ))
  long$SeqId <- substring(long$SeqId, 7)
  setkeyv(long, "SeqId")
  # long <- long[attr(x, "SequenceData")]
#   structure(
#     long,
#     Metadata     = attr(x, "Metadata"),
#     Checksum     = attr(x, "Checksum"),
#     class        = c("LongSomaLogicData", "data.table", "data.frame")
#   )
  setattr(long, "Metadata", attr(x, "Metadata"))
  setattr(long, "Checksum", attr(x, "Checksum"))
  setattr(long, "class", c("LongSomaLogicData", "data.table", "data.frame"))
  long
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
#' @importFrom assertive assert_is_inherited_from
#' @importFrom data.table setattr
#' @export
setSequenceData <- function(x, value)
{
  assert_is_inherited_from(value, "data.table")
  setattr(x, "SequenceData", value)
}

#' @rdname WideSomaLogicDataAttributes
#' @importFrom assertive assert_is_inherited_from
#' @export
setSequenceInfo <- function(x, value)
{
  .Deprecated("setSequenceData")
  setSequenceData(x, value)
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
#' Wrapper to \code{[.data.table}, ensuring that the \code{SequenceData},
#' \code{Metadata} and \code{Checksum} attributes are preserved.
#' @param x A \code{WideSomaLogicData} object.
#' @param ... Passed to \code{[.data.table}.
#' @return If the indexing returns a A \code{WideSomaLogicData} object.
#' @seealso \code{\link[data.table]{data.table}}
#' @importFrom data.table is.data.table
#' @examples
#' # Indexing returns a data.table, so the WideSomaLogicClass is preserved
#' sl1[1:5, list(`SeqId.3896-5_2`)]
#' # Indexing simplifies to a numeric vector, so the class is lost
#' sl1[1:5, `SeqId.3896-5_2`]
#' @export
`[.WideSomaLogicData` <- function(x, ...)
{
  sequenceData <- getSequenceData(x)
  metadata     <- getMetadata(x)
  checksum     <- getChecksum(x)
  class(x) <- c("data.table", "data.frame")
  y <- x[..., drop = FALSE]
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
