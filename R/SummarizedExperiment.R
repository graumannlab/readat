#' Convert objects into SummarizedExperiments
#'
#' Convert objects into \code{SummarizedExperiments}.
#' @param x An object to transform.  Currently only \code{WideSomaLogicData}
#' objects are supported.
#' @param ... Arguments passed between methods. Currently unused.
#' @return An object of classs SummarizedExperiment.
#' @examples
#' somaFile <- extractSampleData()
#' wideSomaData <- readAdat(somaFile)
#' as.SummarizedExperiment(wideSomaData)
#' unlink(somaFile)
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
as.SummarizedExperiment <- function(x, ...)
{
  UseMethod("as.SummarizedExperiment")
}

#' @rdname as.SummarizedExperiment
#' @export
as.SummarizedExperiment.WideSomaLogicData <- function(x, ...)
{
  SummarizedExperiment(
    list(intensities = getIntensities(x, rowsContain = "sequences")),
    rowData = getSequenceData(x),
    colData = getSampleData(x),
    metadata = getMetadata
  )
}

