#' Extract the sample datasets
#'
#' Extracts the sample datasets from the ZIP files in the \code{extdata}
#' directory.
#' @param somascanMenuVersion String describing which data file to extract.
#' Either "1.3k" or "1.1k".
#' @param intermediateDir A path to a directory to extract the ADATY file into.
#' @return A path to the extracted ADAT file.
#' @examples
#' extractSampleData() ##1.3k dataset by default
#' extractSampleData("1.1k")
#' @importFrom utils unzip
#' @export
extractSampleData <- function(somascanMenuVersion = c("1.3k", "1.1k"),
  intermediateDir = tempdir())
{
  somascanMenuVersion <- match.arg(somascanMenuVersion)
  if(!file.exists(intermediateDir))
  {
    dir.create(intermediateDir, recursive = TRUE)
  }
  zipFile <- system.file(
    "extdata",
    paste("PLASMA", somascanMenuVersion, "20151030.adat.zip", sep = "."),
    package = "readat"
  )
  unzip(zipFile, exdir = intermediateDir)
  file.path(
    intermediateDir,
    paste(
      "PLASMA",
      somascanMenuVersion,
      "HybNorm.MedNorm.Cal.20151030.adat",
      sep = "."
    )
  )
}
