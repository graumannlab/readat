# Project: readat
#
# Author: Aditya Bhagwat
###############################################################################


#' Convert objects into ExpressionSets or MSnSets.
#'
#' Converts objects into \code{ExpressionSet} or \code{MSnSet} instances.
#' @param x An object to transform.  Currently only \code{WideSomaLogicData}
#' objects are supported.
#' @param somaObj A \code{WideSomaLogicData} object to transform.
#' @param log2Transform whether to log2 transform intensities or not
#' @param ... Passed between methods.
#' @return ExpressionSet or MSnSet object
#' @author Aditya Bhagwat
#' @examples
#' somaFile <- extractSampleData()
#' wideSomaData <- readAdat(somaFile)
#' as.ExpressionSet(wideSomaData)
#' if(requireNamespace("MSnbase"))
#' {
#'   as.MSnSet(wideSomaData)
#' }
#'
#' unlink(somaFile)
#' @importFrom Biobase annotation<-
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase exprs exprs<-
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @export
as.ExpressionSet <- function(x, log2Transform = TRUE, ...)
{
  UseMethod("as.ExpressionSet")
}

#' @rdname as.ExpressionSet
#' @export
as.ExpressionSet.WideSomaLogicData <- function(x, log2Transform = TRUE, ...)
{
  soma2eset(x, log2Transform = log2Transform)
}

#' @rdname as.ExpressionSet
#' @importFrom Biobase   annotation<-   experimentData<-   preproc<-
#' @export
soma2eset <- function(somaObj, log2Transform = TRUE){

  # assayData
  myIntensities <- getIntensities(
    somaObj, rowsContain = "sequences", reorder = TRUE
  )
  myMeta <- getMetadata(somaObj)

  # fData
  featureDF <- getSequenceData(somaObj)
  featureDF <- data.frame(featureDF, row.names = featureDF$SeqId)

  # pData
  sampleDF <- getSampleData(somaObj)
  sampleDF <- data.frame(sampleDF, row.names = somaObj$ExtIdentifier)

  # forge eset
  myEset <- ExpressionSet(myIntensities)
  Biobase::pData(myEset) <- sampleDF
  Biobase::fData(myEset) <- featureDF
  Biobase::annotation(myEset) <- 'somascan'

  # log2 transform
  if (log2Transform){
    Biobase::exprs(myEset) <- log2(Biobase::exprs(myEset))
  }

  # Add preprocessing and annotation info
  parameters <- attributes(somaObj)$Metadata
  Biobase::preproc(Biobase::experimentData(myEset)) <-
    list(assay    = 'somascan',
         entity   = 'epitope',
         quantity = 'abundance',
         software = 'somalogic',
         parameters = parameters
         )
  Biobase::annotation(myEset) <- parameters$StudyOrganism

  # return eset
  return(myEset)
}

#' Is object a soma eset?
#' @param esetObj eSet
#' @return logical
#' @importFrom Biobase   experimentData   preproc
#' @export
isSomaEset <- function(esetObj){
   Biobase::preproc(Biobase::experimentData(esetObj))$assay == 'somascan'
}

#' Get name of sample id variable
#' @param somaEset eset with soma data
#' @return character
#' @export
getSampleIdVar <- function(somaEset){
   'SampleId'
}

#' @rdname as.ExpressionSet
#' @export
as.MSnSet.WideSomaLogicData <- function(x, log2Transform = FALSE, ...)
{
  if(!requireNamespace("MSnbase", quietly = TRUE))
  {
    stop('MSnbase is not available; try running\nsource("https://bioconductor.org/biocLite.R")\nbiocLite("MSnbase")')
  }
  e <- as.ExpressionSet(x, log2Transform = log2Transform)
  MSnbase::as.MSnSet.ExpressionSet(e)

}

#' @rdname as.ExpressionSet
#' @export
as.MSnSet <- function(x, log2Transform = FALSE, ...)
{
  UseMethod("as.MSnSet")
}
