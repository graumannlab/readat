# Project: readat
#
# Author: Aditya Bhagwat
###############################################################################


#' Convert WideSomaLogicData into ExpressionSet
#' @param somaObj WideSomaLogicData object
#' @param log2Transform whether to log2 transform intensities or not
#' @return ExpressionSet object
#' @author Aditya Bhagwat
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase exprs exprs<-
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @export
soma2eset <- function(somaObj, log2Transform = TRUE){

  # assayData
  myIntensities <- getIntensities(somaObj, rowsContain = "sequences", reorder = TRUE)
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

  # log2 transform
  if (log2Transform){
    Biobase::exprs(myEset) <- log2(Biobase::exprs(myEset))
  }

  # return eset
  return(myEset)
}
