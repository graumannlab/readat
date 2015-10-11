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
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom Biobase exprs
#' @importFrom Biobase exprs<-
#' @export
soma2eset <- function(somaObj, log2Transform = TRUE){

   # assayData
   myIntensities <- t(getIntensities(somaObj))
   myMeta <- getMetadata(somaObj)

   # fData
   featureDF <- getSequenceData(somaObj)
   featureDF <- data.frame(featureDF)
   rownames(featureDF) <- featureDF$SeqId

   # pData
   sampleDF <- data.frame(
         soma_sample_id = somaObj$SampleId,
         soma_plate_id = somaObj$PlateId,
         soma_ext_id = somaObj$ExtIdentifier,
         row.names = somaObj$SampleId)

   # forge eset
   myEset <- ExpressionSet(myIntensities, AnnotatedDataFrame(sampleDF), AnnotatedDataFrame(featureDF))

   # log2 transform
   if (log2Transform){
      exprs(myEset) <- log2(exprs(myEset))
   }

   # return eset
   return(myEset)
}
