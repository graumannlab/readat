# Project: somalogic
# 
# Author: Aditya Bhagwat
###############################################################################


#' Convert WideSomaLogicData into ExpressionSet
#' @param somaObj WideSomaLogicData object
#' @param log2Transform whether to log2 transform intensities or not
#' @return ExpressionSet object
#' @author Aditya Bhagwat
#' @import Biobase
#' @export
soma2eset <- function(somaObj, log2Transform = TRUE){
   
   # assayData
   myIntensities <- t(getIntensities(somaObj))
   myMeta <- getMetadata(somaObj)
   
   # fData
   featureDF <- getSequenceInfo(somaObj)
   featureDF <- data.frame(featureDF)
   rownames(featureDF) <- featureDF$SeqId
   
   # pData
   sampleDF <- data.frame(
         soma_sample_id = somaObj$SampleId, 
         soma_plate_id = somaObj$PlateId, 
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
