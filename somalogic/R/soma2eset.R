# Project: somalogic
# 
# Author: Aditya Bhagwat
###############################################################################


#' Convert WideSomaLogicData into ExpressionSet
#' @param somaObj WideSomaLogicData object
#' @param sampleDF dataframe with sample annotations
#' @param log2Transform whether to log2 transform intensities or not
#' @return ExpressionSet object
#' @author Aditya Bhagwat
#' @import Biobase
#' @export
soma2eset <- function(somaObj, 
      sampleDF = data.frame(sampleId = somaObj$SampleId, row.names = somaObj$SampleId), 
      log2Transform = TRUE){
   
   # assayData
   myIntensities <- t(getIntensities(somaObj))
   myMeta <- getMetadata(somaObj)
   
   # fData
   featureDF <- getSequenceInfo(somaObj)
   featureDF <- data.frame(featureDF)
   rownames(featureDF) <- featureDF$SeqId
   
   # forge eset
   myEset <- ExpressionSet(myIntensities, AnnotatedDataFrame(sampleDF), AnnotatedDataFrame(featureDF))
   
   # log2 transform
   if (log2Transform){
      exprs(myEset) <- log2(exprs(myEset))
   }
   
   # return eset
   return(myEset)
}
