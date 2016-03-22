#' Get sequences with the largest between group variation
#'
#' Get the sequences with the largest between group variation in mean log intensity.
#' @param x An object of class LongSomaLogicData.
#' @param n An integer of the number of sequences to return, passed to
#' \code{\link[dplyr]{top_n}}.
#' @param group A formula, string or quoted column name of the column that
#' defines the groups to split by.
#' @param ... Passed to and from methods, but currently unused.
#' @export
getSequencesWithLargestBetweenGroupVariation <- function(x, n = 10, group = ~ SampleGroup, ...)
{
  UseMethod("getSequencesWithLargestBetweenGroupVariation")
}

#' @importFrom dplyr group_by_
#' @importFrom dplyr summarize_
#' @importFrom dplyr arrange_
#' @importFrom dplyr desc
#' @importFrom utils head
#' @importFrom magrittr %>%
#' @export
getSequencesWithLargestBetweenGroupVariation.LongSomaLogicData <- function(x, n = 10, group = ~ SampleGroup, ...)
{
  intensityByGroup <- x %>%
    group_by_(~ SeqId, group) %>%
    summarize_(MeanLogIntensity = ~ mean(log(Intensity), na.rm = TRUE))

  bigBetweenGroupVariation <- intensityByGroup %>%
    group_by_(~ SeqId) %>%
    summarize_(
      VariationBetweenGroups = ~ max(MeanLogIntensity) / min(MeanLogIntensity)
    ) %>%
    arrange_(~ desc(VariationBetweenGroups)) %>%
    head(n)

  x[SeqId %in% bigBetweenGroupVariation$SeqId]
}
