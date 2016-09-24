#' Get sequences with the largest between group variation
#'
#' Get the sequences with the largest between group variation in mean log
#' intensity.
#' @param x An object of class LongSomaLogicData.
#' @param n An integer of the number of sequences to return, passed to
#' \code{\link[dplyr]{top_n}}.
#' @param group A formula, string or quoted column name of the column that
#' defines the groups to split by.
#' @param ... Passed to and from methods, but currently unused.
#' @return A data table with \code{n} rows and the following columns.
#' \describe{
#'   \item{SeqId}{SomaLogic sequence identifier.}
#'   \item{VariationBetweenGroups}{The largest mean log intensity within a
#'   \code{group} divided by the smallest mean log intensity within a
#'   \code{group}.}
#' }
#' @export
getSequencesWithLargestBetweenGroupVariation <- function(x, n = 10,
  group = ~ SampleGroup, ...)
{
  UseMethod("getSequencesWithLargestBetweenGroupVariation")
}

#' @importFrom dplyr group_by_
#' @importFrom dplyr summarize_
#' @importFrom dplyr arrange_
#' @importFrom dplyr filter_
#' @importFrom dplyr desc
#' @importFrom utils head
#' @importFrom magrittr %>%
#' @export
getSequencesWithLargestBetweenGroupVariation.LongSomaLogicData <- function(x,
  n = 10, group = ~ SampleGroup, ...)
{
  sequenceData <- getSequenceData(x)
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

  # sequenceData[SeqId %in% bigBetweenGroupVariation$SeqId]
  sequenceData %>% filter_(~ SeqId %in% bigBetweenGroupVariation$SeqId)
}
