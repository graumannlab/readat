#' @importFrom assertive.base assert_all_are_not_false
#' @importFrom assertive.numbers assert_all_are_less_than_or_equal_to
#' @importFrom data.table fread
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom stats setNames
read_controls <- function(file = "controls.csv")
{
  controls <- fread(file, sep = ",", nrows = 96, header = FALSE, na.strings = "")
  controls %<>%
    setNames(c("PlatePosition", "BarCode"))
  n_controls <- sum(!is.na(controls$BarCode))
  assert_all_are_less_than_or_equal_to(n_controls, 12, severity = "warning")
  assert_all_are_not_false(
    stri_detect_regex(controls$BarCode, "^I[0-9]{6}$"),
    severity = "warning"
  )
  controls
}

#' @importFrom dplyr select_
read_samples <- function(file = "samples.csv")
{
  samples <- fread(file, sep = ",", nrows = 96, header = FALSE, na.strings = "NO READ")
  samples %<>%
    setNames(c("PlatePosition", "SampleId", "AlternateSampleId")) %>%
    select_(~ PlatePosition, ~ SampleId)
  n_controls <- sum(is.na(samples$SampleId))
  assert_all_are_less_than_or_equal_to(n_controls, 12, severity = "warning")
  assert_all_are_not_false(
    stri_detect_regex(samples$SampleId, "^[0-9]{9}$"),
    severity = "warning"
  )
  samples
}

read_slides <- function(file = "slides.csv")
{
  slides <- fread(file, sep = ",", header = FALSE, colClasses = "character")
  slides %<>%
    setNames("SlideId")
  assert_all_are_less_than_or_equal_to(nrow(slides), 12, severity = "warning")
  assert_all_are_not_false(
    stri_detect_regex(samples$SampleId, "^[0-9]{12}$"),
    severity = "warning"
  )
  pp <- expand.grid(Subarray = LETTERS[1:8], Slide = 1:12) %$%
    paste0(Subarray, Slide)
  data.frame(
    PlatePosition = pp,
    SlideId = rep(slides$SlideId, each = 8)
  )
}

read_comments <- function(file = "comments.csv")
{
  comments <- fread(file, sep = ",", header = FALSE, na.strings = "", colClasses = "character")
  comments %<>%
    setNames(c("PlatePosition", "SampleNotes", "AssayNotes"))
  assert_all_are_less_than_or_equal_to(nrow(comments), 96, severity = "warning")
  not_na <- !is.na(comments$SampleNotes)
  comments$SampleNotes[not_na] <- comments$SampleNotes[not_na] %>%
    tolower %>%
    strsplit("[, ]+") %>%     # standardize separator
    vapply(paste, collapse = ", ", character(1))
  assert_all_are_not_false(
    comments$SampleNotes %in% c("red", "yellow", "turbid", "red, turbid", "yellow, turbid"),
    severity = "warning"
  )
  pp <- expand.grid(Subarray = LETTERS[1:8], Slide = 1:12) %$%
    paste0(Subarray, Slide)
  data.frame(PlatePosition = pp, stringsAsFactors = FALSE) %>%
    left_join(comments, by = "PlatePosition")
}



