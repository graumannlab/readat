#' @importFrom dplyr mutate_
PLATE_POSITIONS <- expand.grid(Subarray = 1:8, Slide = 1:12) %>%
    mutate_(PlatePosition = ~ paste0(LETTERS[Subarray], Slide))

#' Read SomaLogic Sample Submission Slides File
#'
#' @param file A string denoting the path to an input CSV file.  See Input file
#' specification section.
#' @return A \code{data.table} with 96 rows and 5 columns.
#' \describe{
#' \item{SampleNumber}{The integers 1 to 96.}
#' \item{SlideId}{The slide IDs from the input file.}
#' \item{Subarray}{Integers from 1 to 8 denoting the sample position for the
#' slide.}
#' \item{PlatePosition}{A letter followed by a number, constructed from the
#' Subarray ("A" for 1, "B" for 2, etc.) and the slide number from 1 to 12.}
#' \item{PercentDilution}{Always 40.}
#' }
#' @section Input file specification:
#' A CSV file without a header line containing up to twelve rows and one
#' column as follows.
#' \enumerate{
#' \item{Slide IDs, each 12 digits long.}
#' }
#' @seealso \code{\link{readControls}}, \code{\link{readComments}}, and
#' \code{\link{readSamples}} for reading other submission forms and
#' \code{\link{writeSampleSubmissionForm}} for usage examples.
#' @importFrom assertive.base assert_all_are_not_false
#' @importFrom data.table fread
#' @importFrom data.table data.table
#' @importFrom stats setNames
#' @export
readSlides <- function(file = "slides.csv")
{
  slides <- fread(file, sep = ",", header = FALSE, colClasses = "character")[[1]]
  length(slides) <- 12
  assert_all_are_not_false(
    stri_detect_regex(slides, "^[0-9]{12}$"),
    severity = "warning"
  )
  data.table(
    SampleNumber = 1:96,
    SlideId = rep(slides, each = 8),
    Subarray = PLATE_POSITIONS$Subarray,
    PlatePosition = PLATE_POSITIONS$PlatePosition,
    PercentDilution = 40
  )
}

#' Read SomaLogic Sample Submission Slides File
#'
#' @param file A string denoting the path to an input CSV file.  See Input file
#' specification section.
#' @return A \code{data.table} with 96 rows and 2 columns.
#' \describe{
#' \item{PlatePosition}{A letter followed by a number, constructed from the
#' Subarray ("A" for 1, "B" for 2, etc.) and the slide number from 1 to 12.}
#' \item{BarCode}{Sample barcode for QC, Calibrator, and Buffer samples, in the
#' form "I" followed by 6 digits.}
#' }
#' @section Input file specification:
#' A CSV file without a header line containing up to 96 rows and two
#' columns as follows.
#' \enumerate{
#' \item{Plate positions from A1, A2, through to H12.}
#' \item{Barcodes in the form "I" followed by 6 digits.}
#' }
#' @seealso \code{\link{readSlides}}, \code{\link{readComments}}, and
#' \code{\link{readSamples}} for reading other submission forms and
#' \code{\link{writeSampleSubmissionForm}} for usage examples.
#' @importFrom magrittr %<>%
#' @importFrom assertive.base assert_all_are_less_than_or_equal_to
#' @export
readControls <- function(file = "controls.csv")
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

#' @rdname readSlides
#' @importFrom assertive.sets assert_is_subset
#' @importFrom dplyr left_join
#' @importFrom magrittr %>%
#' @importFrom magrittr %$%
#' @importFrom stringi stri_trim_both
#' @export
readComments <- function(file = "comments.csv")
{
  comments <- fread(file, sep = ",", header = FALSE, na.strings = "", colClasses = "character")
  comments %<>%
    setNames(c("PlatePosition", "SampleNotes", "AssayNotes"))
  assert_all_are_less_than_or_equal_to(nrow(comments), 96, severity = "warning")
  not_na <- !is.na(comments$SampleNotes)
  comments$SampleNotes[not_na] <- comments$SampleNotes[not_na] %>%
    stri_trim_both %>%
    tolower %>%
    strsplit("[, ]+") %>%     # standardize separator
    vapply(paste, collapse = ", ", character(1))
  assert_is_subset(
    comments$SampleNotes, c("red", "yellow", "turbid", "red, turbid", "yellow, turbid", NA),
    severity = "warning"
  )
  data.table(PlatePosition = PLATE_POSITIONS$PlatePosition) %>%
    left_join(comments, by = "PlatePosition")
}

#' @rdname readSlides
#' @importFrom dplyr select_
#' @export
readSamples <- function(file = "samples.csv")
{
  samples <- fread(file, sep = ",", nrows = 96, header = FALSE, na.strings = "NO READ")[1:2]
  samples %<>%
    setNames(c("PlatePosition", "SampleId"))
  n_controls <- sum(is.na(samples$SampleId))
  assert_all_are_less_than_or_equal_to(n_controls, 12, severity = "warning")
  assert_all_are_not_false(
    stri_detect_regex(samples$SampleId, "^[0-9]{9}$"),
    severity = "warning"
  )
  samples
}

#' Create a SomaLogic Sample Submission Form
#'
#' Creates a data table of contents for a SomaLogic sample submission form.
#' @param slides A data frame of slide data, as imported by
#' \code{\link{readSlides}}.
#' @param controls A data frame of control data, as imported by
#' \code{\link{readControls}}.
#' @param comments A data frame of comment data, as imported by
#' \code{\link{readComments}}.
#' @param samples A data frame of sample data, as imported by
#' \code{\link{readSamples}}.
#' @param sampleMatrix A string giving the type of samples (plamsa or serum).
#' @param siteId A string giving the SomaLogic ID of your site or institution.
#' @param studyName A string giving your institution's name for the study.
#' @param studyId A string giving the SomaLogic ID of the study
#' @param runName A string naming the plate.
#' @return A data table with the 96 rows (one plate worth, with each row
#' representing a sample). It contains the following columns:
#' \describe{
#' \item{SampleNumber}{One to ninety six.}
#' \item{SlideId}{From \code{slides.csv}.}
#' \item{Subarray}{Position in slide, from 1 to 8.}
#' \item{PlatePosition}{Position in plate, from A1 to H12}
#' \item{PercentDilution}{Always 40.}
#' \item{BarCode}{From \code{controls.csv}.}
#' \item{SampleNotes}{From \code{controls.csv}.}
#' \item{AssayNotes}{From \code{controls.csv}.}
#' \item{SampleId}{From \code{samples.csv}.}
#' \item{SampleMatrix}{\code{sampleMatrix} for samples; blank for controls.}
#' \item{SiteId}{\code{siteId} for samples; blank for controls.}
#' \item{StudyName}{\code{studyName} for samples; blank for controls.}
#' \item{TimePoint}{Currently blank.}
#' \item{SampleGroup}{Currently blank.}
#' \item{SampleDescription}{Currently blank.}
#' \item{StudyId}{\code{studyId} for samples; blank for controls.}
#' \item{RunName}{\code{runName} for samples; blank for controls.}
#' }
#' @seealso \code{\link{writeSampleSubmissionForm}} for usage examples.
#' @importFrom dplyr inner_join
#' @importFrom dplyr arrange_
#' @importFrom dplyr bind_cols
#' @importFrom magrittr extract
#' @export
createSampleSubmission <- function(slides, controls, comments, samples, sampleMatrix = c("EDTA-Plasma", "Sodium Citrate Plasma", "Serum"), siteId = "WCQ", studyName = "", studyId = "", runName = "Set A")
{
  sampleMatrix <- match.arg(sampleMatrix)
  submission <- slides %>%
    inner_join(controls, by = "PlatePosition") %>%
    inner_join(comments, by = "PlatePosition") %>%
    inner_join(samples, by = "PlatePosition") %>%
    arrange_(~ SampleNumber)
  isSample <- is.na(submission$BarCode)
  details <- data.table(
    SampleMatrix = ifelse(isSample, sampleMatrix, NA_character_),
    SiteId = ifelse(isSample, siteId, NA_character_),
    StudyName = ifelse(isSample, studyName, NA_character_),
    TimePoint = NA_character_,  # TODO
    SampleGroup = NA_character_,
    SampleDescription = NA_character_,
    StudyId = ifelse(isSample, studyId, NA_character_),
    RunName = ifelse(isSample, runName, NA_character_)
  )
  submission %>%
    bind_cols(details) %>%
    extract(, c(2:4, 1, 5:17)) # Move PlatePosition column
}

#' Write a SomaLogic sample submission form
#'
#' Writes a SomaLogic sample submission form to Excel XLSX file.
#' @param submission A \code{data.frame} created by
#' \code{\link{createSampleSubmission}}.
#' @param outdir A string denoting the path to the directory where the output
#' file should be written.
#' @return The name of the XLSX file is invisibly returned, but the function is
#' mostly called for the side effect of writing this file.
#' @seealso \code{\link{readSlides}} and \code{\link{createSampleSubmission}}
#' @examples
#' \dontrun{
#' # Tests not run because data does not exist in readat yet!
#' # Import the input files
#' withr::with_dir(
#'   system.file("extdata", package = "readat"),
#'   {
#'     slides <- readSlides()
#'     controls <- readControls()
#'     comments <- readComments()
#'     samples <- readSamples()
#'   }
#' )
#'
#' # Create the sample submission form and write to Excel spreadsheet
#' submission <- createSampleSubmission(
#'   slides, controls, comments, samples,
#'   studyName = "Taheri01", studyId = "WCQ-16-002"
#' )
#' writeSampleSubmissionForm(submission)
#' }
#' @importFrom openxlsx write.xlsx
#' @importFrom pathological create_dirs
#' @export
writeSampleSubmissionForm <- function(submission, outdir = ".")
{
  create_dirs(outdir)
  filename <- paste(
      format(Sys.Date(), "%Y%m%d"),
      valuesOf(submission$StudyId),
      valuesOf(submission$RunName),
      sep = "_"
    )
  outfile <- file.path(outdir, paste0(filename, ".xlsx"))
  message("Writing to ", outfile)
  write.xlsx(submission, outfile, row.names = FALSE)
  invisible(outfile)
}

valuesOf <- function(x)
{
  toString(unique(x[!is.na(x)]))
}
