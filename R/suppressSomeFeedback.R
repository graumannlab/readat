#' Suppress some messages and warnings
#'
#' Suppress messages and warnings that match a regular expression.
#' @param expr Code to evaluate.
#' @param msgRegex A regular expression describing messages to suppress.
#' \code{NULL} means show all messages.
#' @param warnRegex A regular expression describing warnings to suppress.
#' \code{NULL} means show all warnings.
#' @references This uses copies of \code{evaluate_promise} and \code{with_sink}
#' from the \code{testthat} package.
#' @examples
#' \donttest{
#' suppressSomeFeedback(log(-1))
#' suppressSomeFeedback(log(-1), warnRegex = "NaN")
#' }
#' @importFrom stringi stri_detect_regex
#' @importFrom testthat evaluate_promise
#' @noRd
suppressSomeFeedback <- function(expr, msgRegex = NULL, warnRegex = NULL)
{
  # Could use pander::evals or set options(warn = 1) + capture.output. See
  # https://stat.ethz.ch/pipermail/r-devel/2015-November/072046.html
  evaluated <- testthat::evaluate_promise(expr)
  if(!is.null(msgRegex))
  {
    evaluated$messages <- evaluated$messages[
      !stri_detect_regex(evaluated$messages, msgRegex)
    ]
  }
  lapply(evaluated$messages, message)
  if(!is.null(warnRegex))
  {
    evaluated$warnings <- evaluated$warnings[
      !stri_detect_regex(evaluated$warnings, warnRegex)
    ]
  }
  lapply(evaluated$warnings, warning, call. = FALSE)
  evaluated$result
}

