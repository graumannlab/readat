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
#' @noRd
suppressSomeFeedback <- function(expr, msgRegex = NULL, warnRegex = NULL)
{
  # Could use pander::evals or set options(warn = 1) + capture.output. See
  # https://stat.ethz.ch/pipermail/r-devel/2015-November/072046.html
  evaluated <- evaluate_promise(expr)
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


# Borrowed from testthat
evaluate_promise <- function (code, print = FALSE)
{
    warnings <- character()
    wHandler <- function(w) {
        warnings <<- c(warnings, w$message)
        invokeRestart("muffleWarning")
    }
    messages <- character()
    mHandler <- function(m) {
        messages <<- c(messages, m$message)
        invokeRestart("muffleMessage")
    }
    temp <- file()
    on.exit(close(temp))
    result <- with_sink(temp, withCallingHandlers(withVisible(code),
        warning = wHandler, message = mHandler))
    if (result$visible && print) {
        with_sink(temp, print(result$value))
    }
    output <- paste0(readLines(temp, warn = FALSE), collapse = "\n")
    list(result = result$value, output = output, warnings = warnings,
        messages = messages)
}


# Also from testthat
with_sink <- function (connection, code, ...)
{
    sink(connection, ...)
    on.exit(sink())
    code
}

