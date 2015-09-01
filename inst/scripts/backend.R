mGetData <- function(x, envir)
{
  x <- as.character(x)
  mget(x[!is.na(x) & nzchar(x)], envir, ifnotfound = NA_character_) %>%
    listless::list_to_data.frame(stringsAsFactors = FALSE)
}
