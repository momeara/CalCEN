

#' create date code of the form YYYYMMDD
#' inputs:
#'    d (optional): date code in the format of Sys.Date() for which to generate the date code, defaulting to 'today'
#' @export
date_code <- function(d = NA) {
  # reference http://www.r-cookbook.com/node/17
  if (is.na(d)) d <- Sys.Date()
  pattern <- "([[:digit:]]{4})-([[:digit:]]{2})-([[:digit:]]{2})"
  paste(
    sub(pattern, "\\1", d), sub(pattern, "\\2", d), sub(pattern, "\\3", d),
    sep = "")
}
