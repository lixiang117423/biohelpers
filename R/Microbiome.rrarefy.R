#' Rarefaction Species Richness
#'
#' @param data Community data, a matrix-like object or a vector.
#' @param value Subsample size for rarefying community, either a single value or a vector. The default is the minimum value in the sample.
#'
#' @returns A data frame.
#'
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(biohelpers)
#'
#' data(df.permanova.otu)
#'
#' df.permanova.otu %>%
#'   t() %>%
#'   as.data.frame() %>%
#'   Microbiome.rrarefy() -> df.otu.rrarefy
#'
Microbiome.rrarefy <- function(data, value = NULL) {
  if (is.null({{ value }})) {
    {{ data }} %>%
      t() %>%
      vegan::rrarefy(min(colSums({{ data }}))) %>%
      t() %>%
      as.data.frame() -> df.otu.rrarefy
  } else {
    {{ data }} %>%
      t() %>%
      vegan::rrarefy(sample = {{ value }}) %>%
      t() %>%
      as.data.frame() -> df.otu.rrarefy
  }
  return(df.otu.rrarefy)
}
