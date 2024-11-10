#' Identification of Outliers.
#'
#' @param x Numeric Format Vector.
#'
#' @return Vector. If it is an outlier, return 'yes'; otherwise, return 'no'.
#'
#' @export
#'
#' @examples
#' 
#' library(dplyr)
#' library(biohelpers)
#'
#' data(iris)
#'
#' iris %>%
#'   dplyr::group_by(Species) %>%
#'   dplyr::mutate(is.out = find_outliner(Sepal.Length)) %>%
#'   dplyr::ungroup() -> result.find_outlinr
#'
find_outliner <- function(x) {
  return(
    ifelse(
      x < stats::quantile(x, .25) - 1.5 * stats::IQR(x) | x > stats::quantile(x, .75) + 1.5 * stats::IQR(x),
      "yes",
      "no"
    )
  )
}
