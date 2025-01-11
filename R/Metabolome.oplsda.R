#' OPLS-DA.
#'
#' @param data Numerical matrix, (observations x variables; NAs are allowed), data.frame, SummarizedExperiment or ExpressionSet object.
#' @param group A factor (same length as 'x' row number) for (O)PLS-DA
#'
#' @returns A list.
#'
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(biohelpers)
#'
#' data("df.splsda.meta")
#' data("df.splsda.sample")
#'
#' df.splsda.sample %>%
#'   dplyr::filter(day == 180, treatment != "Cover") %>%
#'   dplyr::mutate(treatment = factor(treatment)) -> df.oplsda.sample
#'
#' df.splsda.meta %>%
#'   dplyr::filter(rownames(.) %in% df.oplsda.sample$sample) -> df.oplsda.meta
#'
#' Metabolome.oplsda(df.oplsda.meta, df.oplsda.sample$treatment) -> oplsda.result
#'
Metabolome.oplsda <- function(data, group, VIP = 1) {
  ropls::opls({{ data }}, group) -> oplsda.result

  # parse point
  chemhelper::get_scores(oplsda.result) -> oplsda.point

  # parse VIP
  oplsda.result@vipVn %>%
    as.data.frame() %>%
    magrittr::set_names(c("VIP")) %>%
    dplyr::filter(VIP >= {{ VIP }}) %>%
    tibble::rownames_to_column(var = "meta") %>%
    dplyr::select(meta, VIP) %>%
    dplyr::arrange(-VIP) -> oplsda.vip

  # return
  list(oplsda.result = oplsda.result, oplsda.point = oplsda.point, oplsda.vip) %>% return()
}
