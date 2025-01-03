#' Convert data frame to list for venn plot.
#'
#' @param data A data frame.
#' @param group The column name for group.
#' @param value The column name for value.
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
#' data.frame(group = "A", value = 1:10) -> df.1
#' data.frame(group = "B", value = 3:12) -> df.2
#' data.frame(group = "C", value = 7:20) -> df.3
#'
#' dplyr::bind_rows(df.1, df.2, df.3) %>%
#'   CommonlyUsed.df_to_list(group = "group", value = "value") -> res.list
#'
CommonlyUsed.df_to_list <- function(data, group, value) {
  {{ data }} %>%
    dplyr::select({{ group }}, {{ value }}) %>%
    dplyr::mutate(id = rownames(.)) %>%
    dplyr::select(3, 1, 2) %>%
    dplyr::rename(value = {{ value }}) %>%
    tidyr::pivot_wider(id_cols = 1, names_from = {{ group }}) %>%
    dplyr::select(-id) %>%
    as.list() %>%
    lapply(function(x) x[!is.na(x)]) -> res.list

  return(res.list)
}
