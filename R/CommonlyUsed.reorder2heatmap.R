#' Title
#'
#' @param data Input data in a data frame format, specifically in long data format, tidyr format.
#' @param col Typically, the column name where the sample names are located.
#' @param row Typically, the column name where the feature names are located.
#' @param value The column name where the numerical values needed for plotting the heatmap are located.
#'
#' @return A data frame.
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(biohelpers)
#'
#' data(df.reorder2heatmap)
#'
#' CommonlyUsed.reorder2heatmap(
#'   data = df.reorder2heatmap,
#'   col = "sample",
#'   row = "meta",
#'   value = "value"
#' ) -> df.reorder
#'
CommonlyUsed.reorder2heatmap <- function(data, col, row, value) {
  {{ data }} %>%
    dplyr::mutate(
      col = !!rlang::sym({{ col }}),
      row = !!rlang::sym({{ row }}),
      value = !!rlang::sym({{ value }})
    ) -> data

  my_data_peak_values <- data %>%
    dplyr::group_by(row) %>%
    dplyr::slice_max(order_by = value, n = 1, with_ties = F) %>%
    dplyr::rename(peaked_at = col) %>%
    dplyr::select(-value)
  number_of_peaks <- my_data_peak_values %>%
    dplyr::group_by(peaked_at) %>%
    dplyr::count() %>%
    dplyr::arrange(-n)
  my_data_peak_values_reordered <- my_data_peak_values %>%
    dplyr::inner_join(number_of_peaks, by = "peaked_at") %>%
    dplyr::arrange(-n)
  my_data_reordered <- data %>%
    dplyr::inner_join(number_of_peaks,
      by = c(col = "peaked_at")
    ) %>%
    dplyr::mutate(col = reorder(
      col,
      -n
    )) %>%
    dplyr::select(-n) %>%
    dplyr::inner_join(my_data_peak_values_reordered,
      by = "row"
    ) %>%
    dplyr::mutate(peaked_at = reorder(
      peaked_at,
      n
    )) %>%
    dplyr::mutate(order_rows = as.numeric(peaked_at)) %>%
    dplyr::mutate(row = reorder(row, order_rows)) %>%
    dplyr::select(col, row, value)

  return(my_data_reordered)
}
