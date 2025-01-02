#' Calculate Correlation, Return Correlation Results.
#'
#' @param data.1 First Data Frame. Each column represents a feature value.
#' @param data.2 Second Data Frame. Each column represents a feature value.
#' @param method There are three options for calculating correlation: `pearson`, `kendall`,and `spearman`, with `pearson` being the default.
#' @param cor Correlation threshold, default is 0.5.
#' @param pvalue P-value threshold, default is 0.05.
#'
#' @return A data frame containing the calculation results.
#'
#' @export
#'
#' @examples
#'
#' library(biohelpers)
#'
#' x <- iris[, 1:4]
#' y <- iris[, 1:4]
#'
#' cor_and_plot(data.1 = x, data.2 = y) -> result.cor
#'
cor_and_plot <- function(data.1, data.2 = NULL, method = "pearson", cor = 0.6, pvalue = 0.05) {
  # If there is only one data frame, directly calculate the correlation between the variables.
  if (is.null({{ data.2 }})) {
    WGCNA::corAndPvalue(x = {{ data.1 }}, method = {{ method }}) -> cor.res

    cor.res$cor %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "from") %>%
      tidyr::pivot_longer(cols = 2:ncol(.), names_to = "to", values_to = "cor") -> df.cor

    cor.res$p %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "from") %>%
      tidyr::pivot_longer(cols = 2:ncol(.), names_to = "to", values_to = "pvalue") -> df.pvalue

    dplyr::left_join(df.cor, df.pvalue) %>%
      dplyr::filter(pvalue <= {{ pvalue }}, abs(cor) >= {{ cor }}, from != to) -> result.cor
  } else {
    WGCNA::corAndPvalue(x = {{ data.1 }}, y = {{ data.2 }}, method = {{ method }}) -> cor.res

    cor.res$cor %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "from") %>%
      tidyr::pivot_longer(cols = 2:ncol(.), names_to = "to", values_to = "cor") -> df.cor

    cor.res$p %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "from") %>%
      tidyr::pivot_longer(cols = 2:ncol(.), names_to = "to", values_to = "pvalue") -> df.pvalue

    dplyr::left_join(df.cor, df.pvalue) %>%
      dplyr::filter(pvalue <= {{ pvalue }}, abs(cor) >= {{ cor }}, from != to) -> result.cor
  }

  # return
  return(result.cor)
}
