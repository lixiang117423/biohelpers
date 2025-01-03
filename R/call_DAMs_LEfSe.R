#' Identifying differential abundance microbiomes using LEfSe.
#'
#' @param data Feature table, with rows representing feature values and columns representing sample names.
#' @param sample Sample table, with rows representing sample name and columns representing group information.
#' @param groupCol The column name where the grouping information is located, default is 'group'.
#' @param kruskal.threshold The threshold for the Kruskal test, default is 0.05.
#' @param wilcox.threshold The threshold for the Wilcox test, default is 0.05.
#' @param lda.threshold The threshold for the LDA test, default is 1.
#'
#' @returns A data frame.
#'
#' @export
#'
#' @examples
#'
#' library(biohelpers)
#' library(dplyr)
#'
#' data(df.call_DAMs_LEfSe.otu)
#' data(df.call_DAMs_LEfSe.sample)
#'
#' call_DAMs_LEfSe(df.call_DAMs_LEfSe.otu, df.call_DAMs_LEfSe.sample) -> lefse.result
#'
call_DAMs_LEfSe <- function(data, sample, groupCol = "group", kruskal.threshold = 0.05, wilcox.threshold = 0.05, lda.threshold = 1) {
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = {{ data }}), colData = {{ sample }}) %>%
    lefser::lefser(groupCol = {{ groupCol }}, kruskal.threshold = {{ kruskal.threshold }}, wilcox.threshold = {{ wilcox.threshold }}, lda.threshold = {{ lda.threshold }}) -> lefse.result

  return(lefse.result)
}
