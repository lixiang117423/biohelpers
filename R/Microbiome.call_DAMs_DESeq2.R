#' Title
#'
#' @param data Expression matrix must be integers, with rows representing gene names and columns representing sample names.
#' @param sample Sample information table, with rows representing sample names.
#' @param log2FoldChange The threshold for the absolute value of log2FoldChange is set to 1 by default.
#' @param padj The threshold for the adjusted p-value is set to 0.05 by default.
#'
#' @return A data frame from DESeq2.
#' @export
#'
#' @examples
#'
#' library(biohelpers)
#'
#' data(df.pcoa.otu)
#' data(df.pcoa.sample)
#'
#' Microbiome.call_DAMs_DESeq2(t(df.pcoa.otu), df.pcoa.sample) -> dems.res
#'
Microbiome.call_DAMs_DESeq2 <- function(data, sample, formula = ~ group, log2FoldChange = 1, padj = 0.05) {
  DESeq2::DESeqDataSetFromMatrix(
    countData = {{ data }},
    colData = {{ sample }},
    design = {{ formula }}
  ) %>%
    DESeq2::DESeq() %>%
    DESeq2::results() %>%
    as.data.frame() %>%
    dplyr::mutate(group = dplyr::case_when(
      log2FoldChange > {{ log2FoldChange }} & padj < {{ padj }} ~ "Enriched",
      log2FoldChange < -{{ log2FoldChange }} & padj < {{ padj }} ~ "Depleted",
      TRUE ~ "NS"
    )) %>%
    tibble::rownames_to_column(var = "OTU") -> dams

  return(dams)
}