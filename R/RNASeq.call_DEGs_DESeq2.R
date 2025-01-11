#' Identifying Differentially Expressed Genes Using DESeq2
#'
#' @param data Expression matrix must be integers, with rows representing gene names and columns representing sample names.
#' @param sample Sample information table, with rows representing sample names.
#' @param log2FoldChange The threshold for the absolute value of log2FoldChange is set to 1 by default.
#' @param padj The threshold for the adjusted p-value is set to 0.05 by default.
#'
#' @return A data frame.
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(biohelpers)
#'
#' data(df.rnaseq.gene)
#' data(df.rnaseq.sample)
#' RNASeq.call_DEGs_DESeq2(
#'   data = df.rnaseq.gene,
#'   sample = df.rnaseq.sample
#' ) -> degs
#'
RNASeq.call_DEGs_DESeq2 <- function(data, sample,  formula = ~ group, log2FoldChange = 1, padj = 0.05) {
  DESeq2::DESeqDataSetFromMatrix(
    countData = {{ data }},
    colData = {{ sample }},
    design = {{ formula }}
  ) %>%
    DESeq2::DESeq() %>%
    DESeq2::results() %>%
    as.data.frame() %>%
    dplyr::mutate(group = dplyr::case_when(
      log2FoldChange > {{ log2FoldChange }} & padj < {{ padj }} ~ "Up_regulated",
      log2FoldChange < -{{ log2FoldChange }} & padj < {{ padj }} ~ "Down_regulated",
      TRUE ~ "NS"
    )) %>%
    tibble::rownames_to_column(var = "gene") -> degs

  return(degs)
}
