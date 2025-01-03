#' Create a volcano plot.
#'
#' @param data The data frame used to create the scatter plot is typically the output of DESeq2, specifically the output from call_DEGs_DESeq2.
#' @param x X axis.
#' @param y Y axis.
#' @param size Column name used to map the size of the points, default is baseMean.
#' @param color Column name used to map the color of the points.
#' @param xintercept.abs Intercept of X axes. It is usually log2FoldChange, which includes both positive and negative values.
#' @param yintercept Intercept of Y axes. It is usually -log10(padj).
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#'
#' library(biohelpers)
#'
#' data(df.rnaseq.plot_volcano)
#'
#' RNASeq.plot_volcano(data = df.rnaseq.plot_volcano) -> p
#'
RNASeq.plot_volcano <- function(data, x = "log2FoldChange", y = "padj", size = "baseMean", color = "group", xintercept.abs = 1, yintercept = 0.05) {
  {{ data }} %>%
    ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym({{ x }}), y = -log10(!!rlang::sym({{ y }})), color = !!rlang::sym({{ color }}))) +
    ggplot2::geom_point(ggplot2::aes(size = !!rlang::sym({{ size }})), alpha = 0.8) +
    ggplot2::geom_hline(yintercept = -log10({{ yintercept }}), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c({{ xintercept.abs }}, -{{ xintercept.abs }}), linetype = "dashed") +
    ggsci::scale_color_d3() +
    biohelpers::plot_theme() -> p

  return(p)
}
