#'  Plot the manhattan from a GWAS results.
#' @param df A data frame.
#' @param p.threshold The threshold of P value, the defalut value is 0.05 / number of SNPs.
#'
#' @returns A ggplot2 object.
#'
#' @export
#'
#' @examples
#'
#' library(dplyr)
#'
#' # 设置随机数种子
#' set.seed(20250312)
#'
#' # 构造GWAS数据
#' normentR::simulateGWAS(nSNPs = 1e6, nSigCols = 3) %>%
#'   janitor::clean_names() %>%
#'   dplyr::rename(posi = bp) %>%
#'   dplyr::mutate(label = case_when(p < 0.05 / nrow(.) ~ snp, TRUE ~ NA)) -> df.gwas
#'
#' PopGen.manhattan_plot(df.gwas)
#'
PopGen.manhattan_plot <- function(df, p.threshold = NA) {
  {{ df }} %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(max_posi = max(posi)) %>%
    dplyr::mutate(posi_add = dplyr::lag(cumsum(max_posi), default = 0)) %>%
    dplyr::ungroup() %>%
    dplyr::select(chr, posi_add) -> df.chr.posi

  {{ df }} %>%
    dplyr::left_join(df.chr.posi) %>%
    dplyr::mutate(posi_cum = posi + posi_add) -> df.gwas.new

  df.gwas.new %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(chr.center = mean(posi_cum)) %>%
    dplyr::ungroup() -> df.chr.center

  df.gwas.new %>%
    dplyr::filter(p == min(p)) %>%
    dplyr::mutate(ymax = log10(p) %>% floor() %>% abs() + 2) %>%
    dplyr::pull(ymax) -> ymax

  if (is.na({{ p.threshold }})) {
    0.05 / nrow({{ df }}) -> p.threshold.plot
  } else {
    {{ p.threshold }} -> p.threshold.plot
  }

  df.gwas.new %>%
    ggplot2::ggplot(aes(x = posi_cum, y = -log10(p), color = forcats::as_factor(chr), size = -log10(p))) +
    ggplot2::geom_hline(yintercept = -log10(p.threshold.plot), color = "grey40", linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.5) +
    ggrepel::geom_label_repel(ggplot2::aes(label = label)) +
    ggplot2::scale_x_continuous(label = df.chr.center$chr, breaks = df.chr.center$chr.center) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, ymax), breaks = seq(0, ymax, 2)) +
    ggplot2::scale_color_manual(values = rep(c("#1F77B4FF", "#FF7F0EFF"), unique(length(df.chr.center$chr)))) +
    ggplot2::scale_size_continuous(range = c(0.5, 3)) +
    ggplot2::labs(x = NULL, y = "-log<sub>10</sub>(p)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.title.y = ggtext::element_markdown(),
      axis.text.x = ggplot2::element_text(angle = 0, size = 8, vjust = 0.5)
    ) -> p

  return(p)
}
