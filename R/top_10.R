#' Calculate the top ten microorganisms by abundance at different taxonomic levels.
#'
#' @param data Feature table, such as OTU abundance information, where rows represent samples and columns represent feature values.
#' @param sample Sample table, where the first column is the sample name, the second column and onwards contain grouping information, and there must be a column named 'group'.
#' @param taxo Classification information, such as the classification information of OTUs.
#' @param which Which classification information to analyze and plot, with 'phylum' as the default.
#' @param by Calculate the Top 10 by summation or averaging, with 'sum' as the default.
#' @param method Statistical analysis method, with 'wilcox' (Wilcox test) as the default.
#' @param ref The name of the control group used for statistical analysis, with 'CK' as the default.
#' @param pvalue Threshold for the p-value, with a default of 0.05.
#'
#' @return A list containing a data frame and a ggplot2 object.
#' 
#' @export
#'
#' @examples
#' 
#' library(dplyr)
#' library(biohelpers)
#' 
#' data(df.top10.otu)
#' data(df.top10.sample)
#' data(df.top10.class)
#' 
#' top_10(df.top10.otu, df.top10.sample, df.top10.class) -> top10.res
#' 
top_10 <- function(data, sample, taxo, which = "phylum", by = "sum", method = "wilcox", ref = "CK", pvalue = 0.05) {
  # merger data
  {{ data }} %>%
    tibble::rownames_to_column(var = "sample") %>%
    tidyr::pivot_longer(cols = 2:ncol(.), names_to = "OTU", values_to = "value") %>%
    dplyr::left_join({{ sample }}) %>%
    dplyr::left_join({{ taxo }}) %>%
    dplyr::select(sample, OTU, value, group, {{ which }}) -> df.all

  # Top9
  df.all %>%
    dplyr::group_by(!!rlang::sym({{ which }})) %>%
    dplyr::summarise(sum = sum(value)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-sum) %>%
    dplyr::slice(1:9) %>%
    dplyr::select(1) %>%
    magrittr::set_names("taxo") -> df.top9

  # merger raw data and top9
  df.all %>%
    dplyr::rename(taxo = {{ which }}) %>%
    dplyr::mutate(taxo = dplyr::case_when(taxo %in% df.top9$taxo ~ taxo, TRUE ~ "Other")) %>%
    dplyr::mutate(taxo = factor(taxo, level = c(df.top9$taxo, "Other"))) -> df.all.new

  # sum by sample and statistics
  df.all.new %>%
    dplyr::group_by(sample, taxo) %>%
    dplyr::summarise(sum = sum(value)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join({{ sample }}) %>%
    dplyr::group_by(taxo) %>%
    rstatix::wilcox_test(sum ~ group, ref.group = {{ ref }}) %>%
    dplyr::ungroup() -> stat.result

  # plot
  stat.result %>%
    dplyr::select(taxo, group2, p) %>%
    dplyr::rename(grooup = group2) %>%
    dplyr::mutate(signif = dplyr::case_when(p < 0.001 ~ "***", p > 0.001 & p < 0.01 ~ "**", p > 0.05 ~ "", TRUE ~ "*")) -> df.stat.new

  df.all.new %>%
    dplyr::group_by(sample, taxo) %>%
    dplyr::summarise(sum = sum(value)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(df.stat.new) %>%
    ggplot2::ggplot(ggplot2::aes(sample, sum, fill = taxo)) +
    ggplot2::geom_bar(stat = "identity", position = "fill", width = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = signif), position = ggplot2::position_fill(vjust = 0.35), size = 10, color = "black") +
    ggplot2::labs(x = NULL, y = "Relative abundance", fill = which) +
    ggplot2::scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
    ggsci::scale_fill_d3() -> plot.result

  # return
  list(stat.top_10 = stat.result, plot.top_10 = plot.result) %>% return()
}

