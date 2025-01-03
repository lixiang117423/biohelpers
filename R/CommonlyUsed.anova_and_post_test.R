#' Title
#'
#' @param data A data frame.
#' @param group The column of group information.
#' @param value The column of value.
#' @param method The default method is `Tukey`.
#' @param level The default level value is 0.95.
#'
#' @returns
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(biohelpers)
#'
#' data(iris)
#'
#' CommonlyUsed.anova_and_post_test(iris, "Species", "Sepal.Length") -> anova.result
#'
CommonlyUsed.anova_and_post_test <- function(data, group, value, method = "Tukey", level = 0.95) {
  {{ data }} %>%
    dplyr::select({{ group }}, {{ value }}) %>%
    magrittr::set_names(c("group.anova", "value")) %>%
    dplyr::mutate(group.anova = factor(group.anova, levels = unique(group.anova))) -> data.new

  stats::aov(value ~ group.anova, data = data.new) -> fit

  if ({{ method }} == "Tukey") {
    multcomp::glht(fit, linfct = multcomp::mcp(group.anova = "Tukey")) %>%
      multcomp::cld(level = {{ level }}, decreasing = TRUE) -> res

    res[["mcletters"]][["Letters"]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      magrittr::set_names(c("group", "Tukey.signif")) %>%
      dplyr::mutate(anova.pvalue = broom::tidy(fit)$p.value[1]) %>%
      dplyr::select(1, 3, 2) -> res
  } else if ({{ method }} == "Duncan") {
    agricolae::duncan.test(fit, "group", alpha = 1 - {{ level }}, console = TRUE)[["groups"]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      dplyr::mutate(anova.pvalue = broom::tidy(fit)$p.value[1]) %>%
      dplyr::select(1, 4, 3) %>%
      magrittr::set_names(c("group", "anova.pvalue", "Duncan.signif")) -> res
  }

  res %>%
    dplyr::mutate(anova.signif = ifelse(anova.pvalue < 0.05 & anova.pvalue > 0.01, "*",
      ifelse(anova.pvalue < 0.01 & anova.pvalue > 0.001, "**",
        ifelse(anova.pvalue < 0.001, "***", "NS")
      )
    )) %>%
    dplyr::select(1, 2, 4, 3) -> anova.res

  return(anova.res)
}
