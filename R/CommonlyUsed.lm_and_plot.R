#' "Perform simple linear regression calculations and return the corresponding plot."
#'
#' @param data The input data frame contains information such as @param x, @param y, @param group, and @param color.
#' @param x Independent variable, which is also used as the X-axis in the plot.
#' @param y Dependent variable, which is also used as the Y-axis in the plot.
#' @param group The column name used to differentiate groups, which also serves as the grouping information for batch calculations of regression coefficients and P-values.
#' @param color The color of the scatter plot in the graph.
#'
#' @return A list containing a data frame and a plot, where the data frame includes regression coefficients and P-values, and the plot is a scatter plot.
#' @export
#'
#' @examples
#'
#' library(biohelpers)
#'
#' df <- iris
#'
#' CommonlyUsed.lm_and_plot(
#'   data = df,
#'   x = "Sepal.Length",
#'   y = "Sepal.Width",
#'   group = "Species",
#'   color = "Species"
#' ) -> results.lm
#'
CommonlyUsed.lm_and_plot <- function(data, x, y, group, color) {
  # Batch calculate R² and P-values.

  get_R2 <- function(fit) {
    R2 <- summary(fit)$adj.r.squared %>% round(3)
    return(R2)
  }

  get_pvalue <- function(fit) {
    pvalue <- stats::anova(fit)$`Pr(>F)`[1]
    return(pvalue)
  }

  {{ data }} %>%
    dplyr::group_by(!!rlang::sym({{ group }})) %>%
    tidyr::nest() %>%
    dplyr::mutate(model = purrr::map(data, ~ lm(as.formula(paste0({{ y }}, "~", {{ x }})), data = .x))) %>%
    dplyr::mutate(R2 = purrr::map_dbl(model, ~ get_R2(.x)), pvalue = purrr::map_dbl(model, ~ get_pvalue(.x))) %>%
    dplyr::select(!!rlang::sym({{ group }}), R2, pvalue) -> df.lm

  # plot
  {{ data }} %>%
    ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym(x), y = !!rlang::sym(y), color = !!rlang::sym(color))) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm") -> plot.lm

  # return rsults
  list(df.lm = df.lm, plot.lm = plot.lm) -> result.lm

  return(result.lm)
}
