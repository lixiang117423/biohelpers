#' Perform principal component analysis (PCA) and return the results.
#'
#' @param data A data frame where rows are samples and columns are indicators.
#' @param sample A data frame containing sample information, where one column is the sample names that match those in @param data description, without the need for row names.
#' @param pca.num The returned results include the number of principal components, with a default value of 10.
#' @param plot Whether to plot, with the default being to draw a scatter plot. The plotted graph is a ggplot2 object, which can be customized and enhanced using ggplot2.
#' @param x Which principal component to use as the X-axis, default is pc1.
#' @param y Which principal component to use as the Y-axis, default is pc2.
#' @param color Which grouping information to use as the color of the points, default is species.
#' @param shape Which grouping information to use as the shape of the points, default is species.
#'
#' @return A list containing the calculation results, scree plot, scatter plot, and COS2 plot.
#' @export
#'
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(magrittr)
#'
#' data <- iris[, 1:4]
#' sample <- iris$Species %>%
#'   as.data.frame() %>%
#'   rownames_to_column(var = "sample") %>%
#'   set_names(c("sample", "species"))
#'
#' pca_in_one(data, sample) -> result.pca
#'
pca_in_one <- function(data, sample, pca.num = 5, plot = TRUE, x = "pc1", y = "pc2", color = "species", shape = "species") {
  # run PCA
  FactoMineR::PCA({{ data }}, ncp = {{ pca.num }}, graph = FALSE) -> pca.res

  pca.res %>%
    factoextra::get_eigenvalue() %>%
    base::as.data.frame() %>%
    dplyr::mutate(pc = base::paste0("PC", 1:nrow(.))) %>%
    dplyr::arrange(-variance.percent) %>%
    dplyr::mutate(pc = base::factor(pc, levels = base::unique(pc))) %>%
    dplyr::mutate(variance.percent = base::round(variance.percent, 2)) %>%
    dplyr::select(pc, variance.percent) -> eig.val

  # point data
  pca.res[["ind"]][["coord"]] %>%
    as.data.frame() %>%
    magrittr::set_names(c(paste0("pc", 1:ncol(.)))) %>%
    tibble::rownames_to_column(var = "sample") -> point.data

  # point plot
  pca.res[["ind"]][["coord"]] %>%
    as.data.frame() %>%
    magrittr::set_names(c(paste0("pc", 1:ncol(.)))) %>%
    tibble::rownames_to_column(var = "sample") %>%
    dplyr::left_join({{ sample }}, by = "sample") %>%
    ggplot2::ggplot(ggplot2::aes(!!rlang::sym({{ x }}), !!rlang::sym({{ y }}), color = !!rlang::sym({{ color }}), shape = !!rlang::sym({{ shape }}))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      x = paste0("PC", stringr::str_replace({{ x }}, "pc", ""), "(", eig.val$variance.percent[stringr::str_replace({{ x }}, "pc", "") %>% as.numeric()], "%)"),
      y = paste0("PC", stringr::str_replace({{ y }}, "pc", ""), "(", eig.val$variance.percent[stringr::str_replace({{ y }}, "pc", "") %>% as.numeric()], "%)")
    ) -> plot.pca

  # prepare the return list
  return.list <- list(result.pca = pca.res, plot.pca = plot.pca, point.data = point.data, eigenvalue.pca = eig.val)

  return(return.list)
}












