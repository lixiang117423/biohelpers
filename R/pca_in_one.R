#' Perform principal component analysis (PCA) and return the results.
#'
#' @param data A data frame where rows are samples and columns are indicators.
#' @param sample A data frame containing sample information, where one column is the sample names that match those in @param data description, without the need for row names.
#' @param pca.num The returned results include the number of principal components, with a default value of 10.
#' @param plot Whether to plot, with the default being to draw a scatter plot. The plotted graph is a ggplot2 object, which can be customized and enhanced using ggplot2.
#'
#' @return  A list containing the calculation results, scree plot, scatter plot, and COS2 plot.
#' @export pca_in_all
#' 
#' @importFrom base as.data.frame
#' @importFrom base factor
#' @importFrom base paste0
#' @importFrom base round
#' @importFrom base unique
#' @importFrom dplyr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom factoextra get_eigenvalue
#' @importFrom FactoMineR PCA
#'
#' @examples
#' x <- iris[, 1:4]
#' pca_in_one(x) -> res.tmp
#' 
utils::globalVariables(c("%>%", ".", "variance.percent", "pc"))

pca_in_one <- function(data, sample, pca.num = 5, plot = TRUE) {

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

  # prepare the return list

  return(eig.val)
}
