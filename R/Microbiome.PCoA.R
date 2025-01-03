#' Perform principal coordinate analysis (PCoA).
#'
#' @param data Feature table, with rows representing samples and columns representing feature values, such as OTUs.
#' @param sample Sample table, with the first column containing sample names. There are no specific requirements for the names of the subsequent columns, but the sample names must match those in the feature table.
#' @param method The method for calculating distances in vegan::vegdist() is defaulted to Bray-Curtis.
#' @param x The principal coordinate used for the X-axis in the plot is defaulted to PCo1.
#' @param y The principal coordinate used for the y-axis in the plot is defaulted to PCo2.
#' @param size The size of point.
#' @param color Column names in the sample data frame used for coloring the points.
#' @param alpha The alpha of point.
#'
#' @return A list containing a data frame and a ggplot2 object.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(biohelpers)
#'
#' data(df.pcoa.otu)
#' data(df.pcoa.sample)
#'
#' Microbiome.PCoA(data = df.pcoa.otu, sample = df.pcoa.sample) -> pcoa.res
#'
Microbiome.PCoA <- function(data, sample, method = "bray", x = "PCo1", y = "PCo2", size = 2, color = "group", alpha = 1) {
  {{ data }} %>%
    vegan::vegdist(method = {{ method }}) %>%
    ape::pcoa() -> pcoa.res

  pcoa.res$values %>%
    as.data.frame() %>%
    dplyr::select(2) %>%
    dplyr::mutate(PCoA = paste0("PCo", rownames(.))) %>%
    dplyr::select(PCoA, 1) %>%
    dplyr::mutate(relative.eig = (Relative_eig * 100) %>% round(2) %>% paste0("%")) %>%
    dplyr::select(PCoA, relative.eig) %>%
    dplyr::mutate(label = paste0(PCoA, " (", relative.eig, ")")) -> pcoa.weight

  # sample point
  pcoa.res$vectors %>%
    as.data.frame() %>%
    dplyr::select(1:5) %>%
    magrittr::set_names(paste0("PCo", 1:ncol(.))) %>%
    tibble::rownames_to_column(var = "sample") %>%
    dplyr::left_join({{ sample }}) -> pcoa.point

  # PCoA plot
  pcoa.weight %>%
    dplyr::filter(PCoA == {{ x }}) %>%
    dplyr::select(label) -> df.x.label

  pcoa.weight %>%
    dplyr::filter(PCoA == {{ y }}) %>%
    dplyr::select(label) -> df.y.label

  pcoa.point %>%
    ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym({{ x }}), y = !!rlang::sym({{ y }}), color = !!rlang::sym({{ color }}))) +
    ggplot2::geom_hline(yintercept = 0, color = "#222222", linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, color = "#222222", linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_point(size = size, alpha = alpha) +
    ggplot2::labs(x = df.x.label$label, y = df.x.label$label) +
    ggsci::scale_color_d3() -> p

  # return list
  list(pcoa.point = pcoa.point, pcoa.weight = pcoa.weight, p.pcoa = p) %>% return()
}
