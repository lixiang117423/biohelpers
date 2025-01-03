#' Title
#'
#' @param data A data frame where rows represent samples and columns represent OTUs.
#' @param physicochemical A physicochemical properties data frame where rows represent samples and columns represent physicochemical properties.
#' @param sample A sample information table where the first column is 'sample'.
#' @param method.decostand Method of vegan::decostand().
#' @param scale Scale or not in vegan::rda().
#' @param color Color for plot.
#'
#' @returns A list contains a ggplot2  object and two data frames.
#'
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(biohelpers)
#'
#' data(df.rda.otu)
#' data(df.rda.chem)
#'
#' df.rda.chem %>%
#'   rownames() %>%
#'   as.data.frame() %>%
#'   magrittr::set_names("sample") %>%
#'   dplyr::mutate(group = stringr::str_split(sample, "_") %>% sapply("[", 1)) -> df.sample
#'
#' Microbiome.RDA(df.rda.otu, df.rda.chem, df.sample) -> rda.result
#'
Microbiome.RDA <- function(data, physicochemical, sample, method.decostand = "hellinger", scale = FALSE, color = "group") {
  vegan::decostand({{ data }}, method = "hellinger") -> data.new
  vegan::rda(data.new ~ ., {{ physicochemical }}, scale = {{ scale }}) -> rda.res

  # point of physicochemical data
  vegan::scores(rda.res, display = c("sp", "wa", "lc", "bp", "cn"))$biplot %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "chem") -> df.chem.point

  # sample point
  vegan::scores(rda.res, display = "site") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample") %>%
    dplyr::left_join({{ sample }}) -> df.sample.point

  # plot
  df.sample.point %>%
    ggplot2::ggplot(ggplot2::aes(RDA1, RDA2)) +
    ggplot2::geom_point(ggplot2::aes(color = !!rlang::sym({{ color }})), size = 5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_segment(
      data = df.chem.point,
      ggplot2::aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
      arrow = ggplot2::arrow(angle = 15, length = ggplot2::unit(0.25, "cm")), color = "#222222", alpha = 0.8
    ) +
    ggrepel::geom_label_repel(data = df.chem.point, ggplot2::aes(x = RDA1, y = RDA2, label = chem)) +
    ggsci::scale_color_d3() +
    biohelpers::CommonlyUsed.plot_theme() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank()
    ) -> p

  # perm
  vegan::envfit(rda.res, {{ physicochemical }}, permutations = 999) -> envfit.res

  data.frame(r2 = envfit.res$vectors$r, pvalue = envfit.res$vectors$pvals) -> df.envfit.res

  # return
  list(rda.sample.point = df.sample.point, rda.physicochemical.point = df.chem.point, rda.envfit = df.envfit.res, rda.plot = p) %>% return()
}
