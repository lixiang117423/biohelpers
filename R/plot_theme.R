#' ggplot2 plotting theme.
#'
#' @param base.size Default font size.
#' @param strip.text.italic Whether the font is italic when displayed in groups.
#'
#' @return A ggplot2 plot theme.
#'
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(ggplot2)
#' library(biohelpers)
#'
#' data(iris)
#'
#' iris %>%
#'   ggplot(aes(Sepal.Length, Sepal.Width)) +
#'   geom_point() +
#'   plot_theme()
#'
plot_theme <- function(base.size = 10, strip.text.italic = FALSE) {
  mytheme <- ggthemes::theme_foundation(
    base_size = base.size
  ) +
    ggplot2::theme(
      line = ggplot2::element_line(colour = "black"),
      rect = ggplot2::element_rect(colour = NA, linetype = 1), text = ggplot2::element_text(colour = "black"),
      axis.line.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text.x = ggplot2::element_text(vjust = 0, margin = ggplot2::margin(t = base.size * 0.1, unit = "pt")),
      axis.text.x.top = ggplot2::element_text(vjust = 0, margin = ggplot2::margin(b = base.size, unit = "pt")),
      axis.text.y = ggplot2::element_text(hjust = 0.5, margin = ggplot2::margin(r = base.size * 0.1, unit = "pt")),
      axis.ticks = ggplot2::element_line(),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = base.size * 1, unit = "pt")),
      axis.title.y = ggplot2::element_text(angle = 90, margin = ggplot2::margin(r = base.size * 1, unit = "pt")),
      axis.ticks.length = ggplot2::unit(base.size * 0.5, "points"),
      legend.background = ggplot2::element_rect(linetype = 0),
      legend.spacing = ggplot2::unit(base.size * 1.5, "points"),
      legend.key = ggplot2::element_rect(linetype = 0),
      legend.key.size = ggplot2::unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.25)),
      legend.text.align = NULL,
      legend.title = ggplot2::element_text(size = ggplot2::rel(1), hjust = 0),
      legend.title.align = NULL,
      legend.direction = NULL,
      legend.justification = "center",
      panel.background = ggplot2::element_rect(linetype = 0),
      panel.border = ggplot2::element_rect(color = "black", size = ggplot2::rel(1), linetype = "solid"),
      panel.spacing = ggplot2::unit(0.5, "lines"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(colour = "black", fill = "#f0f0f0", size = ggplot2::rel(1.2)),
      strip.text.y = ggplot2::element_text(color = "black", size = ggplot2::rel(1.1), angle = 0),
      plot.title = ggplot2::element_text(size = ggplot2::rel(1.5), hjust = 0, vjust = 1.5),
      plot.margin = ggplot2::unit(c(6, 5, 6, 5) * 2, "points"),
      complete = FALSE
    )

  if (strip.text.italic) {
    mytheme <- mytheme + ggplot2::theme(strip.text.x = ggplot2::element_text(
      color = "black",
      size = ggplot2::rel(1.1),
      angle = 0,
      face = "italic"
    ))
  } else {
    mytheme <- mytheme + ggplot2::theme(strip.text.x = ggplot2::element_text(
      color = "black",
      size = ggplot2::rel(1.1),
      angle = 0
    ))
  }
  return(mytheme)
}
