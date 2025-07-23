#' Professional ggplot2 theme for biological data visualization
#'
#' @description
#' A clean and professional ggplot2 theme specifically designed for biological 
#' data visualization. This theme provides publication-ready graphics with 
#' customizable elements for scientific plots including research papers, 
#' presentations, and reports.
#'
#' @param base_size Numeric value specifying the base font size in points. 
#'   All other text elements are scaled relative to this size. Default is 10.
#' @param base_family Character string specifying the base font family. 
#'   Default is "" (system default).
#' @param strip_text_italic Logical indicating whether strip text should be 
#'   displayed in italic font when using facets. Default is FALSE.
#' @param show_minor_ticks_x Logical indicating whether to display minor ticks 
#'   on the X-axis. Default is FALSE.
#' @param show_minor_ticks_y Logical indicating whether to display minor ticks 
#'   on the Y-axis. Default is TRUE.
#' @param show_legend Logical indicating whether to display the legend. 
#'   Default is TRUE.
#' @param panel_background_color Character string specifying the panel 
#'   background color. Default is "white".
#' @param panel_border_color Character string specifying the panel border 
#'   color. Default is "black".
#' @param panel_border_size Numeric value specifying the panel border line 
#'   width. Default is 0.8.
#' @param strip_background_color Character string specifying the strip 
#'   background color for facets. Default is "#f0f0f0".
#' @param grid_lines Logical indicating whether to show major grid lines. 
#'   Default is FALSE for a clean look.
#' @param axis_line Logical indicating whether to show axis lines. Default 
#'   is FALSE (uses panel border instead).
#'
#' @return A ggplot2 theme object that can be added to any ggplot.
#'
#' @details
#' This theme provides:
#' \itemize{
#'   \item Clean, publication-ready appearance
#'   \item Proper font sizing and spacing for readability
#'   \item Professional panel borders and backgrounds
#'   \item Customizable minor ticks for precise data reading
#'   \item Flexible legend and strip text formatting
#'   \item Optimized margins and spacing for various plot types
#' }
#'
#' The theme is built on \code{ggthemes::theme_foundation()} and extends it 
#' with biology-specific design choices. It follows the principle of maximizing 
#' the data-ink ratio while maintaining visual clarity.
#'
#' @note
#' \itemize{
#'   \item This theme works best with a white or light background
#'   \item For presentations, consider increasing base_size to 12-14
#'   \item For publications, base_size of 8-10 typically works well
#'   \item Minor ticks are useful for precise quantitative plots
#' }
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(ggplot2)
#' library(biohelpers)
#'
#' # Load example data
#' data(iris)
#'
#' # Basic scatter plot with bio theme
#' iris %>%
#'   ggplot(aes(Sepal.Length, Sepal.Width, color = Species)) +
#'   geom_point(size = 2, alpha = 0.8) +
#'   labs(
#'     x = "Sepal Length (cm)",
#'     y = "Sepal Width (cm)", 
#'     title = "Iris Sepal Measurements"
#'   ) +
#'   theme_bio()
#'
#' # Customized theme for presentations
#' iris %>%
#'   ggplot(aes(Species, Sepal.Length, fill = Species)) +
#'   geom_boxplot(alpha = 0.7) +
#'   stat_summary(fun = mean, geom = "point", shape = 23, size = 3) +
#'   labs(
#'     x = "Species",
#'     y = "Sepal Length (cm)",
#'     title = "Sepal Length by Species"
#'   ) +
#'   theme_bio(
#'     base_size = 12,
#'     show_legend = FALSE,
#'     show_minor_ticks_y = TRUE
#'   )
#'
#' # Faceted plot with italic strip text
#' iris %>%
#'   tidyr::pivot_longer(
#'     cols = c(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width),
#'     names_to = "measurement",
#'     values_to = "value"
#'   ) %>%
#'   ggplot(aes(Species, value, fill = Species)) +
#'   geom_boxplot(alpha = 0.7) +
#'   facet_wrap(~measurement, scales = "free_y") +
#'   labs(
#'     x = "Species",
#'     y = "Measurement (cm)",
#'     title = "Iris Morphological Measurements"
#'   ) +
#'   theme_bio(
#'     strip_text_italic = TRUE,
#'     show_legend = FALSE
#'   )
#'
#' # Time series plot with minor ticks
#' data.frame(
#'   time = 1:24,
#'   expression = sin(1:24 * pi/6) + rnorm(24, 0, 0.1) + 2,
#'   gene = "Gene A"
#' ) %>%
#'   ggplot(aes(time, expression)) +
#'   geom_line(linewidth = 1, color = "steelblue") +
#'   geom_point(size = 2, color = "steelblue") +
#'   labs(
#'     x = "Time (hours)",
#'     y = "Expression Level",
#'     title = "Gene Expression Over Time"
#'   ) +
#'   theme_bio(
#'     show_minor_ticks_x = TRUE,
#'     show_minor_ticks_y = TRUE
#'   )
#'
theme_bio <- function(base_size = 10,
                      base_family = "",
                      strip_text_italic = FALSE,
                      show_minor_ticks_x = FALSE,
                      show_minor_ticks_y = TRUE,
                      show_legend = TRUE,
                      panel_background_color = "white",
                      panel_border_color = "black", 
                      panel_border_size = 0.8,
                      strip_background_color = "#f0f0f0",
                      grid_lines = FALSE,
                      axis_line = FALSE) {
  
  # Input validation
  if (!is.numeric(base_size) || base_size <= 0) {
    stop("'base_size' must be a positive number")
  }
  
  if (!is.character(base_family)) {
    stop("'base_family' must be a character string")
  }
  
  if (!is.character(panel_background_color)) {
    stop("'panel_background_color' must be a character string")
  }
  
  if (!is.character(panel_border_color)) {
    stop("'panel_border_color' must be a character string")
  }
  
  if (!is.numeric(panel_border_size) || panel_border_size < 0) {
    stop("'panel_border_size' must be a non-negative number")
  }
  
  # Build base theme
  bio_theme <- ggthemes::theme_foundation(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(
      # Line and text aesthetics
      line = ggplot2::element_line(colour = "black"),
      rect = ggplot2::element_rect(colour = NA, linetype = 1),
      text = ggplot2::element_text(colour = "black"),
      
      # Axis elements
      axis.line.x = if (axis_line) ggplot2::element_line() else ggplot2::element_blank(),
      axis.line.y = if (axis_line) ggplot2::element_line() else ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text.x = ggplot2::element_text(
        vjust = 0, 
        margin = ggplot2::margin(t = base_size * 0.1, unit = "pt")
      ),
      axis.text.x.top = ggplot2::element_text(
        vjust = 0, 
        margin = ggplot2::margin(b = base_size, unit = "pt")
      ),
      axis.text.y = ggplot2::element_text(
        hjust = 0.5, 
        margin = ggplot2::margin(r = base_size * 0.1, unit = "pt")
      ),
      axis.ticks = ggplot2::element_line(),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = base_size * 1, unit = "pt")
      ),
      axis.title.y = ggplot2::element_text(
        angle = 90, 
        margin = ggplot2::margin(r = base_size * 1, unit = "pt")
      ),
      axis.ticks.length = ggplot2::unit(base_size * 0.5, "points"),
      
      # Legend elements
      legend.background = ggplot2::element_rect(fill = "transparent", linetype = 0),
      legend.spacing = ggplot2::unit(base_size * 1.5, "points"),
      legend.key = ggplot2::element_rect(fill = "transparent", linetype = 0),
      legend.key.size = ggplot2::unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.0)),
      legend.text.align = NULL,
      legend.title = ggplot2::element_text(size = ggplot2::rel(1.0), hjust = 0),
      legend.title.align = NULL,
      legend.direction = NULL,
      legend.justification = "center",
      
      # Panel elements
      panel.background = ggplot2::element_rect(
        fill = panel_background_color, 
        linetype = 0
      ),
      panel.border = ggplot2::element_rect(
        color = panel_border_color, 
        linewidth = panel_border_size, 
        linetype = "solid",
        fill = NA
      ),
      panel.spacing = ggplot2::unit(0.5, "lines"),
      panel.grid.major = if (grid_lines) ggplot2::element_line(color = "grey90", linewidth = 0.5) else ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      
      # Strip elements (for facets)
      strip.background = ggplot2::element_rect(
        colour = panel_border_color, 
        fill = strip_background_color, 
        linewidth = ggplot2::rel(1.0)
      ),
      strip.text.y = ggplot2::element_text(
        color = "black", 
        size = ggplot2::rel(1.0), 
        angle = 0
      ),
      
      # Plot elements
      plot.title = ggplot2::element_text(
        size = ggplot2::rel(1.3), 
        hjust = 0.5, 
        vjust = 1.5,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        size = ggplot2::rel(1.1),
        hjust = 0.5
      ),
      plot.caption = ggplot2::element_text(
        size = ggplot2::rel(0.8),
        hjust = 1
      ),
      plot.margin = ggplot2::unit(c(6, 5, 6, 5) * 2, "points"),
      complete = FALSE
    )
  
  # Add strip text styling based on italic preference
  if (strip_text_italic) {
    bio_theme <- bio_theme + ggplot2::theme(
      strip.text.x = ggplot2::element_text(
        color = "black",
        size = ggplot2::rel(1.0),
        angle = 0,
        face = "italic"
      )
    )
  } else {
    bio_theme <- bio_theme + ggplot2::theme(
      strip.text.x = ggplot2::element_text(
        color = "black",
        size = ggplot2::rel(1.0),
        angle = 0,
        face = "plain"
      )
    )
  }
  
  # Add minor ticks if requested
  if (show_minor_ticks_x) {
    bio_theme <- bio_theme + ggplot2::theme(
      axis.minor.ticks.length.x = ggplot2::unit(base_size * 0.25, "points")
    )
  }
  
  if (show_minor_ticks_y) {
    bio_theme <- bio_theme + ggplot2::theme(
      axis.minor.ticks.length.y = ggplot2::unit(base_size * 0.25, "points")
    )
  }
  
  # Handle legend display
  if (!show_legend) {
    bio_theme <- bio_theme + ggplot2::theme(legend.position = "none")
  }
  
  return(bio_theme)
}