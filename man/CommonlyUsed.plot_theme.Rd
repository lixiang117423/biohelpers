% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CommonlyUsed.plot_theme.R
\name{CommonlyUsed.plot_theme}
\alias{CommonlyUsed.plot_theme}
\title{ggplot2 plotting theme.}
\usage{
CommonlyUsed.plot_theme(
  base.size = 10,
  strip.text.italic = FALSE,
  show.minor.tricks.x = FALSE,
  show.minor.tricks.y = TRUE,
  show.legend = TRUE
)
}
\arguments{
\item{base.size}{Default font size.}

\item{strip.text.italic}{Whether the font is italic when displayed in groups.}

\item{show.minor.tricks.x}{Whether to display the minor ticks on the X-axis.}

\item{show.minor.tricks.y}{Whether to display the minor ticks on the y-axis.}

\item{show.legend}{Whether to display the legend.}
}
\value{
A ggplot2 plot theme.
}
\description{
ggplot2 plotting theme.
}
\examples{

library(dplyr)
library(ggplot2)
library(biohelpers)

data(iris)

iris \%>\%
  ggplot(aes(Sepal.Length, Sepal.Width)) +
  geom_point() +
  CommonlyUsed.plot_theme()

}
