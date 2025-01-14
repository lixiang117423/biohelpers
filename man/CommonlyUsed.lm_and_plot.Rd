% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CommonlyUsed.lm_and_plot.R
\name{CommonlyUsed.lm_and_plot}
\alias{CommonlyUsed.lm_and_plot}
\title{"Perform simple linear regression calculations and return the corresponding plot."}
\usage{
CommonlyUsed.lm_and_plot(data, x, y, group, color)
}
\arguments{
\item{data}{The input data frame contains information such as @param x, @param y, @param group, and @param color.}

\item{x}{Independent variable, which is also used as the X-axis in the plot.}

\item{y}{Dependent variable, which is also used as the Y-axis in the plot.}

\item{group}{The column name used to differentiate groups, which also serves as the grouping information for batch calculations of regression coefficients and P-values.}

\item{color}{The color of the scatter plot in the graph.}
}
\value{
A list containing a data frame and a plot, where the data frame includes regression coefficients and P-values, and the plot is a scatter plot.
}
\description{
"Perform simple linear regression calculations and return the corresponding plot."
}
\examples{

library(biohelpers)

df <- iris

CommonlyUsed.lm_and_plot(
  data = df,
  x = "Sepal.Length",
  y = "Sepal.Width",
  group = "Species",
  color = "Species"
) -> results.lm

}
