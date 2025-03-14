% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CommonlyUsed.PCA.R
\name{CommonlyUsed.PCA}
\alias{CommonlyUsed.PCA}
\title{Perform principal component analysis (PCA) and return the results.}
\usage{
CommonlyUsed.PCA(
  data,
  sample,
  pca.num = 5,
  plot = TRUE,
  x = "pc1",
  y = "pc2",
  color = "species",
  shape = "species"
)
}
\arguments{
\item{data}{A data frame where rows are samples and columns are indicators.}

\item{sample}{A data frame containing sample information, where one column is the sample names that match those in @param data description, without the need for row names.}

\item{pca.num}{The returned results include the number of principal components, with a default value of 10.}

\item{plot}{Whether to plot, with the default being to draw a scatter plot. The plotted graph is a ggplot2 object, which can be customized and enhanced using ggplot2.}

\item{x}{Which principal component to use as the X-axis, default is pc1.}

\item{y}{Which principal component to use as the Y-axis, default is pc2.}

\item{color}{Which grouping information to use as the color of the points, default is species.}

\item{shape}{Which grouping information to use as the shape of the points, default is species.}
}
\value{
A list containing the calculation results, scree plot, scatter plot, and COS2 plot.
}
\description{
Perform principal component analysis (PCA) and return the results.
}
\examples{
library(dplyr)
library(tibble)
library(magrittr)

data <- iris[, 1:4]
sample <- iris$Species \%>\%
  as.data.frame() \%>\%
  rownames_to_column(var = "sample") \%>\%
  set_names(c("sample", "species"))

CommonlyUsed.PCA(data, sample) -> result.pca

}
