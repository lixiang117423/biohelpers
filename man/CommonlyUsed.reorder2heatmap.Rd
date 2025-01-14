% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CommonlyUsed.reorder2heatmap.R
\name{CommonlyUsed.reorder2heatmap}
\alias{CommonlyUsed.reorder2heatmap}
\title{Title}
\usage{
CommonlyUsed.reorder2heatmap(data, col, row, value)
}
\arguments{
\item{data}{Input data in a data frame format, specifically in long data format, tidyr format.}

\item{col}{Typically, the column name where the sample names are located.}

\item{row}{Typically, the column name where the feature names are located.}

\item{value}{The column name where the numerical values needed for plotting the heatmap are located.}
}
\value{
A data frame.
}
\description{
Title
}
\examples{

library(dplyr)
library(biohelpers)

data(df.reorder2heatmap)

CommonlyUsed.reorder2heatmap(
  data = df.reorder2heatmap,
  col = "sample",
  row = "meta",
  value = "value"
) -> df.reorder

}
