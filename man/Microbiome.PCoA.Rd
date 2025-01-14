% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Microbiome.PCoA.R
\name{Microbiome.PCoA}
\alias{Microbiome.PCoA}
\title{Perform principal coordinate analysis (PCoA).}
\usage{
Microbiome.PCoA(
  data,
  sample,
  method = "bray",
  x = "PCo1",
  y = "PCo2",
  size = 2,
  color = "group",
  alpha = 1
)
}
\arguments{
\item{data}{Feature table, with rows representing samples and columns representing feature values, such as OTUs.}

\item{sample}{Sample table, with the first column containing sample names. There are no specific requirements for the names of the subsequent columns, but the sample names must match those in the feature table.}

\item{method}{The method for calculating distances in vegan::vegdist() is defaulted to Bray-Curtis.}

\item{x}{The principal coordinate used for the X-axis in the plot is defaulted to PCo1.}

\item{y}{The principal coordinate used for the y-axis in the plot is defaulted to PCo2.}

\item{size}{The size of point.}

\item{color}{Column names in the sample data frame used for coloring the points.}

\item{alpha}{The alpha of point.}
}
\value{
A list containing a data frame and a ggplot2 object.
}
\description{
Perform principal coordinate analysis (PCoA).
}
\examples{
library(dplyr)
library(biohelpers)

data(df.pcoa.otu)
data(df.pcoa.sample)

Microbiome.PCoA(data = df.pcoa.otu, sample = df.pcoa.sample) -> pcoa.res

}
