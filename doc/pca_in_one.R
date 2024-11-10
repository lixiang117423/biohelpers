## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  #   Sepal.Length Sepal.Width Petal.Length Petal.Width
#  # 1          5.1         3.5          1.4         0.2
#  # 2          4.9         3.0          1.4         0.2
#  # 3          4.7         3.2          1.3         0.2
#  # 4          4.6         3.1          1.5         0.2
#  # 5          5.0         3.6          1.4         0.2
#  # 6          5.4         3.9          1.7         0.4

## ----eval=FALSE---------------------------------------------------------------
#  #   sample species
#  # 1      1  setosa
#  # 2      2  setosa
#  # 3      3  setosa
#  # 4      4  setosa
#  # 5      5  setosa
#  # 6      6  setosa

## ----eval=FALSE---------------------------------------------------------------
#  # devtools::install_github("lixiang117423/biohelpers")
#  
#  library(tidyverse)
#  library(biohelpers)
#  
#  data <- iris[, 1:4]
#  
#  sample <- iris$Species %>%
#    as.data.frame() %>%
#    rownames_to_column(var = "sample") %>%
#    set_names(c("sample", "species"))
#  
#  pca_in_one(data, sample) -> result.pca

