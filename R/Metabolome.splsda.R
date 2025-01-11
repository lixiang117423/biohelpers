#' Sparse Partial Least Squares Discriminant Analysis (sPLS-DA)
#'
#' @param data Numeric matrix of predictors with the rows as individual observations. missing values (NAs) are allowed.
#' @param group A factor or a class vector for the discrete outcome.
#' @param ncomp Positive Integer. The number of components to include in the model. Default to 3.
#'
#' @returns A list.
#'
#' @export
#'
#' @examples
#'
#' library(tidyverse)
#' library(biohelpers)
#'
#' data("df.splsda.meta")
#' data("df.splsda.sample")
#'
#' Metabolome.splsda(df.splsda.meta, as.factor(df.splsda.sample$day)) -> splsda.result
#'
Metabolome.splsda <- function(data, group, ncomp = 3) {
  mixOmics::splsda({{ data }}, {{ group }}, ncomp = {{ ncomp }}) -> splsda.result

  # parse sample point3
  splsda.result[["variates"]][["X"]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample") -> splsda.result.point

  # prop_expl_var
  splsda.result[["prop_expl_var"]][["X"]] %>%
    as.data.frame() %>%
    magrittr::set_names("prop") %>%
    dplyr::mutate(prop = (prop * 100) %>% round(2) %>% paste0("%")) %>%
    tibble::rownames_to_column(var = "comp") %>%
    dplyr::mutate(label = sprintf("X-variate %s: %s expl. var", stringr::str_remove(comp, "comp"), prop)) -> splsda.result.prop_expl_var

  # return
  list(splsda.result = splsda.result, splsda.result.point = splsda.result.point, splsda.result.prop_expl_var = splsda.result.prop_expl_var) %>% return()
}
