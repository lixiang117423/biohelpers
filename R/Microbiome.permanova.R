#' Title
#'
#' @param formula A formula, such as df.permanova.otu ~ location * variety * group.
#' @param sample Sample information table.
#' @param method The name of any method used in vegdist to calculate pairwise distances if the left hand side of the formula was a data frame or a matrix.
#' @param permutations a list of control values for the permutations as returned by the function how, or the number of permutations required, or a permutation matrix where each row gives the permuted indices.
#'
#' @returns A data frame.
#' 
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(biohelpers)
#'
#' data(df.permanova.otu)
#' data(df.permanova.sample)
#'
#' Microbiome.permanova(df.permanova.otu ~ location * variety * group, 
#'                      df.permanova.sample) -> result.permanova
#'
Microbiome.permanova <- function(formula, sample, method = "bray", permutations = 999) {
  vegan::adonis({{ formula }},
    data = {{ sample }},
    permutations = {{ permutations }},
    method = {{ method }}
  ) -> permanova.res

  permanova.res$aov.tab %>%
    as.data.frame() -> result.permanova

  return(result.permanova)
}
