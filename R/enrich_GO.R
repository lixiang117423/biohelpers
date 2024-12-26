#' Perform GO enrichment analysis of differentially expressed genes.
#'
#' @param gene The names of differentially expressed genes are a vector.
#' @param go.db The GO database used for GO enrichment analysis needs to
#' be in data frame format and must contain the following columns: gene, go.id, go.term.
#' The go.ontology column is optional.
#' @param pAdjustMethod Methods for adjusting p-values can be selected from "holm", "hochberg",
#' "hommel", "bonferroni", "BH", "BY", "fdr", and "none", with "BH" as the default.
#' @param p.adjust Threshold for the adjusted P value.
#'
#' @return A data frame.
#' @export
#'
#' @examples
#'
#' library(biohelpers)
#'
#' data(df.rnaseq.go)
#' data(df.rnaseq.degs)
#'
#' enrich_GO(gene = df.rnaseq.degs$gene, go.db = df.rnaseq.go) -> go.res
#'
enrich_GO <- function(gene, go.db, pAdjustMethod = "BH", p.adjust = 0.05) {
  {{ go.db }} %>%
    dplyr::select(go.id, gene) -> df.term2gene

  {{ go.db }} %>%
    dplyr::select(go.id, go.term) -> df.term2name

  clusterProfiler::enricher(
    gene = {{ gene }},
    TERM2GENE = df.term2gene,
    TERM2NAME = df.term2name,
    qvalueCutoff = {{ p.adjust }},
    pAdjustMethod = {{ pAdjustMethod }},
    maxGSSize = length({{ gene }})
  ) -> go.rich

  go.rich@result$a <- 1
  go.rich@result$b <- 2

  for (i in 1:nrow(go.rich@result)) {
    go.rich@result$a[i] <- as.numeric(stringr::str_split(go.rich@result$GeneRatio[i], "/")[[1]][1])
    go.rich@result$b[i] <- as.numeric(stringr::str_split(go.rich@result$GeneRatio[i], "/")[[1]][2])
  }

  go.rich@result$GeneRatio <- go.rich@result$a / go.rich@result$b
  go.rich@result -> go.res

  return(go.res)
}
