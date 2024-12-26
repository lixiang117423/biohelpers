#' Perform KEGG enrichment analysis of differentially expressed genes.
#'
#' @param gene The names of differentially expressed genes are a vector.
#' @param kegg.db The KEGG database used for KEGG enrichment analysis needs to
#' be in data frame format and must contain the following columns: gene, kegg.id, kegg.term.
#' The kegg.ontology column is optional.
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
#' data(df.rnaseq.kegg)
#' data(df.rnaseq.degs)
#'
#' enrich_KEGG(gene = df.rnaseq.degs$gene, kegg.db = df.rnaseq.kegg) -> kegg.res
#'
enrich_KEGG <- function(gene, kegg.db, pAdjustMethod = "BH", p.adjust = 0.05) {
  {{ kegg.db }} %>%
    dplyr::select(kegg.id, gene) -> df.term2gene

  {{ kegg.db }} %>%
    dplyr::select(kegg.id, kegg.term) -> df.term2name

  clusterProfiler::enricher(
    gene = {{ gene }},
    TERM2GENE = df.term2gene,
    TERM2NAME = df.term2name,
    qvalueCutoff = {{ p.adjust }},
    pAdjustMethod = {{ pAdjustMethod }},
    maxGSSize = length({{ gene }})
  ) -> kegg.rich

  kegg.rich@result$a <- 1
  kegg.rich@result$b <- 2

  for (i in 1:nrow(kegg.rich@result)) {
    kegg.rich@result$a[i] <- as.numeric(stringr::str_split(kegg.rich@result$GeneRatio[i], "/")[[1]][1])
    kegg.rich@result$b[i] <- as.numeric(stringr::str_split(kegg.rich@result$GeneRatio[i], "/")[[1]][2])
  }

  kegg.rich@result$GeneRatio <- kegg.rich@result$a / kegg.rich@result$b
  kegg.rich@result -> kegg.res

  return(kegg.res)
}
