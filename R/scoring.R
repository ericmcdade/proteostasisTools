#' Score pathway modules using mean Z of gene-level expression
#'
#' Z-scores each gene across subjects, then computes rowMeans for each module.
#'
#' @param expr_gene_mat Numeric matrix: rows = subjects, cols = gene SYMBOLs.
#' @param modules Named list of character vectors (gene SYMBOLs).
#' @param include_global If TRUE, includes a global score over union of module genes.
#' @return Tibble with module scores (and optionally global score).
#' @export
score_modules_mean_z <- function(expr_gene_mat, modules, include_global = TRUE) {

  if (!is.matrix(expr_gene_mat)) expr_gene_mat <- as.matrix(expr_gene_mat)

  Z <- scale(expr_gene_mat)

  # measurable genes
  modules_measured <- lapply(modules, function(v) intersect(v, colnames(expr_gene_mat)))

  out <- dplyr::tibble(.row = rownames(expr_gene_mat))

  for (nm in names(modules_measured)) {
    genes <- modules_measured[[nm]]
    out[[paste0(nm, "_z")]] <- rowMeans(Z[, genes, drop = FALSE], na.rm = TRUE)
  }

  if (include_global) {
    union_genes <- unique(unlist(modules_measured, use.names = FALSE))
    out[["proteostasis_extended_z"]] <- rowMeans(Z[, union_genes, drop = FALSE], na.rm = TRUE)
  }

  out
}
