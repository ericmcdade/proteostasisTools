#' Collapse probe/analyte-level expression to gene-level expression
#'
#' Given a subject x probe matrix and a probe->SYMBOL mapping,
#' collapses multiple probes per gene by rowMeans (na.rm=TRUE).
#'
#' @param expr_mat Numeric matrix/data.frame. Rows = subjects, cols = probe IDs.
#' @param probe_symbol_map Tibble/data.frame with columns `probe_id`, `SYMBOL`.
#' @return Numeric matrix: rows = subjects, cols = gene SYMBOLs.
#' @export
collapse_to_gene <- function(expr_mat, probe_symbol_map) {

  # Ensure matrix
  if (!is.matrix(expr_mat)) {
    expr_mat <- as.matrix(expr_mat)
  }

  # Keep only probes present
  m <- probe_symbol_map |>
    dplyr::filter(.data$probe_id %in% colnames(expr_mat)) |>
    dplyr::distinct(.data$probe_id, .data$SYMBOL)

  # Group probes by gene
  probe_groups <- split(m$probe_id, m$SYMBOL)

  expr_gene_mat <- vapply(names(probe_groups), function(sym) {
    probes <- probe_groups[[sym]]
    if (length(probes) == 1) {
      expr_mat[, probes]
    } else {
      rowMeans(expr_mat[, probes, drop = FALSE], na.rm = TRUE)
    }
  }, FUN.VALUE = numeric(nrow(expr_mat)))

  expr_gene_mat <- as.matrix(expr_gene_mat)
  rownames(expr_gene_mat) <- rownames(expr_mat)
  expr_gene_mat
}
