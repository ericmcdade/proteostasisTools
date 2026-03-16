#' Fix SomaScan analyte IDs
#'
#' Converts analyte IDs of the form `X10000_28` to `X10000.28`.
#' Leaves other strings unchanged.
#'
#' @param x Character vector of analyte IDs.
#' @return Character vector with `_` converted to `.` for analyte-like tokens.
#' @export
fix_analyte_names <- function(x) {
  x <- as.character(x)
  sub("^X(\\d+)_(\\d+)$", "X\\1.\\2", x)
}

#' Standardize a SomaScan annotation table to a common schema
#'
#' Many SomaScan annotation files use different column names.
#' This helper returns a standardized tibble with columns:
#' `analyte`, `gene_symbol`, `uniprot`.
#'
#' @param annot data.frame/tibble with annotation columns.
#' @param analyte_col Column name for analyte IDs.
#' @param symbol_col Column name for gene symbols.
#' @param uniprot_col Column name for UniProt IDs.
#' @return Tibble with `analyte`, `gene_symbol`, `uniprot`.
#' @export
standardize_soma_annot <- function(annot, analyte_col, symbol_col, uniprot_col) {
  dplyr::tibble(
    analyte     = annot[[analyte_col]],
    gene_symbol = annot[[symbol_col]],
    uniprot     = annot[[uniprot_col]]
  )
}

#' Build a probe/analyte -> gene SYMBOL map from SomaScan annotation
#'
#' This function parses gene symbols and UniProt IDs (optionally via Bioconductor)
#' to generate a mapping from analyte/probe IDs to HGNC-like gene symbols.
#'
#' Use `use_bioc = TRUE` if you have Bioconductor installed, which improves robustness
#' by validating symbols and mapping UniProt IDs to symbols.
#'
#' @param soma_annot Annotation table. Can be standardized (recommended) or raw.
#' @param probe_id_col Column name for analyte/probe IDs (default "analyte").
#' @param gene_symbol_col Column name for gene symbols (default "gene_symbol").
#' @param uniprot_col Column name for UniProt IDs (default "uniprot").
#' @param use_bioc Logical. If TRUE, uses AnnotationDbi + org.Hs.eg.db (Suggests).
#' @return Tibble with columns `probe_id`, `SYMBOL` (one row per mapping).
#' @export
build_probe_symbol_map <- function(soma_annot,
                                   probe_id_col = "analyte",
                                   gene_symbol_col = "gene_symbol",
                                   uniprot_col = "uniprot",
                                   use_bioc = TRUE) {

  # ---- internal helpers ----
  clean_token <- function(x) {
    x |>
      as.character() |>
      stringr::str_replace_all("\u00A0", " ") |>
      stringr::str_squish() |>
      toupper()
  }

  split_multi <- function(x) {
    if (is.na(x) || x == "") return(character(0))
    clean_token(x) |>
      stringr::str_split("[;|,/]+") |>
      (\(z) z[[1]])() |>
      stringr::str_squish() |>
      purrr::discard(~ .x == "" | .x == "NA")
  }

  strip_uniprot_isoform <- function(x) stringr::str_replace(x, "-\\d+$", "")

  # ---- Extract all candidate keys ----
  symbol_tokens <- soma_annot[[gene_symbol_col]] |>
    purrr::map(split_multi) |>
    unlist() |>
    unique()

  uniprot_tokens <- soma_annot[[uniprot_col]] |>
    purrr::map(split_multi) |>
    unlist() |>
    strip_uniprot_isoform() |>
    unique()

  valid_symbols <- symbol_tokens
  valid_uniprot <- uniprot_tokens
  map_uniprot <- NULL

  # ---- Optional Bioconductor augmentation ----
  if (use_bioc) {
    if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("use_bioc=TRUE requires AnnotationDbi and org.Hs.eg.db (install via Bioconductor).")
    }

    valid_symbols <- intersect(
      symbol_tokens,
      AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
    )

    valid_uniprot <- intersect(
      uniprot_tokens,
      AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "UNIPROT")
    )

    map_uniprot <- suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = valid_uniprot,
        keytype = "UNIPROT",
        columns = "SYMBOL"
      )
    )
  }

  # ---- Build mapping ----
  soma_annot |>
    dplyr::rowwise() |>
    dplyr::mutate(
      probe_id = .data[[probe_id_col]],
      symbol_list = list(split_multi(.data[[gene_symbol_col]])),
      uniprot_list = list(split_multi(.data[[uniprot_col]]) |> strip_uniprot_isoform())
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      SYMBOL = purrr::map2(symbol_list, uniprot_list, function(sym, uni) {

        sym_hits <- sym[sym %in% valid_symbols]

        uni_hits <- character(0)
        if (use_bioc && !is.null(map_uniprot)) {
          uni_hits <- map_uniprot |>
            dplyr::filter(.data$UNIPROT %in% uni) |>
            dplyr::pull(.data$SYMBOL)
        }

        unique(c(sym_hits, uni_hits))
      })
    ) |>
    dplyr::select("probe_id", "SYMBOL") |>
    tidyr::unnest("SYMBOL") |>
    dplyr::filter(!is.na(.data$SYMBOL), .data$SYMBOL != "") |>
    dplyr::distinct()
}
