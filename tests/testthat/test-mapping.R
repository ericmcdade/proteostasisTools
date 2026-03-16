test_that("fix_analyte_names converts underscore format", {
  expect_equal(fix_analyte_names("X10000_28"), "X10000.28")
  expect_equal(fix_analyte_names("ABC"), "ABC")
})

test_that("standardize_soma_annot creates expected columns", {
  annot <- data.frame(A="X1.1", B="GENE", C="P12345")
  out <- standardize_soma_annot(annot, "A", "B", "C")
  expect_true(all(c("analyte","gene_symbol","uniprot") %in% names(out)))
})

test_that("build_probe_symbol_map works without bioc", {
  annot <- data.frame(
    analyte = c("X1.1", "X2.2"),
    gene_symbol = c("ABC1|ABC2", "DEF1"),
    uniprot = c("P11111", "P22222"),
    stringsAsFactors = FALSE
  )
  out <- build_probe_symbol_map(annot, use_bioc = FALSE)
  expect_true(all(c("probe_id","SYMBOL") %in% names(out)))
  expect_true(any(out$SYMBOL %in% c("ABC1","ABC2","DEF1")))
})
