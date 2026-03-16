##creating a fixed version of the generated GO panels.
##--- this can be updated at any time with V# to be changed---""
library(proteostasisTools)
library(readr)

dir.create("derived", showWarnings = FALSE)

modules_v1 <- proteostasisTools::build_proteostasis_modules_fixed(
  include_reactome = TRUE,
  go_keytype = "GOALL"
)

# Save versioned artifact
proteostasisTools::save_modules_fixed(modules_v1, "derived/modules_fixed_v1.rds")

# Save provenance for methods/supplement
write_csv(modules_v1$provenance, "derived/modules_fixed_v1_provenance.csv")

# quick sanity check
names(modules_v1$modules)
sapply(modules_v1$modules, length)
length(modules_v1$core_genes)
