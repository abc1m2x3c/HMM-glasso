# Run from the repository root. Packages stay in this project.
lib <- normalizePath(".Rlib", mustWork = FALSE)
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))
packages <- c("Rcpp", "glasso", "mvtnorm", "glmnet")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, lib = lib, repos = "https://cloud.r-project.org")
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
cat("Dependencies ready. Run: Rscript demo.R\n")
