# Run from the repository root. No replicate simulation study is performed.
.libPaths(c(".Rlib", .libPaths()))
source("R/graphs.R")
source("R/model.R")
source("R/evaluate.R")
legacy <- new.env(parent = globalenv())
sys.source("reference/src.R", legacy)
for (type in c("chain", "random")) for (p in c(10L, 20L, 30L)) {
  graph <- make_graphs(type, p)
  old <- with_seed(7, legacy$Graph.gen.multiK(if (type == "chain") 5 else 8, p))
  for (k in 1:3) {
    stopifnot(isTRUE(all.equal(graph$covariance[[k]], old[[paste0("Gr", k, ".C")]], tolerance = 1e-10)),
              isSymmetric(graph$adjacency[[k]]), !any(diag(graph$adjacency[[k]])),
              sum(graph$adjacency[[k]]) / 2 == p - 1 - if (k == 1) 0 else 2,
              max(abs(graph$covariance[[k]] %*% graph$precision[[k]] - diag(p))) < 1e-8)
  }
}
set.seed(77); before <- .Random.seed
invisible(make_graphs("chain"))
stopifnot(identical(before, .Random.seed))
truth <- rep(1:3, each = 4)
permuted <- c(3, 1, 2)[truth]
split_labels <- rep(1:6, each = 2)
stopifnot(adjusted_rand(permuted, truth) == 1, purity(permuted, truth) == 1,
          purity(split_labels, truth) == 1, adjusted_rand(split_labels, truth) < 1,
          adjusted_rand(rep(1, 12), truth) == 0)
for (type in c("chain", "random")) {
  path <- file.path("outputs", type, "demo.rds")
  if (!file.exists(path)) next
  saved <- readRDS(path)
  data <- saved$data
  for (method in names(saved$fits)) {
    fit <- saved$fits[[method]]
    recomputed <- posterior_hmm(data$X, fit$weights, fit$means, fit$covariance, fit$component_group)
    stopifnot(max(abs(rowSums(fit$responsibility) - 1)) < 1e-10,
      max(abs(rowSums(fit$group_probability) - 1)) < 1e-10,
      max(abs(recomputed$responsibility - fit$responsibility)) < 1e-10,
      abs(recomputed$loglik - fit$loglik) < 1e-8,
      all(is.finite(fit$loglik_history)), abs(sum(fit$weights) - 1) < 1e-10)
    for (k in seq_len(fit$K))
      stopifnot(max(abs(fit$covariance[[k]] %*% fit$precision[[k]] - diag(ncol(data$X)))) < 1e-3)
    metrics <- evaluate_fit(fit, data, method)
    stopifnot(metrics$network_purity >= 1/3, metrics$subcluster_purity >= 1/6,
              metrics$assigned_edge_sensitivity >= 0, metrics$assigned_edge_sensitivity <= 1,
              metrics$assigned_edge_specificity >= 0, metrics$assigned_edge_specificity <= 1)
  }
}
cat("All checks passed: original graph parity, graph support, RNG isolation, label invariance, and available demo fits.\n")
