#!/usr/bin/env Rscript
# One dataset per requested graph type. No Monte Carlo replication loop.
cli <- commandArgs(trailingOnly = TRUE)
defaults <- list(graph = "both", n = "600", p = "10", mean_shift = "1",
                 seed = "2026", graph_seed = "123", K = "auto", rho = "auto",
                 max_iter = "1000", bic = "legacy", out = "outputs")
usage <- function() {
  cat("HMM-glasso: one dataset per graph type\n\n",
      "Rscript demo.R [--graph=both|chain|random] [--n=600] [--p=10]\n",
      "  [--mean-shift=1] [--seed=2026] [--graph-seed=123]\n",
      "  [--K=auto|3] [--rho=auto|0.1] [--max-iter=1000]\n",
      "  [--bic=legacy|full] [--out=outputs]\n", sep = "")
}
if ("--help" %in% cli) { usage(); quit(status = 0) }
for (arg in cli) {
  if (!grepl("^--[A-Za-z-]+=.+$", arg)) stop("Use --name=value; see --help.")
  pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
  name <- gsub("-", "_", pair[1])
  if (!name %in% names(defaults)) stop("Unknown option: ", pair[1])
  defaults[[name]] <- paste(pair[-1], collapse = "=")
}
cfg <- defaults
for (name in c("n", "p", "seed", "graph_seed", "max_iter")) {
  value <- suppressWarnings(as.numeric(cfg[[name]]))
  if (length(value) != 1L || !is.finite(value) || value != as.integer(value) || value < 1L)
    stop(name, " must be a positive integer.")
  cfg[[name]] <- as.integer(value)
}
cfg$mean_shift <- as.numeric(cfg$mean_shift)
stopifnot(cfg$graph %in% c("both", "chain", "random"), cfg$n %% 6L == 0,
          cfg$n >= 180L, cfg$p >= 8L, cfg$mean_shift >= 0, is.finite(cfg$mean_shift),
          cfg$bic %in% c("legacy", "full"))
K <- if (cfg$K == "auto") NULL else as.numeric(cfg$K)
if (!is.null(K)) stopifnot(is.finite(K), K >= 1, K == as.integer(K), cfg$n >= 6 * K)
rho <- if (cfg$rho == "auto") NULL else as.numeric(cfg$rho)
if (!is.null(rho)) stopifnot(is.finite(rho), rho > 0)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
for (package in c("glasso", "mvtnorm", "glmnet"))
  if (!requireNamespace(package, quietly = TRUE))
    stop("Missing ", package, ". Run Rscript scripts/install_dependencies.R from the repository root.")
for (file in c("graphs.R", "model.R", "evaluate.R", "plots.R")) source(file.path(root, "R", file))
dir.create(cfg$out, showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(cfg$out, "sessionInfo.txt"))
graphs <- if (cfg$graph == "both") c("chain", "random") else cfg$graph
all_metrics <- list()
for (graph in graphs) {
  # Graph-specific seed makes --graph=random identical to its run in --graph=both.
  data_seed <- cfg$seed + if (graph == "random") 100L else 0L
  data <- simulate_hmm(graph, cfg$n, cfg$p, cfg$mean_shift, data_seed, cfg$graph_seed)
  message("[", graph, "] Simulated one dataset: n=", cfg$n, ", p=", cfg$p, ", true K=3, true L=2.")
  tuning <- if (is.null(rho)) choose_rho(data$X, data_seed + 1L) else list(rho = rho, rule = "user supplied")
  message("[", graph, "] rho=", signif(tuning$rho, 5), "; fitting MM and MM-glasso (6 components each).")
  mm <- fit_mm(data$X, 6L, rho = 0, max_iter = cfg$max_iter, seed = data_seed + 2L)
  mm_glasso <- fit_mm(data$X, 6L, rho = tuning$rho, max_iter = cfg$max_iter, seed = data_seed + 3L)
  message("[", graph, "] Fitting HMM-glasso; K=", cfg$K, "; comparing L=2 and L=3.")
  selected <- select_hmm(data$X, tuning$rho, K = K, bic = cfg$bic,
                         seed = data_seed + 4L, max_iter = cfg$max_iter)
  fits <- list(MM = mm, `MM-glasso` = mm_glasso, `HMM-glasso` = selected$fit)
  metrics <- do.call(rbind, lapply(names(fits), function(name) evaluate_fit(fits[[name]], data, name)))
  metrics <- cbind(graph = graph, metrics)
  all_metrics[[graph]] <- metrics
  folder <- file.path(cfg$out, graph)
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  write.csv(metrics, file.path(folder, "metrics.csv"), row.names = FALSE)
  write.csv(selected$candidates, file.path(folder, "model_selection.csv"), row.names = FALSE)
  assignments <- data.frame(observation = seq_len(cfg$n), true_network = data$network,
    true_subcluster = data$subcluster, MM = mm$labels, MM_glasso = mm_glasso$labels,
    HMM_network = selected$fit$labels, HMM_subcluster = selected$fit$component_labels)
  write.csv(assignments, file.path(folder, "assignments.csv"), row.names = FALSE)
  write.csv(data.frame(observation = seq_len(cfg$n), data$X), file.path(folder, "simulated_data.csv"), row.names = FALSE)
  saveRDS(list(data = data, fits = fits, selection = selected, tuning = tuning, config = cfg), file.path(folder, "demo.rds"))
  save_network_plot(data, selected$fit, file.path(folder, "networks.png"))
  save_diagnostics(data, selected, file.path(folder, "diagnostics.png"))
  for (name in names(fits)) if (!fits[[name]]$converged)
    message("[", graph, "] ", name, " reached max_iter; its convergence flag is recorded as FALSE.")
  message("[", graph, "] Selected K=", selected$fit$K, ", L=", paste(unique(selected$fit$L), collapse = ","), ".")
  print(metrics[, c("method", "fitted_network_groups", "network_ARI", "network_purity", "converged")], row.names = FALSE)
}
write.csv(do.call(rbind, all_metrics), file.path(cfg$out, "metrics.csv"), row.names = FALSE)
message("Finished. Results saved in ", normalizePath(cfg$out), ". Single-run results are illustrative, not comparative performance claims.")
