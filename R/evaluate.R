# Label-invariant metrics. Purity permits many-to-one mappings; ARI penalizes splits.
adjusted_rand <- function(predicted, truth) {
  stopifnot(length(predicted) == length(truth), length(truth) >= 2L)
  tab <- table(predicted, truth)
  pairs <- function(x) x * (x - 1) / 2
  observed <- sum(pairs(tab))
  a <- sum(pairs(rowSums(tab))); b <- sum(pairs(colSums(tab)))
  expected <- a * b / pairs(sum(tab))
  denominator <- (a + b) / 2 - expected
  if (abs(denominator) < 1e-14) return(as.numeric(abs(observed - (a + b) / 2) < 1e-14))
  (observed - expected) / denominator
}

purity <- function(predicted, truth) {
  stopifnot(length(predicted) == length(truth), length(truth) > 0L)
  sum(apply(table(predicted, truth), 1L, max)) / length(truth)
}

precision_edges <- function(precision, edge_fraction = 0.1) {
  # Same threshold family as the original simulation analysis; diagonals excluded.
  stopifnot(is.finite(edge_fraction), edge_fraction >= 0)
  out <- abs(precision) > edge_fraction * min(abs(diag(precision)))
  diag(out) <- FALSE
  out
}

evaluate_fit <- function(fit, data, method, edge_fraction = 0.1) {
  predicted <- fit$labels
  component_predicted <- fit$component_labels
  network_tab <- table(factor(predicted, levels = seq_len(fit$K)),
                       factor(data$network, levels = 1:3))
  counts <- c(TP = 0, TN = 0, FP = 0, FN = 0)
  p <- ncol(data$X)
  ut <- upper.tri(matrix(0, p, p))
  # Each observation receives its assigned group's estimated network.
  # Pool observation/edge pairs: no label matching or overwritten duplicate groups.
  for (j in seq_len(fit$K)) for (k in 1:3) {
    estimate <- precision_edges(fit$precision[[j]], edge_fraction)[ut]
    truth <- data$graphs$adjacency[[k]][ut]
    counts <- counts + network_tab[j, k] * c(
      TP = sum(estimate & truth), TN = sum(!estimate & !truth),
      FP = sum(estimate & !truth), FN = sum(!estimate & truth))
  }
  data.frame(method = method, fitted_network_groups = fit$K,
    fitted_mean_components = sum(fit$L),
    network_ARI = adjusted_rand(predicted, data$network),
    network_purity = purity(predicted, data$network),
    subcluster_ARI = adjusted_rand(component_predicted, data$subcluster),
    subcluster_purity = purity(component_predicted, data$subcluster),
    assigned_edge_sensitivity = counts['TP'] / (counts['TP'] + counts['FN']),
    assigned_edge_specificity = counts['TN'] / (counts['TN'] + counts['FP']),
    iterations = fit$iterations, converged = fit$converged,
    row.names = NULL)
}
