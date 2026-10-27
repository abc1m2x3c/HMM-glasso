# The chain (mod 5) and random graph (mod 8) constructions in the author code.
# Keep graph randomness separate from simulation and fitting randomness.
with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)
  force(expr)
}

make_graphs <- function(type = c("chain", "random"), p = 10L, graph_seed = 123L) {
  type <- match.arg(type)
  stopifnot(length(p) == 1L, is.finite(p), p >= 8L, p == as.integer(p))
  with_seed(graph_seed, {
    if (type == "chain") {
      locations <- c(1, 1 + cumsum(runif(p - 1L, 0.5, 1)))
      covariance <- exp(-0.1 * abs(outer(locations, locations, "-")))
      base_precision <- solve(covariance)
      removed <- list(rbind(c(2, 3), c(3, 4)), rbind(c(5, 6), c(7, 8)))
    } else {
      edge_indices <- sample(choose(p, 2), p - 1L)
      adjacency <- matrix(0, p, p)
      adjacency[upper.tri(adjacency)][edge_indices] <- runif(p - 1L, 0.3, 0.8)
      adjacency <- adjacency + t(adjacency)
      base_precision <- -adjacency
      diag(base_precision) <- rowSums(adjacency) + 0.01
      partial <- -base_precision / sqrt(outer(diag(base_precision), diag(base_precision)))
      ut <- upper.tri(partial)
      strongest <- which(ut, arr.ind = TRUE)[order(abs(partial[ut]), decreasing = TRUE)[1:4], , drop = FALSE]
      removed <- list(strongest[c(1, 4), , drop = FALSE], strongest[c(2, 3), , drop = FALSE])
    }
    precision <- rep(list(base_precision), 3L)
    for (k in 2:3) for (j in seq_len(nrow(removed[[k - 1L]]))) {
      edge <- removed[[k - 1L]][j, ]
      precision[[k]][edge[1], edge[2]] <- precision[[k]][edge[2], edge[1]] <- 0
    }
    invisible(lapply(precision, chol))
    covariance <- lapply(precision, solve)
    adjacency <- lapply(precision, function(x) {
      out <- abs(x) > 1e-8
      diag(out) <- FALSE
      out
    })
    list(type = type, p = p, precision = precision, covariance = covariance,
         adjacency = adjacency, removed_edges = removed, graph_seed = graph_seed)
  })
}

simulate_hmm <- function(type = c("chain", "random"), n = 600L, p = 10L,
                         mean_shift = 1, seed = 2026L, graph_seed = 123L) {
  type <- match.arg(type)
  stopifnot(n >= 6L, n %% 6L == 0L, is.finite(mean_shift), mean_shift >= 0)
  graphs <- make_graphs(type, p, graph_seed)
  X <- with_seed(seed, do.call(rbind, lapply(0:5, function(j)
    mvtnorm::rmvnorm(n / 6L, mean = rep(j * mean_shift, p),
                    sigma = graphs$covariance[[j %/% 2L + 1L]]))))
  list(X = X, network = rep(1:3, each = n / 3L),
       subcluster = rep(1:6, each = n / 6L), graphs = graphs,
       settings = list(n = n, p = p, mean_shift = mean_shift, seed = seed,
                       graph_seed = graph_seed, K_true = 3L, L_true = 2L))
}
