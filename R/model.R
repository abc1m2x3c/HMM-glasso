# Portable EM-style updates adapted from EM.glasso and hierarchical_gmm3.
# Each network group shares one covariance across its mean components.
# rho is the penalty passed to glasso on the normalized weighted scatter.
posterior_hmm <- function(X, weights, means, covariance, component_group) {
  log_density <- vapply(seq_along(weights), function(j)
    log(weights[j]) + mvtnorm::dmvnorm(X, means[j, ],
      covariance[[component_group[j]]], log = TRUE), numeric(nrow(X)))
  log_density <- matrix(log_density, nrow(X))
  row_max <- apply(log_density, 1L, max)
  if (any(!is.finite(row_max))) stop("Non-finite mixture density.")
  shifted <- exp(log_density - row_max)
  normalizer <- rowSums(shifted)
  list(responsibility = shifted / normalizer,
       loglik = sum(row_max + log(normalizer)))
}

fit_hmm <- function(X, K, L = 2L, rho = 0.1, max_iter = 1000L,
                    tol = 1e-6, seed = 1L, ridge = 1e-8) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  stopifnot(is.numeric(X), all(is.finite(X)), p >= 2L,
            length(K) == 1L, is.finite(K), K >= 1L, K == as.integer(K),
            length(rho) == 1L, is.finite(rho), rho >= 0,
            length(max_iter) == 1L, max_iter >= 1L, max_iter == as.integer(max_iter),
            is.finite(tol), tol > 0, is.finite(ridge), ridge > 0)
  if (length(L) == 1L) L <- rep(L, K)
  stopifnot(length(L) == K, all(is.finite(L)), all(L >= 1L),
            all(L == as.integer(L)), n >= 2L * sum(L))
  component_group <- rep(seq_len(K), L)
  J <- sum(L)
  initial <- with_seed(seed, {
    group <- sample(rep(seq_len(K), length.out = n))
    means <- matrix(0, J, p)
    covariance <- vector("list", K)
    for (k in seq_len(K)) {
      rows <- which(group == k)
      covariance[[k]] <- stats::cov(X[rows, , drop = FALSE]) + diag(ridge, p)
      subgroup <- rep(seq_len(L[k]), length.out = length(rows))
      indices <- which(component_group == k)
      for (l in seq_len(L[k]))
        means[indices[l], ] <- colMeans(X[rows[subgroup == l], , drop = FALSE])
    }
    list(means = means, covariance = covariance)
  })
  means <- initial$means; covariance <- initial$covariance
  weights <- (1 / K) / L[component_group]
  state <- posterior_hmm(X, weights, means, covariance, component_group)
  history <- state$loglik
  converged <- FALSE
  precision <- vector("list", K)
  for (iteration in seq_len(max_iter)) {
    responsibility <- state$responsibility
    mass <- colSums(responsibility)
    if (any(mass < 1e-10)) stop("A mean component has collapsed. Try another initialization seed or a larger rho.")
    weights <- mass / n
    means <- crossprod(responsibility, X) / mass
    for (k in seq_len(K)) {
      indices <- which(component_group == k)
      scatter <- matrix(0, p, p)
      for (j in indices) {
        centered <- sweep(X, 2L, means[j, ], "-")
        scatter <- scatter + crossprod(centered, centered * responsibility[, j])
      }
      scatter <- scatter / sum(mass[indices])
      scatter <- (scatter + t(scatter)) / 2 + diag(ridge, p)
      if (rho > 0) {
        estimate <- glasso::glasso(scatter, rho = rho, penalize.diagonal = FALSE, thr = 1e-6)
        if (!is.null(estimate$errflag) && estimate$errflag != 0)
          stop("glasso failed with error flag ", estimate$errflag)
        covariance[[k]] <- estimate$w
        precision[[k]] <- estimate$wi
      } else {
        covariance[[k]] <- scatter
        precision[[k]] <- chol2inv(chol(scatter))
      }
    }
    # Recompute both log likelihood and responsibilities at the returned parameters.
    updated <- posterior_hmm(X, weights, means, covariance, component_group)
    history <- c(history, updated$loglik)
    converged <- abs(updated$loglik - state$loglik) < tol
    state <- updated
    if (converged) break
  }
  group_probability <- vapply(seq_len(K), function(k)
    rowSums(state$responsibility[, component_group == k, drop = FALSE]), numeric(n))
  group_probability <- matrix(group_probability, n, K)
  list(K = K, L = L, rho = rho, weights = weights, means = means,
       covariance = covariance, precision = precision, component_group = component_group,
       responsibility = state$responsibility, group_probability = group_probability,
       labels = max.col(group_probability, ties.method = "first"),
       component_labels = max.col(state$responsibility, ties.method = "first"),
       loglik = state$loglik, loglik_history = history, iterations = iteration,
       converged = converged, seed = seed, ridge = ridge)
}

fit_mm <- function(X, components = 6L, rho = 0, ...)
  fit_hmm(X, K = components, L = 1L, rho = rho, ...)

# Retain the rho heuristic used by the submitted simulation scripts.
# This is NOT the conventional one-standard-error CV rule.
choose_rho <- function(X, seed = 1L) {
  n <- nrow(X); p <- ncol(X)
  folds <- with_seed(seed, cut(sample(n), breaks = 10L, labels = FALSE))
  rho_max <- max(vapply(seq_len(p), function(j)
    max(abs(crossprod(X[, -j, drop = FALSE], X[, j])) / n), numeric(1)))
  if (!is.finite(rho_max) || rho_max <= 0) stop("Degenerate rho grid.")
  grid <- exp(seq(log(rho_max), log(0.0002 * rho_max), length.out = 100L))
  losses <- matrix(NA_real_, length(grid), p)
  for (j in seq_len(p)) {
    cv <- glmnet::cv.glmnet(X[, -j, drop = FALSE], X[, j], foldid = folds, lambda = grid)
    losses[, j] <- cv$cvm
  }
  rho_cv <- grid[which.min(rowMeans(losses))]
  list(rho = rho_cv + mean(cv$cvsd), rho_cv = rho_cv,
       last_variable_cv_sd = mean(cv$cvsd), grid = grid, losses = losses,
       rule = "legacy: rho_cv + mean(last-variable cvsd)")
}

select_hmm <- function(X, rho, K = NULL, max_components = 10L,
                       cut_fraction = 0.6, L_candidates = c(2L, 3L),
                       bic = c("legacy", "full"), seed = 1L,
                       max_iter = 1000L, tol = 1e-6) {
  bic <- match.arg(bic)
  stopifnot(max_components >= 2L, max_components == as.integer(max_components),
            is.finite(cut_fraction), cut_fraction > 0, cut_fraction <= 1,
            all(L_candidates >= 1L), all(L_candidates == as.integer(L_candidates)))
  preliminary <- NULL; tree <- NULL; height <- NA_real_; merged <- NULL
  if (is.null(K)) {
    preliminary <- fit_mm(X, max_components, rho, seed = seed,
                          max_iter = 2L * max_iter, tol = tol)
    distance <- matrix(0, max_components, max_components)
    for (i in seq_len(max_components - 1L)) for (j in (i + 1L):max_components)
      distance[i, j] <- distance[j, i] <- sqrt(sum((preliminary$covariance[[i]] - preliminary$covariance[[j]])^2))
    distance <- stats::as.dist(distance)
    tree <- stats::hclust(distance, method = "ward.D")
    height <- max(distance) * cut_fraction
    merged <- stats::cutree(tree, h = height)
    K <- length(unique(merged))
  }
  fits <- lapply(seq_along(L_candidates), function(i)
    fit_hmm(X, K, L_candidates[i], rho, seed = seed + i,
            max_iter = max_iter, tol = tol))
  p <- ncol(X); n <- nrow(X)
  selection <- do.call(rbind, lapply(seq_along(fits), function(i) {
    J <- sum(fits[[i]]$L)
    legacy_df <- K * p * (p + 1) / 2 + J * p
    data.frame(K = K, L = L_candidates[i], loglik = fits[[i]]$loglik,
               legacy_df = legacy_df, full_df = legacy_df + J - 1L,
               BIC_legacy = -2 * fits[[i]]$loglik + legacy_df * log(n),
               BIC_full = -2 * fits[[i]]$loglik + (legacy_df + J - 1L) * log(n),
               iterations = fits[[i]]$iterations, converged = fits[[i]]$converged)
  }))
  selected <- which.min(selection[[paste0("BIC_", bic)]])
  selection$selected <- seq_len(nrow(selection)) == selected
  list(fit = fits[[selected]], candidates = selection, bic = bic,
       preliminary = preliminary, dendrogram = tree, cut_height = height,
       merged_components = merged)
}
