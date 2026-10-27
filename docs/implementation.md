# Implementation and reproducibility notes

## Model updates

For component j = (k,l), responsibilities are proportional to `pi[k] * eta[k,l] * Normal(x; mu[k,l], Sigma[k])`. The implementation stores their joint weights and normalizes the log densities with a row-wise log-sum-exp calculation. Updating these joint weights is equivalent to updating the top-level and conditional weights separately.

For each network k, update the component means and form

```
S[k] = sum_l sum_i gamma[i,k,l] *
       (x[i] - mu[k,l]) (x[i] - mu[k,l])' / sum_i sum_l gamma[i,k,l]
```

Pass the normalized scatter to `glasso::glasso(S, rho, penalize.diagonal=FALSE)`. All L means in that group then use the same returned covariance. Setting L = 1 gives the flat mixture; setting rho = 0 skips glasso. Rho here is the penalty on the normalized scatter, as in the original implementation. No global penalized-likelihood monotonicity or statistical consistency is asserted by the code.

The runnable implementation has a small diagonal ridge (1e-8) for numerical stability, errors on collapsed mean components, uses vectorized log densities, and recomputes responsibilities after the final parameter update. The glasso inner-solver threshold is tightened to 1e-6 to improve consistency between its returned covariance and precision matrices. A fixed seed controls each fit. Initial group allocations are balanced random allocations, followed by within-group mean-subcluster initialization. This differs slightly from the original flat-mixture random allocation with replacement. Convergence uses the original absolute observed-log-likelihood change tolerance of 1e-6; a maximum-iteration fit is flagged, not reported as converged.

## Graph construction

Both graph generators use seed 123 by default and preserve the caller's RNG state. Data and initialization have separate seeds. This replaces the original `set.seed(NULL)` reset and makes the demo repeatable.

- **Chain (original mod 5):** start at location s[1] = 1, add independent Uniform(0.5,1) gaps, and set covariance[i,j] = exp(-0.1 * abs(s[i]-s[j])). Its inverse has chain support. Group 2 removes precision edges (2,3) and (3,4); group 3 removes (5,6) and (7,8). Recompute covariance by inversion.
- **Random (original mod 8):** sample exactly p-1 unordered pairs uniformly, draw weights Uniform(0.3,0.8), and set Omega = diag(rowSums(W) + 0.01) - W. This is a fixed-edge-count random graph, rather than independent Bernoulli edges. Rank existing edges by absolute partial correlation. Group 2 removes ranks 1 and 4; group 3 removes ranks 2 and 3. Keep diagonal entries and invert the altered precisions. For p=10 and seed 123, the removals are (2,3)/(3,9) and (1,10)/(7,10), matching the displayed paper example.

The demo's adjacency matrices are derived directly from the precision matrices, with a 1e-8 tolerance for inversion roundoff. The original random generator has a duplicated adjacency assignment that can leave the stored group-3 adjacency asymmetric; its generated covariance is unaffected. Using precision support fixes that plotting/truth-bookkeeping issue. Tests compare all generated covariances with the original generator for p = 10, 20, and 30.

Each graph has three true network groups and two mean subclusters per group, with n/6 observations per subcluster. All coordinates in the six means equal 0, delta, 2*delta, 3*delta, 4*delta, and 5*delta. The demo requires n divisible by six and p >= 8. It uses unstandardized simulated observations, as in the original scripts.

## Selection rules

The initial MM-glasso has 10 components. Pairwise distances are between **covariances**, not precisions. `ward.D` and the 0.6*max(pairwise distance) cut are retained exactly as selection rules. Only the number of resulting groups is used to initialize the HMM dimensions. The tree's grouping is not used as an HMM membership prior.

For this one dataset, choose common L from {2,3}, conditional on that selected K. The default `legacy` BIC counts

```
df_legacy = K * p * (p+1)/2 + K * L * p
BIC_legacy = -2 * observed_loglik + df_legacy * log(n)
```

The optional `full` score adds K*L-1 free joint mixing weights. Both scores count dense covariance parameters; neither is a sparse degrees-of-freedom criterion or a proof of selection consistency. Both values, candidate convergence flags, and the selected candidate are exported. Candidates are compared using their available fits, even if a candidate reaches its iteration limit; inspect these flags before interpreting such a selection.

The rho heuristic computes a 100-value log-spaced grid from rho_max to 0.0002*rho_max, where rho_max is the maximum absolute uncentered cross-product divided by n. Use shared 10-fold splits for nodewise `cv.glmnet`; choose the grid value with minimum mean nodewise CV loss, then add the mean `cvsd` from the **last node's** CV fit. This unusual final addition is deliberately preserved for traceability and is not the conventional one-standard-error rule. The tuning object is saved in each RDS.

## Evaluation

The primary assignment display is ARI against the three true network groups. Six-mean-subcluster metrics are separate columns for every method. Purity is sum of the largest count in each predicted-cluster row divided by n; it allows many fitted groups to match one true group and should not be called optimally permuted one-to-one accuracy.

For edge metrics, assign each observation the precision matrix of its fitted network group. Count TP, TN, FP, and FN over observations and unique off-diagonal pairs. Sensitivity = TP/(TP+FN) and specificity = TN/(TN+FP). This rule needs no label mapping, handles extra fitted groups, and avoids overwriting scores when multiple fitted groups correspond to one true group. It measures both membership and network error; it is not an oracle comparison of three correctly matched precision estimates.

The flat methods each receive the true total component count (6), whereas HMM uses the documented K/L heuristics unless K is supplied. This small showcase is not an equally tuned benchmark. Results from one dataset cannot substantiate claims of improved average assignment or edge recovery.

## Scope of validation

The default demo fits one chain dataset and one random dataset. The graph tests use deterministic matrix construction without Monte Carlo fitting. The remaining checks verify label-invariant metrics and reuse saved fits when available. No 200-replicate study is run or included.
