// gravity_model.stan
// Hierarchical multi-lineage gravity coupling on log-growth transitions.
// Pooling: theta1, theta3 SHARED across lineages (full pooling);
//          r, theta0 lineage-specific; sigma, nu shared.
// Kernel: exp(theta1*(p_j + p_k) - theta3*d_jk), p/d centered (04).
// Likelihood: Student-t (robust to heavy tails from near-floor transitions).

data {
  int<lower=1> n_obs;
  int<lower=1> n_loc;
  int<lower=1> n_lineage;
  array[n_obs] int<lower=1, upper=n_lineage> L;   // lineage per obs
  array[n_obs] int<lower=1, upper=n_loc> j;        // focal location per obs
  vector[n_obs] log_growth;                        // response
  matrix[n_obs, n_loc] neighbor_abund;             // donor abundance at t; self/unsampled = 0
  vector[n_loc] log_pop_c;                         // centered log-population
  matrix[n_loc, n_loc] log_dist_c;                 // centered log-distance (offdiag)
}

parameters {
  vector[n_lineage] r;                             // lineage-specific baseline
  vector[n_lineage] theta0;                        // lineage-specific coupling magnitude
  real theta1;                                     // SHARED mass exponent
  real theta3;                                     // SHARED distance-decay exponent
  real<lower=0> sigma;                             // shared Student-t scale
  real<lower=1> nu;                                // shared Student-t dof
}

transformed parameters {
  vector[n_obs] mu;
  for (n in 1:n_obs) {
    int jn = j[n];
    int ln = L[n];
    real gsum = 0;
    for (k in 1:n_loc) {
      if (neighbor_abund[n, k] > 0) {              // skips self (=0) and unsampled (=0)
        real logK = theta1 * (log_pop_c[jn] + log_pop_c[k])
                    - theta3 * log_dist_c[jn, k];
        gsum += exp(logK) * neighbor_abund[n, k];
      }
    }
    mu[n] = r[ln] + theta0[ln] * gsum;
  }
}

model {
  // Priors (centered covariates make unit-scale priors weakly informative).
  r      ~ normal(0, 1);                           // vectorized over lineages
  theta0 ~ normal(0, 1);                           // vectorized over lineages
  theta1 ~ normal(0, 1);
  theta3 ~ normal(0, 1);
  sigma  ~ exponential(1);
  nu     ~ gamma(2, 0.1);

  log_growth ~ student_t(nu, mu, sigma);
}

generated quantities {
  // Pointwise log-lik for LOO/PSIS model checking later.
  vector[n_obs] log_lik;
  for (n in 1:n_obs)
    log_lik[n] = student_t_lpdf(log_growth[n] | nu, mu[n], sigma);
}
