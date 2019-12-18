// K: the ploidy of the species.
// N: the number of individuals in the dataset
// x: The vector of reference counts.
// n: The vector of total counts.
data {
  int<lower=0> K;
  int<lower=0> N;
  int x[N];
  int n[N];
}

// logit_eps: The logit-sequencing error rate.
// logit_tau: The logit-overdispersion parameter.
// log_h: The log-bias parameter
// pivec: pivec[k] is the probability of genotype k.
parameters {
  real<lower=0,upper=1> alpha;
  real<upper=0> logit_eps;
  real<upper=0> logit_tau;
  real log_h;
}

// logpostprobmat: logpostprobmat[i, k + 1] is the log of the *unnormalized*
//                 posterior probability that individual i has genotype k.
transformed parameters {
  matrix[N, K + 1] logpostprobmat;
  for (i in 1:N) {
    for (k in 0:K) {
      real logprobk = binomial_lpmf(k | K, alpha);
      real epsilon  = inv_logit(logit_eps);
      real tau      = inv_logit(logit_tau);
      real h        = exp(log_h);
      real p        = (k * 1.0) / K;
      real fi       = p * (1.0 - epsilon) + (1.0 - p) * epsilon;
      real xii      = fi / (h * (1.0 - fi) + fi);
      real alpha_bb = xii * (1.0 - tau) / tau;
      real beta_bb  = (1.0 - xii) * (1.0 - tau) / tau;
      logpostprobmat[i, k + 1] = logprobk + beta_binomial_lpmf(x[i] | n[i], alpha_bb, beta_bb);
    }
  }

}

model {
  // priors
  logit_eps ~ normal(-4.5, 1) T[, 0]; // these upper bounds help identify the model
  logit_tau ~ normal(-5.5, 1) T[, 0];
  log_h     ~ normal(0, 1);
  alpha     ~ uniform(0, 1);

  // likelihood. Integrate out mixing indicators.
  for (i in 1:N) {
    target += log_sum_exp(logpostprobmat[i]);
  }
}

generated quantities {
  matrix<lower=0,upper=1>[N, K + 1] postprobmat;

  for (i in 1:N) {
    postprobmat[i] = softmax(logpostprobmat[i]')';
  }
}
