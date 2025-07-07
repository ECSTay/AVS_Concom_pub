data {
  int<lower=1> N;                         // number of responders
  int<lower=1> N_strat;                   // number of strategies
  int<lower=1> N_sched;                   // number of schedules
  array[N] int<lower=1, upper=N_strat> s; // strategy
  array[N] int<lower=1, upper=N_sched> t; // schedule
  array[N] int<lower=0, upper=1> w;       // sex
  array[N] int<lower=0, upper=1> x;       // Indigenous status
  array[N] int<lower=0, upper=1> z;       // comorbidity
  array[N] int<lower=0, upper=1> y;       // outcome
  real a;                                 // prior distribution mean
  real<lower=0> b;                        // prior distribution standard deviation
}

parameters {
  matrix[N_strat, N_sched] mu; // strategy/schedule log odds
  real beta;                   // sex parameter
  real gamma;                  // Indigenous status parameter
  real delta;                  // comorbidity parameter
}

model {
  for(n_sched in 1:N_sched) target += normal_lpdf(mu[,n_sched] | a, b);
  target += std_normal_lpdf(beta);
  target += std_normal_lpdf(gamma);
  target += std_normal_lpdf(delta);
  for(i in 1:N) target += bernoulli_logit_lpmf(y[i] | mu[s[i], t[i]] + w[i]*beta + x[i]*gamma + z[i]*delta);
}
