data {
  int<lower=1> N;                         // number of responders
  int<lower=1> N_strat;                   // number of strategies
  int<lower=1> N_sched;                   // number of schedules
  array[N] int<lower=1, upper=N_strat> s; // strategy
  array[N] int<lower=1, upper=N_sched> t; // schedule
  array[N] int<lower=0, upper=1> w;       // sex
  array[N] int<lower=0, upper=1> x;       // Indigenous status
  array[N] int<lower=0, upper=1> z;       // comorbidity
  array[N] int<lower=1, upper=3> y;       // outcomes in {1, 2, 3}, corresponding to categories {0, 1, 2}
  real a;                                 // prior distribution mean for one event
  real<lower=0> b;                        // prior distribution standard deviation for one event
  real c;                                 // prior distribution mean for two events
  real<lower=0> d;                        // prior distribution standard deviation for two events
}

parameters {
  array[N_strat, N_sched, 2] real mu; // strategy/schedule log odds
  vector[2] beta;                     // sex parameter
  vector[2] gamma;                    // Indigenous status parameter
  vector[2] delta;                    // comorbidity parameter
}

transformed parameters {
  array[N] vector[3] eta; // responder and outcome level linear predictors
  array[N] vector[3] p;   // responder and outcome level probabilities

  for(i in 1:N){
    eta[i][1] = 0;
    eta[i][2] = mu[s[i], t[i], 1] + w[i]*beta[1] + x[i]*gamma[1] + z[i]*delta[1];
    if(s[i] == 1){
      eta[i][3] = negative_infinity();
    } else {
      eta[i][3] = mu[s[i], t[i], 2] + w[i]*beta[2] + x[i]*gamma[2] + z[i]*delta[2];
    }
    p[i] = softmax(eta[i,]);
  }
}

model {
  for(n_sched in 1:N_sched) target += normal_lpdf(mu[, n_sched, 1] | a, b);
  for(n_sched in 1:N_sched) target += normal_lpdf(mu[, n_sched, 2] | c, d);
  target += std_normal_lpdf(beta);
  target += std_normal_lpdf(gamma);
  target += std_normal_lpdf(delta);
  for(i in 1:N) target += categorical_lpmf(y[i] | p[i]);
}
