data {  
  int<lower=1> N_R;                          // number of responses
  int<lower=1> N_I;                          // number of responders
  int<lower=1> N_sched;                      // number of schedules = 4
  int<lower=1> N_state;                      // number of states excluding NSW = 6
  int<lower=1> N_clinic;                     // number of clinic types excluding GP = 2
  array[N_R] int<lower=0, upper=N_I> infant; // specific infant
  array[N_R] int<lower=0, upper=1> s;        // strategy
  array[N_R] int<lower=1, upper=N_sched> t;  // schedule timepoint
  array[N_R] int<lower=0, upper=1> w;        // sex
  array[N_R] int<lower=0, upper=1> x;        // Indigenous status
  array[N_R] row_vector[N_state] q;          // jurisdiction
  array[N_R] row_vector[N_clinic] c;         // clinic type  
  array[N_R] int<lower=0, upper=1> z;        // medical conditions
  array[N_R] int<lower=1, upper=3> y;        // outcomes in {1, 2, 3}, corresponding to categories {0, 1, 2}
}

parameters {
  matrix[N_sched, 2] mu;                     // separate strategy log odd ratios compared to no event for each timepoint
  real alpha_star;                           // hierarchical mean log odds ratio of concomitant compared to separate for a single event
  real sigma_alpha;                          // hierarchical standard deviation
  vector[N_sched] alpha_raw;                 // temp variables for alphas
  vector[2] sigma_epsilon;                   // infant random effect standard deviation
  matrix[N_I, 2] epsilon_raw;                // temp variables epsilons
  vector[2] beta;                            // sex
  vector[2] gamma;                           // Indigenous status
  matrix[N_state, 2] rho;                    // clinic state
  matrix[N_clinic, 2] tau;                   // clinic type
  vector[2] delta;                           // medical conditions
}

transformed parameters {
  array[N_sched] real alpha;                 // log odds ratio of concomitant compared to separate for a single event for each timepoint
  matrix[N_I, 2] epsilon;                    // infant random effect
  array[N_R] vector[3] eta;                  // linear predictors
  array[N_R] vector[3] p;                    // probabilities transformed from linear predictors
  
  for(tt in 1:N_sched) alpha[tt] = alpha_star + sigma_alpha*alpha_raw[tt];
  for(j in 1:2) epsilon[,j] = sigma_epsilon[j]*epsilon_raw[,j];
  for(r in 1:N_R){
    eta[r][1] = 0;
    eta[r][2] = mu[t[r], 1] + s[r]*alpha[t[r]] + epsilon[infant[r],1] + w[r]*beta[1] + x[r]*gamma[1] + q[r]*rho[,1] + c[r]*tau[,1] + z[r]*delta[1];
    if(s[r] == 0){
      eta[r][3] = mu[t[r], 2] + epsilon[infant[r],2] + w[r]*beta[2] + x[r]*gamma[2] + q[r]*rho[,2] + c[r]*tau[,2] + z[r]*delta[2];
    } else {
      eta[r][3] = negative_infinity();
    }
    p[r] = softmax(eta[r]);
  }
}

model {
  for(tt in 1:N_sched) target += normal_lpdf(mu[tt,] | 0, 2);
  target += std_normal_lpdf(alpha_star);
  target += inv_gamma_lpdf(sigma_alpha | 3, 1);
  target += std_normal_lpdf(alpha_raw);
  target += inv_gamma_lpdf(sigma_epsilon | 3, 1);
  for(j in 1:2) target += std_normal_lpdf(epsilon_raw[,j]);
  target += std_normal_lpdf(beta);
  target += std_normal_lpdf(gamma);
  for(j in 1:2) target += std_normal_lpdf(rho[,j]);
  for(j in 1:2) target += std_normal_lpdf(tau[,j]);
  target += std_normal_lpdf(delta);
  for(r in 1:N_R) target += categorical_lpmf(y[r] | p[r]);
}
