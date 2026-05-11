data {  
  int<lower=1> N_R;                          // number of responses
  int<lower=1> N_I;                          // number of responders
  int<lower=1> N_sched;                      // number of schedules = 4
  array[N_R] int<lower=0, upper=N_I> infant; // specific infant
  array[N_R] int<lower=0, upper=1> s;        // strategy
  array[N_R] int<lower=1, upper=N_sched> t;  // schedule timepoint
  array[N_R] int<lower=1, upper=3> y;        // outcomes in {1, 2, 3}, corresponding to categories {0, 1, 2}
}

parameters {
matrix[N_sched, 2] mu;                     // separate strategy log odd ratios compared to no event for each timepoint
  //real alpha_star;                           // hierarchical mean log odds ratio of concomitant compared to separate for a single event
  //real sigma_alpha;                          // hierarchical standard deviation
  vector[N_sched] alpha;                  // log odds ratio of concomitant compared to separate for a single event for each timepointlog odds ratio of concomitant compared to separate for a single event for each timepoint
 //vector[2] sigma_epsilon;                   // infant random effect standard deviation
 //matrix[N_I, 2] epsilon_raw;                // temp variables epsilons
}

transformed parameters {
  array[N_sched] real alpha;                 // log odds ratio of concomitant compared to separate for a single event for each timepoint
 //matrix[N_I, 2] epsilon;                    // infant random effect
  array[N_R] vector[3] eta;                  // linear predictors
  array[N_R] vector[3] p;                    // probabilities transformed from linear predictors
  
  for(tt in 1:N_sched) alpha[tt] = alpha_star + sigma_alpha*alpha_raw[tt];
  //for(j in 1:2) epsilon[,j] = sigma_epsilon[j]*epsilon_raw[,j];
  for(r in 1:N_R){
    eta[r][1] = 0;
    eta[r][2] = mu[t[r], 1] + s[r]*alpha[t[r]];
    if(s[r] == 0){
      eta[r][3] = mu[t[r], 2];
    } else {
      eta[r][3] = negative_infinity();
    }
    p[r] = softmax(eta[r]);
  }
}

model {
  for(tt in 1:N_sched) target += normal_lpdf(mu[tt,] | 0, 2);
  target += std_normal_lpdf(alpha);
  target += inv_gamma_lpdf(sigma_alpha | 3, 1);
 target += std_normal_lpdf(alpha_raw);
  //target += inv_gamma_lpdf(sigma_epsilon | 3, 1);
  //for(j in 1:2) target += std_normal_lpdf(epsilon_raw[,j]);
  for(r in 1:N_R) target += categorical_lpmf(y[r] | p[r]);
}
