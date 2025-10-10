data {
  int<lower=1> 			    N; 	      // Number of subjects
  int<lower=1> 			    T;        // Number of trials
  array[N] int<lower=1> 	    sid;      // Subject IDs
  array[N] int<lower=1> 	    Tsubj;    // Number of trials for a subject
  array[N] real<lower=0> 	    minRT;    // Minimum RT per subject + small value to restrict tau
  real 				    RTbound;  // Lower bound or RT across all subjects (e.g., 0.1 second)
  array[N, T] real<lower=0> 	    RT;       // Reaction times
  array[N, T] int<lower=0, upper=1> choice;   // Binary choices (1=play, 0=pass)
}

parameters {
  // Group-level hyperparameters
  array[6] real mu_pr;                // Group means for all parameters
  array[6] real<lower=0> sigma;       // Group standard deviations

  // Subject-level raw parameters
  array[N] real<lower=-5, upper=5> boundary1_pr;  // Boundary separation (T1 a)
  array[N] real<lower=-5, upper=5> boundary_pr;   // Boundary separation (a)
  array[N] real tau1_pr;       // Non-decision time (T1 tau)
  array[N] real tau_pr;        // Non-decision time (tau)
  array[N] real beta_pr;       // Starting point
  array[N] real drift;         // Drift rate
}

transformed parameters {
  array[N] real<lower=0, upper=6> 		      boundary1;
  array[N] real<lower=0, upper=6> 		      boundary;
  array[N] real<lower=RTbound, upper=max(minRT)> tau1;
  array[N] real<lower=RTbound, upper=max(minRT)> tau;
  array[N] real<lower=0, upper=1> 	      beta;

  // Hierarchical transformation - explicit loops
  for (n in 1:N) {
    boundary1[n] = inv_logit(mu_pr[1] + sigma[1] * boundary1_pr[n]) * 5 + 0.01;
    boundary[n]  = inv_logit(mu_pr[2] + sigma[2] * boundary_pr[n]) * 5 + 0.01;
    tau1[n]      = inv_logit(mu_pr[3] + sigma[3] * tau1_pr[n]) * (minRT[n] - RTbound - 1e-6) * 0.99 + RTbound;
    tau[n]       = inv_logit(mu_pr[4] + sigma[4] * tau_pr[n]) * (minRT[n] - RTbound - 1e-6) * 0.99 + RTbound;
    beta[n]      = inv_logit(mu_pr[5] + sigma[5] * beta_pr[n]);
  }
}

model {
  // Hyperpriors
  for (i in 1:6) {
    mu_pr[i] ~ normal(0, 1);
    sigma[i] ~ normal(0, 2);
  }

  // Subject-level priors
  for (n in 1:N) {
    boundary1_pr[n] ~ normal(0, 1);
    boundary_pr[n]  ~ normal(0, 1);
    tau1_pr[n]      ~ normal(0, 1);
    tau_pr[n]       ~ normal(0, 1);
    beta_pr[n]      ~ normal(0, 1);
    drift[n]        ~ normal(mu_pr[6], sigma[6]);
  }
  
  // For each subject
  for (n in 1:N) {
    // Accumulate play/pass info with valid RTs
    array[Tsubj[n]] int play_indices;
    array[Tsubj[n]] int pass_indices;
    int play_count = 0;
    int pass_count = 0;
    
    // Separate by block and choice
    vector[Tsubj[n]] boundaries;
    vector[Tsubj[n]] nondt;
    
    for (t in 1:Tsubj[n]) {
      // Set block-specific boundary and tau
      if (t <= 20) {
        boundaries[t] = boundary1[n];
        nondt[t] = tau1[n];
      } else {
        boundaries[t] = boundary[n];
        nondt[t] = tau[n];
      }
      
      // Store indices for play and pass - ONLY for valid RTs
      if (RT[n, t] != 999) {
        if (choice[n, t] == 1) {
          play_count += 1;
          play_indices[play_count] = t;
        } else {
          pass_count += 1;
          pass_indices[pass_count] = t;
        }
      }
    }
    
    // Likelihood
    RT[n, play_indices[:play_count]] ~ wiener(boundaries[play_indices[:play_count]], nondt[play_indices[:play_count]], beta[n], drift[n]);
    RT[n, pass_indices[:pass_count]] ~ wiener(boundaries[pass_indices[:pass_count]], nondt[pass_indices[:pass_count]], 1-beta[n], -drift[n]);
  }
}

generated quantities {
  // Group-level parameters in interpretable scale
  real<lower=0, upper=6>     mu_boundary1;
  real<lower=0, upper=6>     mu_boundary;
  real<lower=0>              mu_tau1;
  real<lower=0>              mu_tau;
  real<lower=0, upper=1>     mu_beta;
  real                       mu_drift;
  
  mu_boundary1 = inv_logit(mu_pr[1]) * 5 + 0.01;
  mu_boundary  = inv_logit(mu_pr[2]) * 5 + 0.01;
  mu_tau1      = inv_logit(mu_pr[3]) * (mean(minRT) - RTbound) * 0.99 + RTbound;
  mu_tau       = inv_logit(mu_pr[4]) * (mean(minRT) - RTbound) * 0.99 + RTbound;
  mu_beta      = inv_logit(mu_pr[5]);
  mu_drift     = mu_pr[6];
}
