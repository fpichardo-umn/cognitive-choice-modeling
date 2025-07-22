functions {
  vector igt_model_lp(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, int Tsub, 
        real sensitivity, real update, real wgt_pun, real wgt_rew
        ) {
    // Define values
    real curUtil;
    vector[4] local_ev = ev;
    
    // For each trial
    for (t in 1:Tsub) {
      // Choice probability using consistency
      target += categorical_logit_lpmf(choice[t] | sensitivity * local_ev);
      
      // Compute utility with separate weights for rewards and punishments
      curUtil = wgt_rew * wins[t] - wgt_pun * losses[t];
      
      // Update expected value of chosen deck
      local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
    }
    
    return local_ev;
  }
}

data {
  int<lower=1> N;                              // Number of subjects
  int<lower=1> T;                              // Maximum number of trials
  array[N] int<lower=1> Tsubj;                 // Number of trials for each subject
  array[N, T] int<lower=1, upper=4> choice;   // Choices made at each trial (1-4)
  array[N, T] real<lower=0> wins;             // Win amount at each trial
  array[N, T] real<lower=0> losses;           // Loss amount at each trial
}

parameters {
  // Group-level hyperparameters - use arrays instead of vectors
  array[4] real mu_pr;                // Group means for all parameters
  array[4] real<lower=0> sigma;       // Group standard deviations

  // Subject-level raw parameters - use arrays instead of vector[N]
  array[N] real con_pr;               // Consistency parameter
  array[N] real wgt_pun_pr;           // Weight for punishments
  array[N] real wgt_rew_pr;           // Weight for rewards
  array[N] real update_pr;            // Updating rate
}

transformed parameters {
  // Transform subject-level raw parameters - use arrays
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=1> wgt_pun;
  array[N] real<lower=0, upper=1> wgt_rew;
  array[N] real<lower=0, upper=1> update;
  
  // Hierarchical transformation - explicit loops instead of vectorized ops
  for (n in 1:N) {
    con[n]     = inv_logit(mu_pr[1] + sigma[1] * con_pr[n]) * 5;
    wgt_pun[n] = inv_logit(mu_pr[2] + sigma[2] * wgt_pun_pr[n]);
    wgt_rew[n] = inv_logit(mu_pr[3] + sigma[3] * wgt_rew_pr[n]);
    update[n]  = inv_logit(mu_pr[4] + sigma[4] * update_pr[n]);
  }
}

model {
  // Hyperpriors - loop instead of vectorized
  for (i in 1:4) {
    mu_pr[i] ~ normal(0, 1);
    sigma[i] ~ cauchy(0, 2.5);
  }

  // Subject-level priors - explicit loops
  for (n in 1:N) {
    con_pr[n]     ~ normal(0, 1);
    wgt_pun_pr[n] ~ normal(0, 1);
    wgt_rew_pr[n] ~ normal(0, 1);
    update_pr[n]  ~ normal(0, 1);
  }

  // Process each subject using the reusable function
  for (n in 1:N) {
    vector[4] ev = rep_vector(0., 4);
    real sensitivity = pow(3, con[n]) - 1;
    
    // Use the same function as single-subject model
    ev = igt_model_lp(choice[n, 1:Tsubj[n]],
		      wins[n, 1:Tsubj[n]], abs(losses[n, 1:Tsubj[n]]), ev, 
                      Tsubj[n], sensitivity, update[n], wgt_pun[n], wgt_rew[n]);
  }
}

generated quantities {
  // Group-level parameters in interpretable scale
  real<lower=0, upper=5> mu_con;
  real<lower=0, upper=1> mu_wgt_pun;
  real<lower=0, upper=1> mu_wgt_rew;
  real<lower=0, upper=1> mu_update;
  
  // Compute interpretable group-level parameters
  mu_con     = inv_logit(mu_pr[1]) * 5;
  mu_wgt_pun = inv_logit(mu_pr[2]);
  mu_wgt_rew = inv_logit(mu_pr[3]);
  mu_update  = inv_logit(mu_pr[4]);
}
