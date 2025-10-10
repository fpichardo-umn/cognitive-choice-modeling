functions {
  vector igt_model_lp(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, int Tsub,
        real sensitivity, real gain, real loss, real update
        ) {
    // Define values
    real curUtil;
    vector[4] local_ev = ev;
    
    // For each trial
    for (t in 1:Tsub) {
      // Choice probability
      target += categorical_logit_lpmf(choice[t] | sensitivity * local_ev);
      
      // Compute utility from wins and losses separately
      curUtil = pow(wins[t], gain) - loss * pow(losses[t], gain);
      
      // Update expected value using delta rule
      local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
    }
    
    return local_ev;
  }
}

data {
  int<lower=1> N;                              // Number of subjects
  int<lower=1> T;                              // Maximum number of trials
  array[N] int<lower=1> sid;      	       // Subject IDs
  array[N] int<lower=1> Tsubj;                 // Number of trials for each subject
  array[N, T] int<lower=1, upper=4> choice;   // Choices made at each trial (1-4)
  array[N, T] real<lower=0> wins;             // Win amount at each trial
  array[N, T] real<lower=0> losses;           // Loss amount at each trial
}

parameters {
  // Group-level hyperparameters
  array[4] real mu_pr;                // Group means for all parameters
  array[4] real<lower=0> sigma;       // Group standard deviations

  // Subject-level raw parameters
  array[N] real con_pr;               // Consistency parameter
  array[N] real gain_pr;              // Value sensitivity
  array[N] real loss_pr;              // Loss aversion
  array[N] real update_pr;            // Updating rate
}

transformed parameters {
  // Transform subject-level raw parameters
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=2> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> update;
  
  // Hierarchical transformation
  for (n in 1:N) {
    con[n]    = inv_logit(mu_pr[1] + sigma[1] * con_pr[n]) * 5;
    gain[n]   = inv_logit(mu_pr[2] + sigma[2] * gain_pr[n]) * 2;
    loss[n]   = inv_logit(mu_pr[3] + sigma[3] * loss_pr[n]) * 10;
    update[n] = inv_logit(mu_pr[4] + sigma[4] * update_pr[n]);
  }
}

model {
  // Hyperpriors
  for (i in 1:4) {
    mu_pr[i] ~ normal(0, 1);
    sigma[i] ~ normal(0, 2);
  }

  // Subject-level priors
  for (n in 1:N) {
    con_pr[n]    ~ normal(0, 1);
    gain_pr[n]   ~ normal(0, 1);
    loss_pr[n]   ~ normal(0, 1);
    update_pr[n] ~ normal(0, 1);
  }

  // Process each subject
  for (n in 1:N) {
    vector[4] ev = rep_vector(0., 4);
    real sensitivity = pow(3, con[n]) - 1;
    
    // Use the same function as single-subject model
    ev = igt_model_lp(choice[n, 1:Tsubj[n]],
                      wins[n, 1:Tsubj[n]], abs(losses[n, 1:Tsubj[n]]), ev, 
                      Tsubj[n], sensitivity, gain[n], loss[n], update[n]);
  }
}

generated quantities {
  // Group-level parameters in interpretable scale
  real<lower=0, upper=5> mu_con;
  real<lower=0, upper=2> mu_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0, upper=1> mu_update;
  
  // Compute interpretable group-level parameters
  mu_con    = inv_logit(mu_pr[1]) * 5;
  mu_gain   = inv_logit(mu_pr[2]) * 2;
  mu_loss   = inv_logit(mu_pr[3]) * 10;
  mu_update = inv_logit(mu_pr[4]);
}
