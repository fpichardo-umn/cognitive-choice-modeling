functions {
  vector igt_model_lp(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, vector pers, int Tsub,
        real sensitivity, real update, real gain, real loss, 
        real epP, real epN, real K, real w, real decay
        ) {
    // Define values
    real curUtil;
    vector[4] local_ev = ev;
    vector[4] local_pers = pers;
    vector[4] V;

    // For each trial
    for (t in 1:Tsub) {
      // Calculate combined value
      V = w * local_ev + (1-w) * local_pers;
      
      // Choice probability using sensitivity
      target += categorical_logit_lpmf(choice[t] | sensitivity * V);
      
      // Decay perseverance
      local_pers = local_pers * K;
      
      // Calculate utility & update perseverance
      // Use separate values for wins and losses
      curUtil = pow(wins[t], gain) - loss * pow(losses[t], gain);
      
      // Update perseverance based on net outcome
      if (wins[t] >= losses[t]) {
        local_pers[choice[t]] += epP;
      } else {
        local_pers[choice[t]] += epN;
      }
      
      // First decay all deck values
      local_ev = local_ev * (1 - decay);
      
      // Then update chosen deck with delta rule
      local_ev[choice[t]] += curUtil + update * (curUtil - local_ev[choice[t]]);
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
  // Group-level hyperparameters
  array[9] real mu_pr;                // Group means for all parameters
  array[9] real<lower=0> sigma;       // Group standard deviations

  // Subject-level raw parameters
  array[N] real con_pr;               // Consistency parameter
  array[N] real update_pr;            // Learning rate parameter
  array[N] real gain_pr;              // Outcome sensitivity parameter
  array[N] real loss_pr;              // Loss aversion parameter
  array[N] real epP_pr;               // Positive perseverance strength
  array[N] real epN_pr;               // Negative perseverance strength
  array[N] real K_pr;                 // Perseverance decay parameter
  array[N] real w_pr;                 // Value weight parameter
  array[N] real decay_pr;             // Decay rate
}

transformed parameters {
  // Transform subject-level raw parameters
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=1> update;
  array[N] real<lower=0, upper=1> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real epP;
  array[N] real epN;
  array[N] real<lower=0, upper=1> K;
  array[N] real<lower=0, upper=1> w;
  array[N] real<lower=0, upper=1> decay;
  
  // Hierarchical transformation
  for (n in 1:N) {
    con[n]    = inv_logit(mu_pr[1] + sigma[1] * con_pr[n]) * 5;
    update[n] = inv_logit(mu_pr[2] + sigma[2] * update_pr[n]);
    gain[n]   = inv_logit(mu_pr[3] + sigma[3] * gain_pr[n]);
    loss[n]   = inv_logit(mu_pr[4] + sigma[4] * loss_pr[n]) * 10;
    epP[n]    = mu_pr[5] + sigma[5] * epP_pr[n];
    epN[n]    = mu_pr[6] + sigma[6] * epN_pr[n];
    K[n]      = inv_logit(mu_pr[7] + sigma[7] * K_pr[n]);
    w[n]      = inv_logit(mu_pr[8] + sigma[8] * w_pr[n]);
    decay[n]  = inv_logit(mu_pr[9] + sigma[9] * decay_pr[n]);
  }
}

model {
  // Hyperpriors
  for (i in 1:9) {
    mu_pr[i] ~ normal(0, 1);
    sigma[i] ~ cauchy(0, 2.5);
  }

  // Subject-level priors
  for (n in 1:N) {
    con_pr[n]    ~ normal(0, 1);
    update_pr[n] ~ normal(0, 1);
    gain_pr[n]   ~ normal(0, 1);
    loss_pr[n]   ~ normal(0, 1);
    epP_pr[n]    ~ normal(0, 1);
    epN_pr[n]    ~ normal(0, 1);
    K_pr[n]      ~ normal(0, 1);
    w_pr[n]      ~ normal(0, 1);
    decay_pr[n]  ~ normal(0, 1);
  }

  // Process each subject
  for (n in 1:N) {
    vector[4] ev = rep_vector(0., 4);
    vector[4] pers = rep_vector(0., 4);
    real sensitivity = pow(3, con[n]) - 1;
    
    // Use the same function as single-subject model
    ev = igt_model_lp(choice[n, 1:Tsubj[n]],
                      wins[n, 1:Tsubj[n]], abs(losses[n, 1:Tsubj[n]]), 
                      ev, pers, Tsubj[n], sensitivity, update[n], gain[n], loss[n], 
                      epP[n], epN[n], K[n], w[n], decay[n]);
  }
}

generated quantities {
  // Group-level parameters in interpretable scale
  real<lower=0, upper=5> mu_con;
  real<lower=0, upper=1> mu_update;
  real<lower=0, upper=1> mu_gain;
  real<lower=0, upper=10> mu_loss;
  real mu_epP;
  real mu_epN;
  real<lower=0, upper=1> mu_K;
  real<lower=0, upper=1> mu_w;
  real<lower=0, upper=1> mu_decay;
  
  // Compute interpretable group-level parameters
  mu_con    = inv_logit(mu_pr[1]) * 5;
  mu_update = inv_logit(mu_pr[2]);
  mu_gain   = inv_logit(mu_pr[3]);
  mu_loss   = inv_logit(mu_pr[4]) * 10;
  mu_epP    = mu_pr[5];
  mu_epN    = mu_pr[6];
  mu_K      = inv_logit(mu_pr[7]);
  mu_w      = inv_logit(mu_pr[8]);
  mu_decay  = inv_logit(mu_pr[9]);
}
