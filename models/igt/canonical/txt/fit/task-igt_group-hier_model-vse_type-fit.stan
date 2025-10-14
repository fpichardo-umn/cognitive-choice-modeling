functions {
  vector igt_model_lp(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev_exploit, vector ev_explore, int Tsub,
        real sensitivity, real gain, real loss, real decay, 
        real explore_alpha, real explore_bonus
        ) {
    // Define values
    real curUtil;
    vector[4] local_ev_exploit = ev_exploit;
    vector[4] local_ev_explore = ev_explore;
    vector[4] combined_value;
    
    // For each trial
    for (t in 1:Tsub) {
      // Combine exploitation and exploration values
      combined_value = local_ev_exploit + local_ev_explore;
      
      // Choice probability using sensitivity
      target += categorical_logit_lpmf(choice[t] | sensitivity * combined_value);
      
      // Calculate utility (using value sensitivity for wins and losses)
      curUtil = pow(wins[t], gain) - loss * pow(losses[t], gain);
      
      // Exploitation: Decay all deck values
      local_ev_exploit = local_ev_exploit * (1 - decay);
      
      // Exploitation: Update chosen deck
      local_ev_exploit[choice[t]] += curUtil;
      
      // Exploration: Reset chosen deck to zero
      local_ev_explore[choice[t]] = 0;
      
      // Exploration: Update unchosen decks (return to exploration bonus)
      for (d in 1:4) {
        if (d != choice[t]) {
          local_ev_explore[d] += explore_alpha * (explore_bonus - local_ev_explore[d]);
        }
      }
    }
    
    return local_ev_exploit;
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
  array[6] real mu_pr;                // Group means for all parameters
  array[6] real<lower=0> sigma;       // Group standard deviations

  // Subject-level raw parameters
  array[N] real con_pr;               // Consistency parameter
  array[N] real gain_pr;              // Value sensitivity parameter
  array[N] real loss_pr;              // Loss aversion
  array[N] real decay_pr;             // Decay parameter for exploitation
  array[N] real explore_alpha_pr;     // Learning rate for exploration
  array[N] real explore_bonus_pr;     // Exploration bonus parameter
}

transformed parameters {
  // Transform subject-level raw parameters
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=1> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> decay;
  array[N] real<lower=0, upper=1> explore_alpha;
  array[N] real<lower=-10, upper=10> explore_bonus;
  
  // Hierarchical transformation
  for (n in 1:N) {
    con[n]           = inv_logit(mu_pr[1] + sigma[1] * con_pr[n]) * 5;
    gain[n]          = inv_logit(mu_pr[2] + sigma[2] * gain_pr[n]);
    loss[n]          = inv_logit(mu_pr[3] + sigma[3] * loss_pr[n]) * 10;  // Fixed bug: was using con_pr
    decay[n]         = inv_logit(mu_pr[4] + sigma[4] * decay_pr[n]);
    explore_alpha[n] = inv_logit(mu_pr[5] + sigma[5] * explore_alpha_pr[n]);
    explore_bonus[n] = -10 + inv_logit(mu_pr[6] + sigma[6] * explore_bonus_pr[n]) * 20;
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
    con_pr[n]           ~ normal(0, 1);
    gain_pr[n]          ~ normal(0, 1);
    loss_pr[n]          ~ normal(0, 1);
    decay_pr[n]         ~ normal(0, 1);
    explore_alpha_pr[n] ~ normal(0, 1);
    explore_bonus_pr[n] ~ normal(0, 1);
  }

  // Process each subject
  for (n in 1:N) {
    vector[4] ev_exploit = rep_vector(0., 4);
    vector[4] ev_explore = rep_vector(0., 4);
    real sensitivity = pow(3, con[n]) - 1;
    
    // Use the same function as single-subject model
    ev_exploit = igt_model_lp(choice[n, 1:Tsubj[n]],
                              wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], 
                              ev_exploit, ev_explore, Tsubj[n], 
				sensitivity, gain[n],loss[n], decay[n],
				explore_alpha[n], explore_bonus[n]);
  }
}

generated quantities {
  // Group-level parameters in interpretable scale
  real<lower=0, upper=5> mu_con;
  real<lower=0, upper=1> mu_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0, upper=1> mu_decay;
  real<lower=0, upper=1> mu_explore_alpha;
  real<lower=-10, upper=10> mu_explore_bonus;
  
  // Compute interpretable group-level parameters
  mu_con           = inv_logit(mu_pr[1]) * 5;
  mu_gain          = inv_logit(mu_pr[2]);
  mu_loss          = inv_logit(mu_pr[3]) * 10;
  mu_decay         = inv_logit(mu_pr[4]);
  mu_explore_alpha = inv_logit(mu_pr[5]);
  mu_explore_bonus = -10 + inv_logit(mu_pr[6]) * 20;
}
