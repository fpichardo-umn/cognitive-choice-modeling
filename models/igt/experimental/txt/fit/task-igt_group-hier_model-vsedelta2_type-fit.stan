// Optimized Hierarchical VSEdelta Model for the Iowa Gambling Task
functions {
  // Returns log-likelihood instead of modifying target
  real igt_subject(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev_exploit, vector ev_explore, int Tsub,
        real sensitivity, real gain, real loss, real update, 
        real explore_alpha, real explore_bonus
        ) {
    real log_lik = 0.0;
    real curUtil;
    vector[4] local_ev_exploit = ev_exploit;
    vector[4] local_ev_explore = ev_explore;
    vector[4] combined_value;
    
    for (t in 1:Tsub) {
      combined_value = local_ev_exploit + local_ev_explore;
      log_lik += categorical_logit_lpmf(choice[t] | sensitivity * combined_value);
      
      real win_component = (wins[t] == 0) ? 0.0 : exp(gain * log(wins[t]));
real loss_component = (losses[t] == 0) ? 0.0 : exp(gain * log(losses[t]));
curUtil = win_component - loss * loss_component;

      local_ev_exploit[choice[t]] += curUtil + update * (curUtil  - local_ev_exploit[choice[t]]);
      local_ev_explore[choice[t]] = 0;
      
      for (d in 1:4) {
        if (d != choice[t]) {
          local_ev_explore[d] += explore_alpha * (explore_bonus - local_ev_explore[d]);
        }
      }
    }
    
    return log_lik;
  }
  
  // Parallelization wrapper
  real partial_sum_func(array[] int slice_n, int start, int end,
                        array[,] int choice, array[,] real wins, array[,] real losses, 
                        array[] int Tsubj,
                        array[] real con, array[] real gain, array[] real loss,
                        array[] real update,
                        array[] real explore_alpha, array[] real explore_bonus) {
    real log_lik = 0.0;
    
    for (n in start:end) {
      vector[4] ev_exploit = rep_vector(0., 4);
      vector[4] ev_explore = rep_vector(0., 4);
      real sensitivity = expm1(log(3) * con[n]);
      
      log_lik += igt_subject(choice[n, 1:Tsubj[n]],
                             wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], 
                             ev_exploit, ev_explore, Tsubj[n], 
                             sensitivity, gain[n], loss[n], update[n],
                             explore_alpha[n], explore_bonus[n]);
    }
    return log_lik;
  }
}

data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1> Tsubj;
  array[N, T] int<lower=1, upper=4> choice;
  array[N, T] real<lower=0> wins;
  array[N, T] real<lower=0> losses;
}

transformed data {
  array[N] int subject_indices;
  for (i in 1:N) {
    subject_indices[i] = i;
  }
}

parameters {
  // Global means for the NCP parameters (now size 5 since loss is separate)
  vector[5] mu_pr_rest; 
  vector<lower=0>[5] sigma_rest;

  // Centered Parameterization for Loss
  real<lower=0, upper=10> mu_loss;      // The group mean directly on the 0-10 scale
  real<lower=0> sigma_loss;             // The group SD
  vector<lower=0, upper=10>[N] loss;    // Individual loss parameters

  // Individual NCP "pr" parameters
  vector[N] con_pr;
  vector[N] gain_pr;
  vector[N] explore_alpha_pr;
  vector[N] explore_bonus_pr;
  vector[N] update_pr;
}

transformed parameters {
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=1> gain;
  array[N] real<lower=0, upper=1> explore_alpha;
  array[N] real<lower=-10, upper=10> explore_bonus;
  array[N] real<lower=0, upper=1> update;
  
  // Non-Centered Parameterization for the rest
  con           = to_array_1d(inv_logit(mu_pr_rest[1] + sigma_rest[1] * con_pr) * 5);
  gain          = to_array_1d(inv_logit(mu_pr_rest[2] + sigma_rest[2] * gain_pr));
  explore_alpha = to_array_1d(inv_logit(mu_pr_rest[3] + sigma_rest[3] * explore_alpha_pr));
  explore_bonus = to_array_1d(-10 + inv_logit(mu_pr_rest[4] + sigma_rest[4] * explore_bonus_pr) * 20);
  update        = to_array_1d(inv_logit(mu_pr_rest[5] + sigma_rest[5] * update_pr));
}

model {
  // Priors for NCP parameters
  mu_pr_rest ~ normal(0, 1);
  sigma_rest ~ normal(0, 1); // Swapped Student-T for Normal for stability
  
  con_pr ~ normal(0, 1);
  gain_pr ~ normal(0, 1);
  explore_alpha_pr ~ normal(0, 1);
  explore_bonus_pr ~ normal(0, 1);
  update_pr ~ normal(0, 1);

  // Priors/Hierarchy for Centered Loss
  mu_loss ~ uniform(0, 10); // Or a Normal(5, 2) truncated to [0,10]
  sigma_loss ~ normal(0, 1);
  loss ~ normal(mu_loss, sigma_loss); 

  // Likelihood
  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       con, gain, loss, update, explore_alpha, explore_bonus);
}

generated quantities {
  // Re-mapping names for your summary table
  real mu_con           = inv_logit(mu_pr_rest[1]) * 5;
  real mu_gain          = inv_logit(mu_pr_rest[2]);
  real mu_explore_alpha = inv_logit(mu_pr_rest[3]);
  real mu_explore_bonus = -10 + inv_logit(mu_pr_rest[4]) * 20;
  real mu_update        = inv_logit(mu_pr_rest[5]);
  
  // mu_loss is already calculated in the parameters block!
}
