// Optimized Hybrid Hierarchical VSE Model for the IGT
functions {
  real igt_subject(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev_exploit, vector ev_explore, int Tsub,
        real sensitivity, real gain, real loss, real decay, 
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

      local_ev_exploit = local_ev_exploit * (1 - decay);
      local_ev_exploit[choice[t]] += curUtil;
      local_ev_explore[choice[t]] = 0;
      
      for (d in 1:4) {
        if (d != choice[t]) {
          local_ev_explore[d] += explore_alpha * (explore_bonus - local_ev_explore[d]);
        }
      }
    }
    return log_lik;
  }
  
  real partial_sum_func(array[] int slice_n, int start, int end,
                        array[,] int choice, array[,] real wins, array[,] real losses, 
                        array[] int Tsubj,
                        array[] real con, array[] real gain, array[] real loss,
                        array[] real decay, array[] real explore_alpha, 
                        array[] real explore_bonus) {
    real log_lik = 0.0;
    for (n in start:end) {
      vector[4] ev_exploit = rep_vector(0., 4);
      vector[4] ev_explore = rep_vector(0., 4);
      real sensitivity = expm1(log(3) * con[n]);
      
      log_lik += igt_subject(choice[n, 1:Tsubj[n]],
                             wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], 
                             ev_exploit, ev_explore, Tsubj[n], 
                             sensitivity, gain[n], loss[n], decay[n],
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
  // NCP Hyperparameters (Indices: 1=con, 2=loss, 3=alpha, 4=bonus)
  vector[4] mu_pr_rest;
  vector<lower=0>[4] sigma_rest;

  // CP Hyperparameters (Gain & Decay)
  real<lower=0, upper=1> mu_gain;
  real<lower=0> sigma_gain;
  real<lower=0, upper=1> mu_decay;
  real<lower=0> sigma_decay;

  // Individual parameters
  vector[N] con_pr;
  vector[N] loss_pr;
  vector[N] explore_alpha_pr;
  vector[N] explore_bonus_pr;
  
  array[N] real<lower=0, upper=1> gain;   // Centered
  array[N] real<lower=0, upper=1> decay;  // Centered
}

transformed parameters {
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> explore_alpha;
  array[N] real<lower=-10, upper=10> explore_bonus;
  
  // Non-Centered reparameterization
  con           = to_array_1d(inv_logit(mu_pr_rest[1] + sigma_rest[1] * con_pr) * 5);
  loss          = to_array_1d(inv_logit(mu_pr_rest[2] + sigma_rest[2] * loss_pr) * 10);
  explore_alpha = to_array_1d(inv_logit(mu_pr_rest[3] + sigma_rest[3] * explore_alpha_pr));
  explore_bonus = to_array_1d(-10 + inv_logit(mu_pr_rest[4] + sigma_rest[4] * explore_bonus_pr) * 20);
}

model {
  // Tighter Priors for stability
  mu_pr_rest ~ normal(0, 1);
  sigma_rest ~ normal(0, 1); 
  
  mu_gain ~ beta(2, 2); // Informed prior for probability-scale mean
  sigma_gain ~ normal(0, 0.5);
  mu_decay ~ beta(2, 2);
  sigma_decay ~ normal(0, 0.5);

  con_pr           ~ normal(0, 1);
  loss_pr          ~ normal(0, 1);
  explore_alpha_pr ~ normal(0, 1);
  explore_bonus_pr ~ normal(0, 1);

  // Centered Likelihood for Gain and Decay
  gain  ~ normal(mu_gain, sigma_gain);
  decay ~ normal(mu_decay, sigma_decay);

  int grainsize = max(1, N %/% 2); // Larger grainsize for N=100
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       con, gain, loss, decay, explore_alpha, explore_bonus);
}

generated quantities {
  real mu_con           = inv_logit(mu_pr_rest[1]) * 5;
  real mu_loss          = inv_logit(mu_pr_rest[2]) * 10;
  real mu_explore_alpha = inv_logit(mu_pr_rest[3]);
  real mu_explore_bonus = -10 + inv_logit(mu_pr_rest[4]) * 20;
  // mu_gain and mu_decay are already defined in parameters!
}
