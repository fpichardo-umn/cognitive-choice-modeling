// Optimized Hybrid Hierarchical VSEBoth Model for the IGT
functions {
  real igt_subject(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev_exploit, vector ev_explore, int Tsub,
        real sensitivity, real gain, real loss, real update, real decay, 
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
      local_ev_exploit[choice[t]] += curUtil + update * (curUtil - local_ev_exploit[choice[t]]);
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
                        array[] real update, array[] real decay, 
                        array[] real explore_alpha, array[] real explore_bonus) {
    real log_lik = 0.0;
    for (n in start:end) {
      vector[4] ev_exploit = rep_vector(0., 4);
      vector[4] ev_explore = rep_vector(0., 4);
      real sensitivity = expm1(log(3) * con[n]);
      
      log_lik += igt_subject(choice[n, 1:Tsubj[n]],
                             wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], 
                             ev_exploit, ev_explore, Tsubj[n], 
                             sensitivity, gain[n], loss[n], update[n], decay[n],
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
  // NCP Hyperparameters for Update
  real mu_pr_update;
  real<lower=0> sigma_update;

  // CP Hyperparameters (Means and SDs)
  real<lower=0, upper=5> mu_con;
  real<lower=0> sigma_con;
  real<lower=0, upper=1> mu_gain;
  real<lower=0> sigma_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0> sigma_loss;
  real<lower=0, upper=1> mu_decay;
  real<lower=0> sigma_decay;
  real<lower=0, upper=1> mu_explore_alpha;
  real<lower=0> sigma_explore_alpha;
  real<lower=-10, upper=10> mu_explore_bonus;
  real<lower=0> sigma_explore_bonus;

  // Individual parameters
  vector[N] update_pr;                           // Non-Centered [cite: 53]
  array[N] real<lower=0, upper=5> con;           // Centered [cite: 54]
  array[N] real<lower=0, upper=1> gain;          // Centered [cite: 54]
  array[N] real<lower=0, upper=10> loss;         // Centered [cite: 54]
  array[N] real<lower=0, upper=1> decay;         // Centered [cite: 54]
  array[N] real<lower=0, upper=1> explore_alpha; // Centered [cite: 54]
  array[N] real<lower=-10, upper=10> explore_bonus; // Centered [cite: 54]
}

transformed parameters {
  array[N] real<lower=0, upper=1> update;
  // Non-Centered reparameterization for Update [cite: 60]
  update = to_array_1d(inv_logit(mu_pr_update + sigma_update * update_pr));
}

model {
  // 1. Tightened Global Priors
  mu_pr_update ~ normal(0, 1);
  sigma_update ~ normal(0, 0.5); // Tightened 
  
  // Weakly informative priors for Centered means
  mu_con    ~ normal(2.5, 1);
  mu_gain   ~ beta(2, 2);
  mu_loss   ~ normal(5, 2);
  mu_decay  ~ beta(2, 2);
  mu_explore_alpha ~ beta(2, 2);
  mu_explore_bonus ~ normal(0, 5);

  sigma_con           ~ normal(0, 1);
  sigma_gain          ~ normal(0, 0.5);
  sigma_loss          ~ normal(0, 1);
  sigma_decay         ~ normal(0, 0.5);
  sigma_explore_alpha ~ normal(0, 0.5);
  sigma_explore_bonus ~ normal(0, 1);

  // 2. Hierarchy [cite: 62-66]
  update_pr ~ normal(0, 1);
  con           ~ normal(mu_con, sigma_con);
  gain          ~ normal(mu_gain, sigma_gain);
  loss          ~ normal(mu_loss, sigma_loss);
  decay         ~ normal(mu_decay, sigma_decay);
  explore_alpha ~ normal(0.5, 0.2); // Alternative informative prior
  explore_bonus ~ normal(mu_explore_bonus, sigma_explore_bonus);

  // 3. Likelihood [cite: 67]
  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       con, gain, loss, update, decay, explore_alpha, explore_bonus);
}

generated quantities {
  // mu_con, mu_gain, mu_loss, mu_decay, mu_explore_alpha, mu_explore_bonus are in parameters
  real mu_update = inv_logit(mu_pr_update);
}
