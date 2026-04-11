// Optimized Hybrid Hierarchical VSEdelta Model for the IGT
functions {
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
  // 1. NCP Hyperparameters (Loss, Explore Alpha, Explore Bonus)
  vector[3] mu_pr_rest;
  vector<lower=0>[3] sigma_rest;

  // 2. CP Hyperparameters (Con, Gain, Update)
  real<lower=0, upper=5> mu_con;
  real<lower=0> sigma_con;
  real<lower=0, upper=1> mu_gain;
  real<lower=0> sigma_gain;
  real<lower=0, upper=1> mu_update;
  real<lower=0> sigma_update;

  // 3. Individual parameters
  vector[N] loss_pr;
  vector[N] explore_alpha_pr;
  vector[N] explore_bonus_pr;
  
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=1> gain;
  array[N] real<lower=0, upper=1> update;
}

transformed parameters {
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> explore_alpha;
  array[N] real<lower=-10, upper=10> explore_bonus;

  // Non-Centered reparameterization for the rest
  loss          = to_array_1d(inv_logit(mu_pr_rest[1] + sigma_rest[1] * loss_pr) * 10);
  explore_alpha = to_array_1d(inv_logit(mu_pr_rest[2] + sigma_rest[2] * explore_alpha_pr));
  explore_bonus = to_array_1d(-10 + inv_logit(mu_pr_rest[3] + sigma_rest[3] * explore_bonus_pr) * 20);
}

model {
  // Tighter Global Priors for stability
  mu_pr_rest ~ normal(0, 1);
  sigma_rest ~ normal(0, 0.5); 
  
  mu_con ~ normal(2.5, 1);
  sigma_con ~ normal(0, 0.5);
  mu_gain ~ beta(2, 2);
  sigma_gain ~ normal(0, 0.5);
  mu_update ~ beta(2, 2);
  sigma_update ~ normal(0, 0.5);

  loss_pr          ~ normal(0, 1);
  explore_alpha_pr ~ normal(0, 1);
  explore_bonus_pr ~ normal(0, 1);

  // CP Hierarchy for the core parameters
  con    ~ normal(mu_con, sigma_con);
  gain   ~ normal(mu_gain, sigma_gain);
  update ~ normal(mu_update, sigma_update);

  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       con, gain, loss, update, explore_alpha, explore_bonus);
}

generated quantities {
  real mu_loss          = inv_logit(mu_pr_rest[1]) * 10;
  real mu_explore_alpha = inv_logit(mu_pr_rest[2]);
  real mu_explore_bonus = -10 + inv_logit(mu_pr_rest[3]) * 20;
}
