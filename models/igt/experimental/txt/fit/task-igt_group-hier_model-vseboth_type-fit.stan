// Optimized Hierarchical VSEBoth Model for the Iowa Gambling Task
functions {
  // Returns log-likelihood instead of modifying target
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
  array[7] real mu_pr;
  array[7] real<lower=0> sigma;

  array[N] real con_pr;
  array[N] real gain_pr;
  array[N] real loss_pr;
  array[N] real decay_pr;
  array[N] real explore_alpha_pr;
  array[N] real explore_bonus_pr;
  array[N] real update_pr;
}

transformed parameters {
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=1> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> decay;
  array[N] real<lower=0, upper=1> explore_alpha;
  array[N] real<lower=-10, upper=10> explore_bonus;
  array[N] real<lower=0, upper=1> update;
  
  con           = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(con_pr)) * 5);
  gain          = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(gain_pr)));
  loss          = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(loss_pr)) * 10);
  decay         = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(decay_pr)));
  explore_alpha = to_array_1d(inv_logit(mu_pr[5] + sigma[5] .* to_vector(explore_alpha_pr)));
  explore_bonus = to_array_1d(-10 + inv_logit(mu_pr[6] + sigma[6] .* to_vector(explore_bonus_pr)) * 20);
  update        = to_array_1d(inv_logit(mu_pr[7] + sigma[7] .* to_vector(update_pr)));
}

model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  con_pr           ~ normal(0, 1);
  gain_pr          ~ normal(0, 1);
  loss_pr          ~ normal(0, 1);
  decay_pr         ~ normal(0, 1);
  explore_alpha_pr ~ normal(0, 1);
  explore_bonus_pr ~ normal(0, 1);
  update_pr        ~ normal(0, 1);

  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       con, gain, loss, update, decay, explore_alpha, explore_bonus);
}

generated quantities {
  real<lower=0, upper=5> mu_con;
  real<lower=0, upper=1> mu_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0, upper=1> mu_decay;
  real<lower=0, upper=1> mu_explore_alpha;
  real<lower=-10, upper=10> mu_explore_bonus;
  real<lower=0, upper=1> mu_update;
  
  mu_con           = inv_logit(mu_pr[1]) * 5;
  mu_gain          = inv_logit(mu_pr[2]);
  mu_loss          = inv_logit(mu_pr[3]) * 10;
  mu_decay         = inv_logit(mu_pr[4]);
  mu_explore_alpha = inv_logit(mu_pr[5]);
  mu_explore_bonus = -10 + inv_logit(mu_pr[6]) * 20;
  mu_update        = inv_logit(mu_pr[7]);
}
