// Optimized Hierarchical LagVSEBoth Model for the Iowa Gambling Task
functions {
  // Returns log-likelihood instead of modifying target
  real igt_subject(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev_exploit, vector choice_lag, int Tsub,
        real sensitivity, real gain, real loss,
        real update, real decay, real phi
        ) {
    real log_lik = 0.0;
    real curUtil;
    vector[4] local_ev_exploit = ev_exploit;
    vector[4] local_choice_lag = choice_lag;
    vector[4] combined_value;
    
    for (t in 1:Tsub) {
      local_choice_lag += 1;
      combined_value = local_ev_exploit + phi * local_choice_lag;
      log_lik += categorical_logit_lpmf(choice[t] | sensitivity * combined_value);
      
      real win_component = (wins[t] == 0) ? 0.0 : exp(gain * log(wins[t]));
real loss_component = (losses[t] == 0) ? 0.0 : exp(gain * log(losses[t]));
curUtil = win_component - loss * loss_component;
      local_ev_exploit = local_ev_exploit * (1 - decay);
      local_ev_exploit[choice[t]] += curUtil + update * (curUtil - local_ev_exploit[choice[t]]);
      local_choice_lag[choice[t]] = 0;
    }
    
    return log_lik;
  }
  
  // Parallelization wrapper
  real partial_sum_func(array[] int slice_n, int start, int end,
                        array[,] int choice, array[,] real wins, array[,] real losses, 
                        array[] int Tsubj,
                        array[] real con, array[] real gain, array[] real loss,
                        array[] real update, array[] real decay, array[] real phi) {
    real log_lik = 0.0;
    
    for (n in start:end) {
      vector[4] ev_exploit = rep_vector(0., 4);
      vector[4] ev_explore = rep_vector(0., 4);
      real sensitivity = expm1(log(3) * con[n]);
      
      log_lik += igt_subject(choice[n, 1:Tsubj[n]], 
                             wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]],
                             ev_exploit, ev_explore, Tsubj[n], 
                             sensitivity, gain[n], loss[n],
                             update[n], decay[n], phi[n]);
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
  array[6] real mu_pr;
  array[6] real<lower=0> sigma;

  array[N] real con_pr;
  array[N] real gain_pr;
  array[N] real loss_pr;
  array[N] real update_pr;
  array[N] real decay_pr;
  array[N] real phi_pr;
}

transformed parameters {
  array[N] real<lower=0, upper=3> con;
  array[N] real<lower=0, upper=1> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> update;
  array[N] real<lower=0, upper=1> decay;
  array[N] real<lower=-1, upper=1> phi;
  
  con    = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(con_pr)) * 3);
  gain   = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(gain_pr)));
  loss   = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(loss_pr)) * 10);
  update = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(update_pr)));
  decay  = to_array_1d(inv_logit(mu_pr[5] + sigma[5] .* to_vector(decay_pr)));
  phi    = to_array_1d(-1 + inv_logit(mu_pr[6] + sigma[6] .* to_vector(phi_pr)) * 2);
}

model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  con_pr    ~ normal(0, 1);
  gain_pr   ~ normal(0, 1);
  loss_pr   ~ normal(0, 1);
  update_pr ~ normal(0, 1);
  decay_pr  ~ normal(0, 1);
  phi_pr    ~ normal(0, 1);

  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       con, gain, loss, update, decay, phi);
}

generated quantities {
  real<lower=0, upper=5> mu_con;
  real<lower=0, upper=1> mu_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0, upper=1> mu_update;
  real<lower=0, upper=1> mu_decay;
  real<lower=-1, upper=1> mu_phi;
  
  mu_con    = inv_logit(mu_pr[1]) * 3;
  mu_gain   = inv_logit(mu_pr[2]);
  mu_loss   = inv_logit(mu_pr[3]) * 10;
  mu_update = inv_logit(mu_pr[4]);
  mu_decay  = inv_logit(mu_pr[5]);
  mu_phi    = -1 + inv_logit(mu_pr[6]) * 2;
}
