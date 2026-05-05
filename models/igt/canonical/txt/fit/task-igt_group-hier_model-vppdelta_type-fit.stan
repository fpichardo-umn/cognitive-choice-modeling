// Optimized Hierarchical VPP-Delta Model for the Iowa Gambling Task
functions {
  // Returns log-likelihood instead of modifying target
  real igt_subject(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, vector pers, int Tsub,
        real sensitivity, real update, real gain, real loss, 
        real epP, real epN, real K, real w
        ) {
    real log_lik = 0.0;
    real curUtil;
    vector[4] local_ev = ev;
    vector[4] local_pers = pers;
    vector[4] V;

    for (t in 1:Tsub) {
      V = w * local_ev + (1-w) * local_pers;
      log_lik += categorical_logit_lpmf(choice[t] | sensitivity * V);
      
      local_pers = local_pers * K;
      real win_component = (wins[t] == 0) ? 0.0 : exp(gain * log(wins[t]));
real loss_component = (losses[t] == 0) ? 0.0 : exp(gain * log(losses[t]));
curUtil = win_component - loss * loss_component;
      
      if (wins[t] >= losses[t]) {
        local_pers[choice[t]] += epP;
      } else {
        local_pers[choice[t]] += epN;
      }
      
      local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
    }
    
    return log_lik;
  }
  
  // Parallelization wrapper
  real partial_sum_func(array[] int slice_n, int start, int end,
                        array[,] int choice, array[,] real wins, array[,] real losses, 
                        array[] int Tsubj,
                        array[] real con, array[] real update, array[] real gain, 
                        array[] real loss, array[] real epP, array[] real epN,
                        array[] real K, array[] real w) {
    real log_lik = 0.0;
    
    for (n in start:end) {
      vector[4] ev = rep_vector(0., 4);
      vector[4] pers = rep_vector(0., 4);
      real sensitivity = expm1(log(3) * con[n]);
      
      log_lik += igt_subject(choice[n, 1:Tsubj[n]],
                             wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], 
                             ev, pers, Tsubj[n], sensitivity, update[n], 
                             gain[n], loss[n], epP[n], epN[n], K[n], w[n]);
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
  array[8] real mu_pr;
  array[8] real<lower=0> sigma;

  array[N] real con_pr;
  array[N] real update_pr;
  array[N] real gain_pr;
  array[N] real loss_pr;
  array[N] real epP_pr;
  array[N] real epN_pr;
  array[N] real K_pr;
  array[N] real w_pr;
}

transformed parameters {
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=1> update;
  array[N] real<lower=0, upper=1> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real epP;
  array[N] real epN;
  array[N] real<lower=0, upper=1> K;
  array[N] real<lower=0, upper=1> w;
  
  con    = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(con_pr)) * 5);
  update = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(update_pr)));
  gain   = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(gain_pr)));
  loss   = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(loss_pr)) * 10);
  epP    = to_array_1d(mu_pr[5] + sigma[5] .* to_vector(epP_pr));
  epN    = to_array_1d(mu_pr[6] + sigma[6] .* to_vector(epN_pr));
  K      = to_array_1d(inv_logit(mu_pr[7] + sigma[7] .* to_vector(K_pr)));
  w      = to_array_1d(inv_logit(mu_pr[8] + sigma[8] .* to_vector(w_pr)));
}

model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  con_pr    ~ normal(0, 1);
  update_pr ~ normal(0, 1);
  gain_pr   ~ normal(0, 1);
  loss_pr   ~ normal(0, 1);
  epP_pr    ~ normal(0, 1);
  epN_pr    ~ normal(0, 1);
  K_pr      ~ normal(0, 1);
  w_pr      ~ normal(0, 1);

  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       con, update, gain, loss, epP, epN, K, w);
}

generated quantities {
  real<lower=0, upper=5> mu_con;
  real<lower=0, upper=1> mu_update;
  real<lower=0, upper=1> mu_gain;
  real<lower=0, upper=10> mu_loss;
  real mu_epP;
  real mu_epN;
  real<lower=0, upper=1> mu_K;
  real<lower=0, upper=1> mu_w;
  
  mu_con    = inv_logit(mu_pr[1]) * 5;
  mu_update = inv_logit(mu_pr[2]);
  mu_gain   = inv_logit(mu_pr[3]);
  mu_loss   = inv_logit(mu_pr[4]) * 10;
  mu_epP    = mu_pr[5];
  mu_epN    = mu_pr[6];
  mu_K      = inv_logit(mu_pr[7]);
  mu_w      = inv_logit(mu_pr[8]);
}
