// Optimized Hybrid Hierarchical VPPBoth Model for the IGT
functions {
  real igt_subject(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, vector pers, int Tsub,
        real sensitivity, real update, real gain, real loss, 
        real epP, real epN, real K, real w, real decay
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
      
      local_ev = local_ev * (1 - decay);
      local_ev[choice[t]] += curUtil + update * (curUtil - local_ev[choice[t]]);
    }
    return log_lik;
  }
  
  real partial_sum_func(array[] int slice_n, int start, int end,
                        array[,] int choice, array[,] real wins, array[,] real losses, 
                        array[] int Tsubj,
                        array[] real con, array[] real update, array[] real gain, 
                        array[] real loss, array[] real epP, array[] real epN,
                        array[] real K, array[] real w, array[] real decay) {
    real log_lik = 0.0;
    for (n in start:end) {
      vector[4] ev = rep_vector(0., 4);
      vector[4] pers = rep_vector(0., 4);
      real sensitivity = expm1(log(3) * con[n]);
      
      log_lik += igt_subject(choice[n, 1:Tsubj[n]],
                             wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], 
                             ev, pers, Tsubj[n], 
                             sensitivity, update[n], gain[n], loss[n], 
                             epP[n], epN[n], K[n], w[n], decay[n]);
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
  // NCP Hyperparameters (epP, epN, K)
  vector[3] mu_pr_rest;
  vector<lower=0>[3] sigma_rest;

  // CP Hyperparameters (Means and SDs for the sluggards)
  real<lower=0, upper=5> mu_con;
  real<lower=0> sigma_con;
  real<lower=0, upper=1> mu_update;
  real<lower=0> sigma_update;
  real<lower=0, upper=1> mu_gain;
  real<lower=0> sigma_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0> sigma_loss;
  real<lower=0, upper=1> mu_w;
  real<lower=0> sigma_w;
  real<lower=0, upper=1> mu_decay;
  real<lower=0> sigma_decay;

  // Individual parameters
  vector[N] epP_pr;
  vector[N] epN_pr;
  vector[N] K_pr;
  
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=1> update;
  array[N] real<lower=0, upper=1> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> w;
  array[N] real<lower=0, upper=1> decay;
}

transformed parameters {
  array[N] real<lower=-5, upper=5> epP;
  array[N] real<lower=-5, upper=5> epN;
  array[N] real<lower=0, upper=1> K;
  
  // Non-Centered reparameterization for persistence parameters
  epP = to_array_1d(-5 + inv_logit(mu_pr_rest[1] + sigma_rest[1] * epP_pr) * 10);
  epN = to_array_1d(-5 + inv_logit(mu_pr_rest[2] + sigma_rest[2] * epN_pr) * 10);
  K   = to_array_1d(inv_logit(mu_pr_rest[3] + sigma_rest[3] * K_pr));
}

model {
  // 1. Tightened Global Priors (Normal instead of Student-T) 
  mu_pr_rest ~ normal(0, 1);
  sigma_rest ~ normal(0, 0.5); 
  
  mu_con    ~ normal(2.5, 1);
  mu_update ~ beta(2, 2);
  mu_gain   ~ beta(2, 2);
  mu_loss   ~ normal(5, 2);
  mu_w      ~ beta(2, 2);
  mu_decay  ~ beta(2, 2);

  sigma_con    ~ normal(0, 0.5);
  sigma_update ~ normal(0, 0.5);
  sigma_gain   ~ normal(0, 0.5);
  sigma_loss   ~ normal(0, 1);
  sigma_w      ~ normal(0, 0.5);
  sigma_decay  ~ normal(0, 0.5);

  // 2. Hierarchy [cite: 138-141]
  epP_pr ~ normal(0, 1);
  epN_pr ~ normal(0, 1);
  K_pr   ~ normal(0, 1);
  
  con    ~ normal(mu_con, sigma_con);
  update ~ normal(mu_update, sigma_update);
  gain   ~ normal(mu_gain, sigma_gain);
  loss   ~ normal(mu_loss, sigma_loss);
  w      ~ normal(mu_w, sigma_w);
  decay  ~ normal(mu_decay, sigma_decay);

  // 3. Likelihood [cite: 142]
  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       con, update, gain, loss, epP, epN, K, w, decay);
}

generated quantities {
  real mu_epP = -5 + inv_logit(mu_pr_rest[1]) * 10;
  real mu_epN = -5 + inv_logit(mu_pr_rest[2]) * 10;
  real mu_K   = inv_logit(mu_pr_rest[3]);
  // Others are in parameters
}
