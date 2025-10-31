// Optimized Hierarchical EV Model for IGT_MOD
functions {
  // Returns log-likelihood instead of modifying target
  real igt_subject(
      array[] int choice, array[] int shown, array[] real outcome,
      vector ev, int Tsub, real sensitivity,
      real update, real wgt_pun, real wgt_rew
      ) {
    real log_lik = 0.0;
    real curUtil;
    int curDeck;
    vector[Tsub] Info;
    vector[4] local_ev = ev;

    for (t in 1:Tsub) {
      curDeck = shown[t];
      Info[t] = sensitivity * local_ev[curDeck];
      curUtil = ((outcome[t] > 0 ? wgt_rew : wgt_pun)) * outcome[t] * choice[t];
      local_ev[curDeck] += (curUtil - local_ev[curDeck]) * update * choice[t];
    }
    
    log_lik += bernoulli_logit_lpmf(choice | Info);
    return log_lik;
  }
  
  // Parallelization wrapper
  real partial_sum_func(array[] int slice_n, int start, int end,
                        array[,] int choice, array[,] int shown, array[,] real outcome,
                        array[] int Tsubj,
                        array[] real con, array[] real update,
                        array[] real wgt_pun, array[] real wgt_rew) {
    real log_lik = 0.0;
    
    for (n in start:end) {
      vector[4] ev = rep_vector(0., 4);
      real sensitivity = expm1(log(3) * con[n]);
      
      log_lik += igt_subject(choice[n, 1:Tsubj[n]], shown[n, 1:Tsubj[n]], 
                             outcome[n, 1:Tsubj[n]], ev, Tsubj[n], sensitivity,
                             update[n], wgt_pun[n], wgt_rew[n]);
    }
    return log_lik;
  }
}

data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1> Tsubj;
  array[N, T] int<lower=0, upper=1> choice;
  array[N, T] int<lower=0, upper=4> shown;
  array[N, T] real outcome;
}

transformed data {
  array[N] int subject_indices;
  for (i in 1:N) {
    subject_indices[i] = i;
  }
}

parameters {
  array[4] real mu_pr;
  array[4] real<lower=0> sigma;

  array[N] real con_pr;
  array[N] real wgt_pun_pr;
  array[N] real wgt_rew_pr;
  array[N] real update_pr;
}

transformed parameters {
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=1> wgt_pun;
  array[N] real<lower=0, upper=1> wgt_rew;
  array[N] real<lower=0, upper=1> update;

  con     = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(con_pr)) * 5);
  wgt_pun = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(wgt_pun_pr)));
  wgt_rew = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(wgt_rew_pr)));
  update  = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(update_pr)));
}

model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  con_pr     ~ normal(0, 1);
  wgt_pun_pr ~ normal(0, 1);
  wgt_rew_pr ~ normal(0, 1);
  update_pr  ~ normal(0, 1);

  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, shown, outcome, Tsubj,
                       con, update, wgt_pun, wgt_rew);
}

generated quantities {
  real<lower=0, upper=5> mu_con;
  real<lower=0, upper=1> mu_wgt_pun;
  real<lower=0, upper=1> mu_wgt_rew;
  real<lower=0, upper=1> mu_update;
  
  mu_con     = inv_logit(mu_pr[1]) * 5;
  mu_wgt_pun = inv_logit(mu_pr[2]);
  mu_wgt_rew = inv_logit(mu_pr[3]);
  mu_update  = inv_logit(mu_pr[4]);
}
