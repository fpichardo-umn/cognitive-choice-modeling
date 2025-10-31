// Optimized Hierarchical PVL-Decay Model for IGT_MOD
functions {
  // Returns log-likelihood instead of modifying target
  real igt_subject(
      array[] int choice, array[] int shown, array[] real outcome,
      vector ev, int Tsub, real sensitivity,
      real gain, real loss, real retend
      ) {
    real log_lik = 0.0;
    real curUtil;
    int curDeck;
    vector[Tsub] Info;
    vector[4] local_ev = ev;

    for (t in 1:Tsub) {
      curDeck = shown[t];
      Info[t] = sensitivity * local_ev[curDeck];

      if (outcome[t] > 0) {
  curUtil = exp(gain * log(outcome[t]));
} else if (outcome[t] < 0) {
  curUtil = -loss * exp(gain * log(-outcome[t]));
} else {
  curUtil = 0;
}

      local_ev *= retend;
      local_ev[curDeck] += curUtil * choice[t];
    }
    
    log_lik += bernoulli_logit_lpmf(choice | Info);
    return log_lik;
  }
  
  // Parallelization wrapper
  real partial_sum_func(array[] int slice_n, int start, int end,
                        array[,] int choice, array[,] int shown, array[,] real outcome,
                        array[] int Tsubj,
                        array[] real con, array[] real gain,
                        array[] real loss, array[] real retend) {
    real log_lik = 0.0;
    
    for (n in start:end) {
      vector[4] ev = rep_vector(0., 4);
      real sensitivity = expm1(log(3) * con[n]);
      
      log_lik += igt_subject(choice[n, 1:Tsubj[n]], shown[n, 1:Tsubj[n]], 
                             outcome[n, 1:Tsubj[n]], ev, Tsubj[n], sensitivity,
                             gain[n], loss[n], retend[n]);
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
  array[N] real gain_pr;
  array[N] real loss_pr;
  array[N] real retend_pr;
}

transformed parameters {
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=2> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> retend;

  con    = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(con_pr)) * 5);
  gain   = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(gain_pr)) * 2);
  loss   = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(loss_pr)) * 10);
  retend = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(retend_pr)));
}

model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  con_pr    ~ normal(0, 1);
  gain_pr   ~ normal(0, 1);
  loss_pr   ~ normal(0, 1);
  retend_pr ~ normal(0, 1);

  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, shown, outcome, Tsubj,
                       con, gain, loss, retend);
}

generated quantities {
  real<lower=0, upper=5> mu_con;
  real<lower=0, upper=2> mu_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0, upper=1> mu_retend;
  
  mu_con    = inv_logit(mu_pr[1]) * 5;
  mu_gain   = inv_logit(mu_pr[2]) * 2;
  mu_loss   = inv_logit(mu_pr[3]) * 10;
  mu_retend = inv_logit(mu_pr[4]);
}
