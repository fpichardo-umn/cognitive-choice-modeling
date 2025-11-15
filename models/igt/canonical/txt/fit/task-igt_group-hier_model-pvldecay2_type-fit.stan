// Optimized Hierarchical PVL-Decay Model for the Iowa Gambling Task
functions {
  // Returns log-likelihood instead of modifying target
  real igt_subject(
        array[] int choice, array[] real outcome,
        vector ev, int Tsub,
        real sensitivity, real gain, real loss, real decay
        ) {
    real log_lik = 0.0;
    real curUtil;
    vector[4] local_ev = ev;
    
    for (t in 1:Tsub) {
      log_lik += categorical_logit_lpmf(choice[t] | sensitivity * local_ev);

      if (outcome[t] >= 0) {
        curUtil = pow(outcome[t], gain);
      } else {
        curUtil = -1 * loss * pow(-1 * outcome[t], gain);
      }

      local_ev = local_ev * (1 - decay);
      local_ev[choice[t]] += curUtil;
    }
    
    return log_lik;
  }
  
  // Parallelization wrapper
  real partial_sum_func(array[] int slice_n, int start, int end,
                        array[,] int choice, array[,] real wins, array[,] real losses, 
                        array[] int Tsubj,
                        array[] real con, array[] real gain, 
                        array[] real loss, array[] real decay) {
    real log_lik = 0.0;
    
    for (n in start:end) {
      vector[4] ev = rep_vector(0., 4);
      real sensitivity = expm1(log(3) * con[n]);
      array[Tsubj[n]] real outcome = wins[n, 1:Tsubj[n]] - losses[n, 1:Tsubj[n]];
      
      log_lik += igt_subject(choice[n, 1:Tsubj[n]],
                             outcome, 
                             ev, Tsubj[n], sensitivity, gain[n], 
                             loss[n], decay[n]);
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
  array[4] real mu_pr;
  array[4] real<lower=0> sigma;

  array[N] real con_pr;
  array[N] real gain_pr;
  array[N] real loss_pr;
  array[N] real decay_pr;
}

transformed parameters {
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=2> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> decay;
  
  con   = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(con_pr)) * 5);
  gain  = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(gain_pr)) * 2);
  loss  = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(loss_pr)) * 10);
  decay = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(decay_pr)));
}

model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  con_pr   ~ normal(0, 1);
  gain_pr  ~ normal(0, 1);
  loss_pr  ~ normal(0, 1);
  decay_pr ~ normal(0, 1);

  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       con, gain, loss, decay);
}

generated quantities {
  real<lower=0, upper=5> mu_con;
  real<lower=0, upper=2> mu_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0, upper=1> mu_decay;
  
  mu_con   = inv_logit(mu_pr[1]) * 5;
  mu_gain  = inv_logit(mu_pr[2]) * 2;
  mu_loss  = inv_logit(mu_pr[3]) * 10;
  mu_decay = inv_logit(mu_pr[4]);
}
