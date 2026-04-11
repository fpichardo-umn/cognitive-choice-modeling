// Refined Hybrid Hierarchical PVL-Delta4 Model
functions {
  real igt_subject(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, int Tsub,
        real sensitivity, real update, real gain, real loss
        ) {
    real log_lik = 0.0;
    real curUtil;
    vector[4] local_ev = ev;

    for (t in 1:Tsub) {
      log_lik += categorical_logit_lpmf(choice[t] | sensitivity * local_ev);

      real win_component = (wins[t] == 0) ? 0.0 : exp(gain * log(wins[t]));
      real loss_component = (losses[t] == 0) ? 0.0 : exp(gain * log(losses[t]));
      curUtil = win_component - loss * loss_component;
      
      local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
    }
    return log_lik;
  }

  real partial_sum_func(array[] int slice_n, int start, int end,
                        array[,] int choice, array[,] real wins, array[,] real losses,
                        array[] int Tsubj,
                        array[] real con, array[] real update, array[] real gain,
                        array[] real loss) {
    real log_lik = 0.0;
    for (n in start:end) {
      vector[4] ev = rep_vector(0., 4);
      real sensitivity = expm1(log(3) * con[n]);

      log_lik += igt_subject(choice[n, 1:Tsubj[n]],
                             wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]],
                             ev, Tsubj[n], sensitivity, update[n],
                             gain[n], loss[n]);
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
  // NCP Hyperparameters (update, con, gain)
  vector[3] mu_pr_rest;
  vector<lower=0>[3] sigma_rest;

  // CP Hyperparameters (loss)
  real<lower=0, upper=10> mu_loss;
  real<lower=0> sigma_loss;

  // Individual parameters
  vector[N] update_pr;
  vector[N] con_pr;
  vector[N] gain_pr;
  
  array[N] real<lower=0, upper=10> loss; // Centered
}

transformed parameters {
  array[N] real<lower=0, upper=1> update;
  array[N] real<lower=0, upper=5> con;
  array[N] real<lower=0, upper=1> gain;

  // NCP reparameterization
  update = to_array_1d(inv_logit(mu_pr_rest[1] + sigma_rest[1] * update_pr));
  con    = to_array_1d(inv_logit(mu_pr_rest[2] + sigma_rest[2] * con_pr) * 5);
  gain   = to_array_1d(inv_logit(mu_pr_rest[3] + sigma_rest[3] * gain_pr));
}

model {
  // Tighter Priors
  mu_pr_rest ~ normal(0, 1);
  sigma_rest ~ normal(0, 0.5); 
  
  mu_loss ~ normal(5, 2);
  sigma_loss ~ normal(0, 1);

  update_pr ~ normal(0, 1);
  con_pr    ~ normal(0, 1);
  gain_pr   ~ normal(0, 1);

  // CP Likelihood for Loss
  loss ~ normal(mu_loss, sigma_loss);

  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       con, update, gain, loss);
}

generated quantities {
  real mu_update = inv_logit(mu_pr_rest[1]);
  real mu_con    = inv_logit(mu_pr_rest[2]) * 5;
  real mu_gain   = inv_logit(mu_pr_rest[3]);
  // mu_loss is in parameters
}
