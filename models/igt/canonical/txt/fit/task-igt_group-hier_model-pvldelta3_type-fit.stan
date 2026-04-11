// Optimized Hybrid Hierarchical PVL-Delta Model for the IGT
functions {
  real igt_subject(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, int Tsub,
        real sensitivity, real gain, real loss, real update
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
                        array[] real con, array[] real gain,
                        array[] real loss, array[] real update) {
    real log_lik = 0.0;
    for (n in start:end) {
      vector[4] ev = rep_vector(0., 4);
      real sensitivity = expm1(log(3) * con[n]);

      log_lik += igt_subject(choice[n, 1:Tsubj[n]],
                             wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]],
                             ev, Tsubj[n], sensitivity, gain[n],
                             loss[n], update[n]);
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
  // CP Hyperparameters (Means and SDs)
  real<lower=0, upper=5> mu_con;
  real<lower=0> sigma_con;
  real<lower=0, upper=2> mu_gain;
  real<lower=0> sigma_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0> sigma_loss;
  real<lower=0, upper=1> mu_update;
  real<lower=0> sigma_update;

  // Individual parameters
  array[N] real<lower=0, upper=5> con;    // Centered
  array[N] real<lower=0, upper=2> gain;   // Centered
  array[N] real<lower=0, upper=10> loss;  // Centered
  array[N] real<lower=0, upper=1> update; // Centered
}

model {
  // Weakly informative priors for Centered means
  mu_con  ~ normal(2.5, 1);   // Midpoint of 0-5
  sigma_con ~ normal(0, 0.5);
  mu_gain ~ normal(1.0, 0.5); // Midpoint of 0-2
  sigma_gain ~ normal(0, 0.5);
  mu_loss ~ normal(5.0, 2);   // Midpoint of 0-10
  sigma_loss ~ normal(0, 1);
  mu_update ~ normal(0.5, 2);   // Midpoint of 0-1
  sigma_update ~ normal(0, 1);

  // 2. Hierarchy
  con    ~ normal(mu_con, sigma_con);
  gain   ~ normal(mu_gain, sigma_gain);
  loss   ~ normal(mu_loss, sigma_loss);
  update ~ normal(mu_update, sigma_update);

  // 3. Likelihood
  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       con, gain, loss, update);
}
