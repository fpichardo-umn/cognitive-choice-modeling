// Optimized Hierarchical PVL-Both-DDM b1p2 Model for IGT_MOD
functions {
  vector igt_model_lp(
      array[] int choice, array[] int shown, array[] real outcome,
      array[] real RT, vector ev, int Tsub,
      real boundary1, real boundary, real tau1, real tau, real beta,
      real sensitivity, real gain, real loss, real update, real decay
      ) {
    real curUtil;
    int curDeck;
    vector[Tsub] drift_rates;
    vector[Tsub] boundaries;
    vector[Tsub] nondt;
    vector[4] local_ev = ev;

    array[Tsub] int play_indices;
    array[Tsub] int pass_indices;
    int play_count = 0;
    int pass_count = 0;

    for (t in 1:Tsub) {
      if (t <= 20) {
        boundaries[t] = boundary1;
        nondt[t] = tau1;
      } else {
        boundaries[t] = boundary;
        nondt[t] = tau;
      }
      
      curDeck = shown[t];
      drift_rates[t] = local_ev[curDeck] * sensitivity;

      if (outcome[t] >= 0) {
        curUtil = pow(outcome[t], gain);
      } else {
        curUtil = -loss * pow(abs(outcome[t]), gain);
      }

      local_ev = local_ev * (1 - decay);
      if (choice[t] == 1) {
        local_ev[curDeck] += update * (curUtil - local_ev[curDeck]);
      }

      if (RT[t] != 999) {
        if (choice[t] == 1) {
          play_count += 1;
          play_indices[play_count] = t;
        } else {
          pass_count += 1;
          pass_indices[pass_count] = t;
        }
      }
    }

    target += wiener_lpdf(RT[play_indices[:play_count]] | boundaries[play_indices[:play_count]], nondt[play_indices[:play_count]], beta, drift_rates[play_indices[:play_count]]);
    target += wiener_lpdf(RT[pass_indices[:pass_count]] | boundaries[pass_indices[:pass_count]], nondt[pass_indices[:pass_count]], 1-beta, -drift_rates[pass_indices[:pass_count]]);

    return local_ev;
  }
}

data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1> Tsubj;
  array[N] real<lower=0> minRT;
  real RTbound;
  array[N, T] real<lower=0> RT;
  array[N, T] int<lower=0, upper=1> choice;
  array[N, T] int<lower=0, upper=4> shown;
  array[N, T] real outcome;
}

parameters {
  array[10] real mu_pr;
  array[10] real<lower=0> sigma;

  array[N] real<lower=-5, upper=5> boundary1_pr;
  array[N] real<lower=-5, upper=5> boundary_pr;
  array[N] real tau1_pr;
  array[N] real tau_pr;
  array[N] real beta_pr;
  array[N] real drift_con_pr;
  array[N] real gain_pr;
  array[N] real loss_pr;
  array[N] real update_pr;
  array[N] real decay_pr;
}

transformed parameters {
  array[N] real<lower=0, upper=6> boundary1;
  array[N] real<lower=0, upper=6> boundary;
  array[N] real<lower=RTbound, upper=max(minRT)> tau1;
  array[N] real<lower=RTbound, upper=max(minRT)> tau;
  array[N] real<lower=0, upper=1> beta;
  array[N] real<lower=0, upper=5> drift_con;
  array[N] real<lower=0, upper=2> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> update;
  array[N] real<lower=0, upper=1> decay;

  boundary1 = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(boundary1_pr)) * 5 + 0.01);
  boundary  = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(boundary_pr)) * 5 + 0.01);
  tau1      = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(tau1_pr)) .* (to_vector(minRT) - RTbound - 1e-6) * 0.99 + RTbound);
  tau       = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(tau_pr)) .* (to_vector(minRT) - RTbound - 1e-6) * 0.99 + RTbound);
  beta      = to_array_1d(inv_logit(mu_pr[5] + sigma[5] .* to_vector(beta_pr)));
  drift_con = to_array_1d(inv_logit(mu_pr[6] + sigma[6] .* to_vector(drift_con_pr)) * 5);
  gain      = to_array_1d(inv_logit(mu_pr[7] + sigma[7] .* to_vector(gain_pr)) * 2);
  loss      = to_array_1d(inv_logit(mu_pr[8] + sigma[8] .* to_vector(loss_pr)) * 10);
  update    = to_array_1d(inv_logit(mu_pr[9] + sigma[9] .* to_vector(update_pr)));
  decay     = to_array_1d(inv_logit(mu_pr[10] + sigma[10] .* to_vector(decay_pr)));
}

model {
  mu_pr ~ normal(0, 1);
  sigma ~ normal(0, 1);

  boundary1_pr ~ normal(0, 1);
  boundary_pr  ~ normal(0, 1);
  tau1_pr      ~ normal(0, 1);
  tau_pr       ~ normal(0, 1);
  beta_pr      ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  gain_pr      ~ normal(0, 1);
  loss_pr      ~ normal(0, 1);
  update_pr    ~ normal(0, 1);
  decay_pr     ~ normal(0, 1);

  for (n in 1:N) {
    vector[4] ev = rep_vector(0., 4);
    real sensitivity = pow(3, drift_con[n]) - 1;
    
    ev = igt_model_lp(choice[n, 1:Tsubj[n]], shown[n, 1:Tsubj[n]], outcome[n, 1:Tsubj[n]],
                      RT[n, 1:Tsubj[n]], ev, Tsubj[n],
                      boundary1[n], boundary[n], tau1[n], tau[n], beta[n],
                      sensitivity, gain[n], loss[n], update[n], decay[n]);
  }
}

generated quantities {
  real<lower=0, upper=6> mu_boundary1;
  real<lower=0, upper=6> mu_boundary;
  real<lower=0> mu_tau1;
  real<lower=0> mu_tau;
  real<lower=0, upper=1> mu_beta;
  real<lower=0, upper=5> mu_drift_con;
  real<lower=0, upper=2> mu_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0, upper=1> mu_update;
  real<lower=0, upper=1> mu_decay;
  
  mu_boundary1 = inv_logit(mu_pr[1]) * 5 + 0.01;
  mu_boundary  = inv_logit(mu_pr[2]) * 5 + 0.01;
  mu_tau1      = inv_logit(mu_pr[3]) * (mean(minRT) - RTbound) * 0.99 + RTbound;
  mu_tau       = inv_logit(mu_pr[4]) * (mean(minRT) - RTbound) * 0.99 + RTbound;
  mu_beta      = inv_logit(mu_pr[5]);
  mu_drift_con = inv_logit(mu_pr[6]) * 5;
  mu_gain      = inv_logit(mu_pr[7]) * 2;
  mu_loss      = inv_logit(mu_pr[8]) * 10;
  mu_update    = inv_logit(mu_pr[9]);
  mu_decay     = inv_logit(mu_pr[10]);
}
