functions {
  vector igt_model_lp(
      array[] int choice, array[] int shown, array[] real outcome,
      array[] real RT, vector ev, int Tsub, 
      real sensitivity, real update, real wgt_pun,
      real wgt_rew, real boundary1, real boundary, real tau1, real tau, real beta
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
      curUtil = (outcome[t] > 0 ? wgt_rew : wgt_pun) * outcome[t] * choice[t];
      local_ev[curDeck] += (curUtil - local_ev[curDeck]) * update * choice[t];

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

    target += wiener_lpdf(RT[play_indices[:play_count]] | boundaries[play_indices[:play_count]], tau, beta, drift_rates[play_indices[:play_count]]);
    target += wiener_lpdf(RT[pass_indices[:pass_count]] | boundaries[pass_indices[:pass_count]], tau, 1-beta, -drift_rates[pass_indices[:pass_count]]);

    return local_ev;
  }
}

data {
  int<lower=1> sid;
  int<lower=1> T;
  real<lower=0> minRT;
  real RTbound;
  array[T] real<lower=0> RT;
  array[T] int<lower=0, upper=1> choice;
  array[T] int<lower=0, upper=4> shown;
  array[T] real outcome;
  vector[9] pr_mu;
  vector[9] pr_sigma;
}

parameters {
  real<lower=-5, upper=5> boundary1_pr;
  real<lower=-5, upper=5> boundary_pr;
  real<lower=-3, upper=3> tau1_pr;
  real<lower=-3, upper=3> tau_pr;
  real<lower=-3, upper=3> beta_pr;
  real<lower=-3, upper=3> drift_con_pr;
  real<lower=-3, upper=3> wgt_pun_pr;
  real<lower=-3, upper=3> wgt_rew_pr;
  real<lower=-3, upper=3> update_pr;
}

transformed parameters {
  real<lower=0, upper=6> boundary1;
  real<lower=0, upper=6> boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0, upper=1> beta;
  real<lower=0, upper=3> drift_con;
  real<lower=0, upper=1> wgt_pun;
  real<lower=0, upper=1> wgt_rew;
  real<lower=0, upper=1> update;

  boundary1 = inv_logit(boundary1_pr) * 5 + 0.01;
  boundary  = inv_logit(boundary_pr) * 5 + 0.01;
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound;
  beta      = inv_logit(beta_pr);
  drift_con = inv_logit(drift_con_pr) * 3;
  wgt_pun   = inv_logit(wgt_pun_pr);
  wgt_rew   = inv_logit(wgt_rew_pr);
  update    = inv_logit(update_pr);
}

model {
  boundary1_pr  ~ normal(pr_mu[1], pr_sigma[1]);
  boundary_pr  ~ normal(pr_mu[2], pr_sigma[2]);
  tau1_pr       ~ normal(pr_mu[3], pr_sigma[3]);
  tau_pr        ~ normal(pr_mu[4], pr_sigma[4]);
  beta_pr      ~ normal(pr_mu[5], pr_sigma[5]);
  drift_con_pr ~ normal(pr_mu[6], pr_sigma[6]);
  update_pr    ~ normal(pr_mu[7], pr_sigma[7]);
  wgt_pun_pr   ~ normal(pr_mu[8], pr_sigma[8]);
  wgt_rew_pr   ~ normal(pr_mu[9], pr_sigma[9]);

  vector[4] ev = rep_vector(0., 4);
  real sensitivity = pow(3, drift_con) - 1;

  {
    ev = igt_model_lp(choice, shown, outcome, RT, ev, T,
                      sensitivity, update, wgt_pun,
                      wgt_rew, boundary1, boundary, tau1, tau, beta);
  }
}
