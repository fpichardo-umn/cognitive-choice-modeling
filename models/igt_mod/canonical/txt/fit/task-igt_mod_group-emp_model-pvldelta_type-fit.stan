functions {
  vector igt_model_lp(
      array[] int choice, array[] int shown, array[] real outcome,
      vector ev, int Tsub, real sensitivity,
      real gain, real loss, real update
      ) {
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

      local_ev[curDeck] += (curUtil - local_ev[curDeck]) * update * choice[t];
    }
    
    target += bernoulli_logit_lpmf(choice | Info);
    return local_ev;
  }
}

data {
  int<lower=1> sid;
  int<lower=1> T;
  array[T] int<lower=0, upper=1> choice;
  array[T] int<lower=0, upper=4> shown;
  array[T] real outcome;
  vector[4] pr_mu;
  vector[4] pr_sigma;
}

parameters {
  real con_pr;
  real gain_pr;
  real loss_pr;
  real update_pr;
}

transformed parameters {
  real<lower=0, upper=5> con;
  real<lower=0, upper=2> gain;
  real<lower=0, upper=10> loss;
  real<lower=0, upper=1> update;
  
  con    = inv_logit(con_pr) * 5;
  gain   = inv_logit(gain_pr) * 2;
  loss   = inv_logit(loss_pr) * 10;
  update = inv_logit(update_pr);
}

model {
  con_pr    ~ normal(pr_mu[1], pr_sigma[1]);
  gain_pr   ~ normal(pr_mu[2], pr_sigma[2]);
  loss_pr   ~ normal(pr_mu[3], pr_sigma[3]);
  update_pr ~ normal(pr_mu[4], pr_sigma[4]);

  vector[4] ev = rep_vector(0., 4);
  real sensitivity = expm1(log(3) * con);

  ev = igt_model_lp(choice, shown, outcome, ev, T, sensitivity,
                    gain, loss, update);
}
