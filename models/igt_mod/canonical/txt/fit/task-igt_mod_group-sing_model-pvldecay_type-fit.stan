functions {
  vector igt_model_lp(
      array[] int choice, array[] int shown, array[] real outcome,
      vector ev, int Tsub, real sensitivity,
      real gain, real loss, real decay
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

      local_ev *= (1 - decay);
      local_ev[curDeck] += curUtil * choice[t];
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
}

parameters {
  real con_pr;
  real gain_pr;
  real loss_pr;
  real decay_pr;
}

transformed parameters {
  real<lower=0, upper=5> con;
  real<lower=0, upper=2> gain;
  real<lower=0, upper=10> loss;
  real<lower=0, upper=1> decay;
  
  con    = inv_logit(con_pr) * 5;
  gain   = inv_logit(gain_pr) * 2;
  loss   = inv_logit(loss_pr) * 10;
  decay  = inv_logit(decay_pr);
}

model {
  con_pr    ~ normal(0, 1);
  gain_pr   ~ normal(0, 1);
  loss_pr   ~ normal(0, 1);
  decay_pr  ~ normal(0, 1);

  vector[4] ev = rep_vector(0., 4);
  real sensitivity = expm1(log(3) * con);

  ev = igt_model_lp(choice, shown, outcome, ev, T, sensitivity,
                    gain, loss, decay);
}
