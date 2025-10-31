functions {
  vector igt_model_lp(
      array[] int choice, array[] int shown, array[] real outcome,
      vector ev, int Tsub, real sensitivity,
      real update, real wgt_pun, real wgt_rew
      ) {
    real curUtil;
    int curDeck;
    vector[Tsub] Info;
    vector[4] local_ev = ev;

    for (t in 1:Tsub) {
      curDeck = shown[t];
      Info[t] = sensitivity * local_ev[curDeck];
      curUtil = ((outcome[t] > 0 ? wgt_rew : wgt_pun)) * outcome[t];
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
}

parameters {
  real con_pr;
  real wgt_pun_pr;
  real wgt_rew_pr;
  real update_pr;
}

transformed parameters {
  real<lower=0, upper=5> con;
  real<lower=0, upper=1> wgt_pun;
  real<lower=0, upper=1> wgt_rew;
  real<lower=0, upper=1> update;
  
  con     = inv_logit(con_pr) * 5;
  wgt_pun = inv_logit(wgt_pun_pr);
  wgt_rew = inv_logit(wgt_rew_pr);
  update  = inv_logit(update_pr);
}

model {
  con_pr     ~ normal(0, 1);
  wgt_pun_pr ~ normal(0, 1);
  wgt_rew_pr ~ normal(0, 1);
  update_pr  ~ normal(0, 1);

  vector[4] ev = rep_vector(0., 4);
  real sensitivity = expm1(log(3) * con);

  ev = igt_model_lp(choice, shown, outcome, ev, T, sensitivity,
                    update, wgt_pun, wgt_rew);
}
