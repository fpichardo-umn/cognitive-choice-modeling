functions {
  vector igt_model_lp(
        array[] int choice,
        array[] int shown,
        array[] real outcome,
        vector ev,
        vector ef,
        int Tsub,
        real Arew,
        real Apun,
        real betaF,
        real betaB
        ) {
    vector[4] local_ev = ev;
    vector[4] local_ef = ef;
    vector[Tsub] Info;
    int curDeck;
    real PEval;
    real PEfreq;
    real efChosen; 
    vector[4] PEfreq_fic;
    real sign_outcome;

    for (t in 1:Tsub) {
      curDeck = shown[t];
      Info[t] = local_ev[curDeck] + local_ef[curDeck] * betaF + betaB;

      if (choice[t] == 1) {
        sign_outcome = (outcome[t] > 0) ? 1.0 : -1.0;
        PEval = outcome[t] - local_ev[curDeck];
        PEfreq = sign_outcome - local_ef[curDeck];
        efChosen = local_ef[choice[t]];
        
        for (d in 1:4) {
          PEfreq_fic[d] = -sign_outcome / 3.0 - local_ef[d];
        }
        
        if (outcome[t] > 0) {
          local_ef += Apun * PEfreq_fic;
        local_ef[choice[t]] = efChosen + Apun * PEfreq; // Correct the chosen deck using the stored value
          local_ev[curDeck] += Arew * PEval;
        } else {
          local_ef += Arew * PEfreq_fic;
        local_ef[choice[t]] = efChosen + Apun * PEfreq; // Correct the chosen deck using the stored value
          local_ev[curDeck] += Apun * PEval;
        }
      }
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
  real Arew_pr;
  real Apun_pr;
  real betaF;
  real betaB;
}

transformed parameters {
  real<lower=0, upper=1> Arew;
  real<lower=0, upper=1> Apun;

  Arew = inv_logit(Arew_pr);
  Apun = inv_logit(Apun_pr);
}

model {
  Arew_pr ~ normal(0, 1);
  Apun_pr ~ normal(0, 1);
  betaF   ~ normal(0, 1);
  betaB   ~ normal(0, 1);

  vector[4] ev = rep_vector(0., 4);
  vector[4] ef = rep_vector(0., 4);
    
  ev = igt_model_lp(choice, shown, outcome, ev, ef, T,
                    Arew, Apun, betaF, betaB);
}
