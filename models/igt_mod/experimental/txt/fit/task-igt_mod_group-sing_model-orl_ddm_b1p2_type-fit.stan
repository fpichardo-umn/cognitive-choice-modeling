functions {
  vector igt_pp_orl_ddm_model_lp(
      array[] int choice, array[] int shown, array[] real outcome,
      array[] real RT, vector ev, vector ef, int Tsub,
      real Arew, real Apun, real betaF, real sensitivity,
      real boundary1, real boundary, real tau1, real tau, real beta
      ) {
    int curDeck;
    vector[Tsub] drift_rates;
    vector[Tsub] boundaries;
    vector[Tsub] nondt;
    vector[4] local_ev = ev;
    vector[4] local_ef = ef;
    real PEval;
    real PEfreq;
    real efChosen; 
    vector[4] PEfreq_fic;
    real sign_outcome;

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
      drift_rates[t] = (local_ev[curDeck] + local_ef[curDeck] * betaF) * sensitivity;

      if (choice[t] == 1) {
        sign_outcome = (outcome[t] >= 0) ? 1.0 : -1.0;
        PEval = outcome[t] - local_ev[curDeck];
        PEfreq = sign_outcome - local_ef[curDeck];
        efChosen = local_ef[curDeck];
        
        for (d in 1:4) {
          PEfreq_fic[d] = -sign_outcome / 3.0 - local_ef[d];
        }
        
        if (outcome[t] >= 0) {
          local_ef += Apun * PEfreq_fic;
          local_ef[curDeck] = efChosen + Arew * PEfreq; // Correct the chosen deck using the stored value
          local_ev[curDeck] += Arew * PEval;
        } else {
          local_ef += Arew * PEfreq_fic;
          local_ef[curDeck] = efChosen + Apun * PEfreq; // Correct the chosen deck using the stored value
          local_ev[curDeck] += Apun * PEval;
        }
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
  int<lower=1> sid;
  int<lower=1> T;
  real<lower=0> minRT;
  real RTbound;
  array[T] real<lower=0> RT;
  array[T] int<lower=0, upper=1> choice;
  array[T] int<lower=1, upper=4> shown;
  array[T] real outcome;
}

parameters {
  real<lower=-5, upper=5> boundary1_pr;
  real<lower=-5, upper=5> boundary_pr;
  real<lower=-3, upper=3> tau1_pr;
  real<lower=-3, upper=3> tau_pr;
  real<lower=-3, upper=3> beta_pr;
  real<lower=-3, upper=3> Apun_pr;
  real<lower=-3, upper=3> Arew_pr;
  real<lower=-3, upper=3> betaF_pr;
}

transformed parameters {
  real<lower=0, upper=6> boundary1;
  real<lower=0, upper=6> boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0, upper=1> beta;
  real<lower=0, upper=1> Apun;
  real<lower=0, upper=1> Arew;
  real betaF;

  boundary1 = inv_logit(boundary1_pr) * 5 + 0.01;
  boundary  = inv_logit(boundary_pr) * 5 + 0.01;
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound;
  beta      = inv_logit(beta_pr);
  Apun      = inv_logit(Apun_pr);
  Arew      = inv_logit(Arew_pr);
  betaF     = betaF_pr;
}

model {
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  beta_pr ~ normal(0, 1);
  Apun_pr ~ normal(0, 1);
  Arew_pr ~ normal(0, 1);
  betaF_pr ~ normal(0, 1);

  vector[4] ev = rep_vector(0., 4);
  vector[4] ef = rep_vector(0., 4);
  real sensitivity = 0.15;
  
  ev = igt_pp_orl_ddm_model_lp(choice, shown, outcome, RT, ev, ef, T,
                                Arew, Apun, betaF, sensitivity,
                                boundary1, boundary, tau1, tau, beta);
}
