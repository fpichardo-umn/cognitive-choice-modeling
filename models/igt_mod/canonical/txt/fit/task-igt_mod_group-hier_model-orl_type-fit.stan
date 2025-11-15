// Optimized Hierarchical ORL Model for IGT_MOD
functions {
  // Returns log-likelihood instead of modifying target
  real igt_subject(
        array[] int choice, array[] int shown, array[] real outcome,
        vector ev, vector ef, int Tsub,
        real Arew, real Apun, real betaF, real betaB
        ) {
    real log_lik = 0.0;
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
        sign_outcome = (outcome[t] >= 0) ? 1.0 : -1.0;
        PEval = outcome[t] - local_ev[curDeck];
        PEfreq = sign_outcome - local_ef[curDeck];
        efChosen = local_ef[choice[t]];
        
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
    }
    
    log_lik += bernoulli_logit_lpmf(choice | Info);
    return log_lik;
  }
  
  // Parallelization wrapper
  real partial_sum_func(array[] int slice_n, int start, int end,
                        array[,] int choice, array[,] int shown, array[,] real outcome,
                        array[] int Tsubj,
                        array[] real Arew, array[] real Apun,
                        array[] real betaF, array[] real betaB) {
    real log_lik = 0.0;
    
    for (n in start:end) {
      vector[4] ev = rep_vector(0., 4);
      vector[4] ef = rep_vector(0., 4);
      
      log_lik += igt_subject(choice[n, 1:Tsubj[n]], shown[n, 1:Tsubj[n]],
                             outcome[n, 1:Tsubj[n]], ev, ef, Tsubj[n],
                             Arew[n], Apun[n], betaF[n], betaB[n]);
    }
    return log_lik;
  }
}

data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1> Tsubj;
  array[N, T] int<lower=0, upper=1> choice;
  array[N, T] int<lower=0, upper=4> shown;
  array[N, T] real outcome;
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

  array[N] real Arew_pr;
  array[N] real Apun_pr;
  array[N] real betaF_pr;
  array[N] real betaB_pr;
}

transformed parameters {
  array[N] real<lower=0, upper=1> Arew;
  array[N] real<lower=0, upper=1> Apun;
  array[N] real betaF;
  array[N] real betaB;

  Arew  = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(Arew_pr)));
  Apun  = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(Apun_pr)));
  betaF = to_array_1d(mu_pr[3] + sigma[3] .* to_vector(betaF_pr));
  betaB = to_array_1d(mu_pr[4] + sigma[4] .* to_vector(betaB_pr));
}

model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  Arew_pr  ~ normal(0, 1);
  Apun_pr  ~ normal(0, 1);
  betaF_pr ~ normal(0, 1);
  betaB_pr ~ normal(0, 1);

  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, shown, outcome, Tsubj,
                       Arew, Apun, betaF, betaB);
}

generated quantities {
  real<lower=0, upper=1> mu_Arew;
  real<lower=0, upper=1> mu_Apun;
  real mu_betaF;
  real mu_betaB;
  
  mu_Arew  = inv_logit(mu_pr[1]);
  mu_Apun  = inv_logit(mu_pr[2]);
  mu_betaF = mu_pr[3];
  mu_betaB = mu_pr[4];
}
