// Optimized Hybrid Hierarchical ORL Model for the IGT
functions {
  real igt_subject(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, vector ef, int Tsub,
        real Arew, real Apun, real K, real betaF, real betaP
        ) {
    real log_lik = 0.0;
    vector[4] local_ev = ev;
    vector[4] local_ef = ef;
    vector[4] pers = rep_vector(0.0, 4);
    vector[4] util;
    real PEval;
    real PEfreq;
    real efChosen;
    vector[4] PEfreq_fic;
    array[Tsub] real sign_outcome;
    real K_tr = pow(3, K) - 1;

    for (t in 1:Tsub) {
      sign_outcome[t] = wins[t] >= losses[t] ? 1.0 : -1.0;
    }

    for (t in 1:Tsub) {
      util = local_ev + local_ef * betaF + pers * betaP;
      log_lik += categorical_logit_lpmf(choice[t] | util);
      
      PEval = wins[t] - losses[t] - local_ev[choice[t]];
      PEfreq = sign_outcome[t] - local_ef[choice[t]];
      efChosen = local_ef[choice[t]];

      for (d in 1:4) {
        PEfreq_fic[d] = -sign_outcome[t]/3.0 - local_ef[d];
      }
      
      if (wins[t] >= losses[t]) {
        local_ef += Apun * PEfreq_fic;
        local_ef[choice[t]] = efChosen + Arew * PEfreq;
        local_ev[choice[t]] = local_ev[choice[t]] + Arew * PEval;
      } else {
        local_ef += Arew * PEfreq_fic;
        local_ef[choice[t]] = efChosen + Apun * PEfreq;
        local_ev[choice[t]] = local_ev[choice[t]] + Apun * PEval;
      }
      
      pers[choice[t]] = 1;
      pers = pers / (1 + K_tr);
    }
    return log_lik;
  }
  
  real partial_sum_func(array[] int slice_n, int start, int end,
                        array[,] int choice, array[,] real wins, array[,] real losses, 
                        array[] int Tsubj,
                        array[] real Arew, array[] real Apun, array[] real K,
                        array[] real betaF, array[] real betaP) {
    real log_lik = 0.0;
    for (n in start:end) {
      vector[4] ev = rep_vector(0., 4);
      vector[4] ef = rep_vector(0., 4);
      log_lik += igt_subject(choice[n, 1:Tsubj[n]],
                             wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]],
                             ev, ef, Tsubj[n], 
                             Arew[n], Apun[n], K[n], 
                             betaF[n], betaP[n]);
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
  // 1. NCP Hyperparameters (Arew & betaP)
  vector[2] mu_pr_rest;
  vector<lower=0>[2] sigma_rest;

  // 2. CP Hyperparameters (Apun, K, betaF)
  real<lower=0, upper=1> mu_Apun;
  real<lower=0> sigma_Apun;
  real<lower=0, upper=5> mu_K;
  real<lower=0> sigma_K;
  real mu_betaF;
  real<lower=0> sigma_betaF;

  // 3. Individual parameters
  vector[N] Arew_pr;
  vector[N] betaP_pr;
  
  array[N] real<lower=0, upper=1> Apun; // Centered
  array[N] real<lower=0, upper=5> K;    // Centered
  array[N] real betaF;                 // Centered
}

transformed parameters {
  array[N] real<lower=0, upper=1> Arew;
  array[N] real betaP;
  
  // Non-Centered reparameterization
  Arew  = to_array_1d(inv_logit(mu_pr_rest[1] + sigma_rest[1] * Arew_pr));
  betaP = to_array_1d(mu_pr_rest[2] + sigma_rest[2] * betaP_pr);
}

model {
  // Tighter Global Priors for stability
  mu_pr_rest ~ normal(0, 1);
  sigma_rest ~ normal(0, 0.5); 
  
  mu_Apun ~ beta(2, 2);
  sigma_Apun ~ normal(0, 0.5);
  mu_K ~ normal(2.5, 1);
  sigma_K ~ normal(0, 0.5);
  mu_betaF ~ normal(0, 2);
  sigma_betaF ~ normal(0, 1);

  Arew_pr  ~ normal(0, 1);
  betaP_pr ~ normal(0, 1);

  // Centered Likelihoods
  Apun  ~ normal(mu_Apun, sigma_Apun);
  K     ~ normal(mu_K, sigma_K);
  betaF ~ normal(mu_betaF, sigma_betaF);

  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       Arew, Apun, K, betaF, betaP);
}

generated quantities {
  real mu_Arew = inv_logit(mu_pr_rest[1]);
  real mu_betaP = mu_pr_rest[2];
  // Other means are already in the parameters block
}
