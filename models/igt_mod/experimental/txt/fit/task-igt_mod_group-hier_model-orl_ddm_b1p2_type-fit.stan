// Optimized Hierarchical ORL-DDM b1p2 Model for IGT_MOD
functions {
  real igt_pp_orl_ddm_model(
      array[] int choice, array[] int shown,
      array[] real outcome, array[] real RT,
      vector ev, vector ef, int Tsub,
      vector boundaries, vector taus,
      real beta, real sensitivity,
      real Arew, real Apun, real betaF
      ) {
    int curDeck;
    real log_lik = 0.0;

    vector[Tsub] drift_rates;
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
      curDeck = shown[t];
      drift_rates[t] = (local_ev[curDeck] + local_ef[curDeck] * betaF) * sensitivity;

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

    // Vectorized likelihood computation
    if (play_count > 0) {
      log_lik += wiener_lpdf(RT[play_indices[:play_count]] | 
                             boundaries[play_indices[:play_count]], 
                             taus[play_indices[:play_count]], 
                             beta, 
                             drift_rates[play_indices[:play_count]]);
    }

    if (pass_count > 0) {
      log_lik += wiener_lpdf(RT[pass_indices[:pass_count]] | 
                             boundaries[pass_indices[:pass_count]], 
                             taus[pass_indices[:pass_count]], 
                             1-beta, 
                             -drift_rates[pass_indices[: pass_count]]);
    }

    return log_lik;
  }

  // Main parallelization function
  real partial_sum(array[] int slice_n, int start, int end,
                   array[] int Tsubj, 
                   array[,] int choice, array[,] int shown, 
                   array[,] real outcome, array[,] real RT,
                   array[] vector boundary_subj,
                   array[] vector tau_subj,
                   array[] real beta, array[] real drift_con,
                   array[] real Arew, array[] real Apun, array[] real betaF) {
    real log_lik = 0.0;
    
    for (n in start:end) {
      vector[4] ev = rep_vector(0., 4);
      vector[4] ef = rep_vector(0., 4);
      real sensitivity = pow(3, drift_con[n]) - 1;
      
      log_lik += igt_pp_orl_ddm_model(
          choice[n, 1:Tsubj[n]], shown[n, 1:Tsubj[n]], 
          outcome[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]], 
          ev, ef, Tsubj[n],
          boundary_subj[n][1:Tsubj[n]],
	  tau_subj[n][1:Tsubj[n]], 
          beta[n], sensitivity,
	  Arew[n], Apun[n], betaF[n]
      );
    }
    return log_lik;
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
  array[N, T] int<lower=1, upper=4> shown;
  array[N, T] real outcome;
}

transformed data {
  int block = 20; 
  array[N] int subject_indices;
  for (i in 1:N) subject_indices[i] = i;
}

parameters {
  array[9] real mu_pr;
  array[9] real<lower=0> sigma;

  array[N] real boundary1_pr;
  array[N] real boundary_pr;
  array[N] real tau1_pr;
  array[N] real tau_pr;
  array[N] real beta_pr;
  array[N] real drift_con_pr;
  array[N] real Apun_pr;
  array[N] real Arew_pr;
  array[N] real betaF_pr;
}

transformed parameters {
  array[N] real<lower=0, upper=6> boundary1;
  array[N] real<lower=0, upper=6> boundary;
  array[N] real<lower=RTbound, upper=max(minRT)> tau1;
  array[N] real<lower=RTbound, upper=max(minRT)> tau;
  array[N] real<lower=0, upper=1> beta;
  array[N] real<lower=0, upper=3> drift_con;
  array[N] real<lower=0, upper=1> Apun;
  array[N] real<lower=0, upper=1> Arew;
  array[N] real betaF;

  boundary1 = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(boundary1_pr)) * 5 + 0.01);
  boundary  = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(boundary_pr)) * 5 + 0.01);
  tau1      = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(tau1_pr)) .* (to_vector(minRT) - RTbound - 1e-6) * 0.99 + RTbound);
  tau       = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(tau_pr)) .* (to_vector(minRT) - RTbound - 1e-6) * 0.99 + RTbound);
  beta      = to_array_1d(inv_logit(mu_pr[5] + sigma[5] .* to_vector(beta_pr)));
  drift_con = to_array_1d(inv_logit(mu_pr[6] + sigma[6] .* to_vector(drift_con_pr)) * 3);
  Apun      = to_array_1d(inv_logit(mu_pr[7] + sigma[7] .* to_vector(Apun_pr)));
  Arew      = to_array_1d(inv_logit(mu_pr[8] + sigma[8] .* to_vector(Arew_pr)));
  betaF     = to_array_1d(mu_pr[9] + sigma[9] .* to_vector(betaF_pr));
}

model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  boundary1_pr ~ normal(0, 1);
  boundary_pr  ~ normal(0, 1);
  tau1_pr      ~ normal(0, 1);
  tau_pr       ~ normal(0, 1);
  beta_pr      ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  Apun_pr      ~ normal(0, 1);
  Arew_pr      ~ normal(0, 1);
  betaF_pr     ~ normal(0, 1);

  // Build per-subject boundary/tau vectors
  array[N] vector[T] boundary_subj;
  array[N] vector[T] tau_subj;
  
  for (n in 1:N) {
    int Tsubj_n = Tsubj[n];
    
    // First block
    boundary_subj[n][1:block] = rep_vector(boundary1[n], block);
    tau_subj[n][1:block]      = rep_vector(tau1[n], block);
    
    // Rest of blocks
    int rest_len = Tsubj_n - block;
    boundary_subj[n][(block+1): Tsubj_n] = rep_vector(boundary[n], rest_len);
    tau_subj[n][(block+1): Tsubj_n]      = rep_vector(tau[n], rest_len);
  }

  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum,
                       subject_indices, grainsize,
                       Tsubj,
		       choice, shown,
		       outcome, RT,
                       boundary_subj, tau_subj,
                       beta, drift_con,
		       Arew, Apun, betaF);
}

generated quantities {
  real<lower=0, upper=6> mu_boundary1;
  real<lower=0, upper=6> mu_boundary;
  real<lower=0> mu_tau1;
  real<lower=0> mu_tau;
  real<lower=0, upper=1> mu_beta;
  real<lower=0, upper=5> mu_drift_con;
  real<lower=0, upper=1> mu_Apun;
  real<lower=0, upper=1> mu_Arew;
  real mu_betaF;
  
  mu_boundary1 = inv_logit(mu_pr[1]) * 5 + 0.01;
  mu_boundary  = inv_logit(mu_pr[2]) * 5 + 0.01;
  mu_tau1      = inv_logit(mu_pr[3]) * (mean(minRT) - RTbound) * 0.99 + RTbound;
  mu_tau       = inv_logit(mu_pr[4]) * (mean(minRT) - RTbound) * 0.99 + RTbound;
  mu_beta      = inv_logit(mu_pr[5]);
  mu_drift_con = inv_logit(mu_pr[6]) * 3;
  mu_Apun      = inv_logit(mu_pr[7]);
  mu_Arew      = inv_logit(mu_pr[8]);
  mu_betaF     = mu_pr[9];
}
