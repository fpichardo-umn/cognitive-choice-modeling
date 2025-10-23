// Individual ORL-ARD Model for the Iowa Gambling Task
// Updated with improved numerical stability
functions {
  // Log PDF for Racing Diffusion/Wald (numerically stable)
  vector race_log_pdf_vec(vector t, vector boundary, vector drift) {
    int n = num_elements(t);
    vector[n] log_t = log(t);
    vector[n] sqrt_t = sqrt(t);
    
    // Compute log(boundary / (sqrt(2*pi) * t^(3/2)))
    vector[n] log_coef = log(boundary) - (0.5*log(2*pi()) + 1.5*log_t);
    
    // Exponent: -0.5 * (drift*t - boundary)^2 / t
    vector[n] exponent = -0.5 * square(drift .* t - boundary) ./ t;
    
    vector[n] log_pdf = log_coef + exponent;
    
    // Soft floor instead of hard threshold
    return fmax(log_pdf, -50);  // fmax is differentiable
  }

  // CDF for Racing Diffusion/Wald (numerically stable)
  vector race_cdf_vec(vector t, vector boundary, vector drift) {
    int n = num_elements(t);
    vector[n] sqrt_t = sqrt(t);
    vector[n] term1 = (drift .* t - boundary) ./ sqrt_t;
    vector[n] term2 = -(drift .* t + boundary) ./ sqrt_t;
    vector[n] expo_arg = 2.0 * drift .* boundary;
    vector[n] result;
    
    for (i in 1:n) {
      // For numerical safety when expo_arg is huge
      if (expo_arg[i] > 30) {
        result[i] = Phi_approx(term1[i]);
      } else {
        result[i] = Phi_approx(term1[i]) + exp(expo_arg[i]) * Phi_approx(term2[i]);
      }
    }
    
    // Clamp using differentiable functions
    return fmin(fmax(result, 1e-10), 1 - 1e-10);
  }

  // ARD likelihood: all 3 accumulators for chosen deck must win
  // Vectorized for better numerical stability
  real ard_win_all(real RT, int choice, real tau, real boundary, vector drift_rates,
                   array[,] int win_indices_all, array[,] int lose_indices_all) {
    real t = fmax(RT - tau, 1e-3);
    if (t <= 1e-5) return negative_infinity();
  
    array[3] int winning_indices = win_indices_all[choice];
    array[9] int losing_indices = lose_indices_all[choice];
  
    // Get log-PDFs for all 3 winners using vectorization
    vector[3] log_pdf_winners = race_log_pdf_vec(
      rep_vector(t, 3), 
      rep_vector(boundary, 3), 
      drift_rates[winning_indices]
    );
    
    // Get CDFs for all 9 losers using vectorization
    vector[9] cdf_losers = race_cdf_vec(
      rep_vector(t, 9), 
      rep_vector(boundary, 9), 
      drift_rates[losing_indices]
    );
  
    // Combine: sum of log-PDFs for winners + sum of log(1-CDF) for losers
    return sum(log_pdf_winners) + sum(log1m(cdf_losers));
  }

  // ORL-ARD trial-level function
  real igt_ard_model(
      array[] int choice, array[] real wins, array[] real losses, array[] real RT,
      vector ev, vector ef, int T,
      real sensitivity, real Arew, real Apun, real K,
      real betaF, real betaP,
      vector boundaries, vector taus, real urgency, real wd, real ws,
      array[,] int win_indices_all,
      array[,] int lose_indices_all,
      array[,] int other_indices) {

    vector[4] local_ev = ev;
    vector[4] local_ef = ef;
    vector[4] pers = rep_vector(0.0, 4);
    real log_lik = 0.0;
    real PEval;
    real PEfreq;
    vector[4] PEfreq_fic;
    array[T] real sign_outcome;
    real K_tr = pow(3, K) - 1;

    real scaled_urgency = urgency * sensitivity;
    real scaled_wswd_plus = (ws + wd) * sensitivity;
    real scaled_wswd_minus = (ws - wd) * sensitivity;
    real scaled_betaF = betaF * sensitivity;
    real scaled_betaP = betaP * sensitivity;

    for (t in 1:T) {
      sign_outcome[t] = wins[t] >= losses[t] ? 1.0 : -1.0;
    }

    for (t in 1:T) {
      vector[12] drift_rates;
      int k = 1;
      for (i in 1:4) {
        real combined_value = local_ev[i] + local_ef[i] * scaled_betaF + pers[i] * scaled_betaP;
        drift_rates[k:k+2] = scaled_urgency +
                             scaled_wswd_plus * combined_value +
                             scaled_wswd_minus * (local_ev[other_indices[i]] + 
                                                   local_ef[other_indices[i]] * scaled_betaF + 
                                                   pers[other_indices[i]] * scaled_betaP);
        k += 3;
      }

      if (RT[t] != 999) {
        log_lik += ard_win_all(RT[t], choice[t], taus[t], boundaries[t], drift_rates,
                               win_indices_all, lose_indices_all);
      }

      PEval = wins[t] - losses[t] - local_ev[choice[t]];
      PEfreq = sign_outcome[t] - local_ef[choice[t]];
      
      for (d in 1:4) {
        PEfreq_fic[d] = -sign_outcome[t]/3.0 - local_ef[d];
      }
      
      if (wins[t] >= losses[t]) {
        local_ef += Apun * PEfreq_fic;
        local_ef[choice[t]] = local_ef[choice[t]] + Arew * PEfreq;
        local_ev[choice[t]] = local_ev[choice[t]] + Arew * PEval;
      } else {
        local_ef += Arew * PEfreq_fic;
        local_ef[choice[t]] = local_ef[choice[t]] + Apun * PEfreq;
        local_ev[choice[t]] = local_ev[choice[t]] + Apun * PEval;
      }
      
      pers[choice[t]] = 1;
      pers = pers / (1 + K_tr);
    }
    return log_lik;
  }

  // Main parallelization function
  real partial_sum(array[] int slice_n, int start, int end,
                   array[] int Tsubj, array[,] int choice,
                   array[,] real wins, array[,] real losses, array[,] real RT,
                   array[,] int win_indices_all, array[,] int lose_indices_all,
                   array[,] int other_indices,
                   array[] vector boundary_subj,
                   array[] vector tau_subj,
                   array[] real drift_con, array[] real Arew, array[] real Apun, array[] real K,
                   array[] real betaF, array[] real betaP,
                   array[] real urgency, array[] real wd, array[] real ws) {
    real log_lik = 0.0;
    vector[4] ev = rep_vector(0.0, 4);
    vector[4] ef = rep_vector(0.0, 4);

    for (n in start:end) {
      real sensitivity = pow(3, drift_con[n]) - 1;

      log_lik += igt_ard_model(
          choice[n, 1:Tsubj[n]], wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]], RT[n, 1:Tsubj[n]],
          ev, ef, Tsubj[n],
          sensitivity, Arew[n], Apun[n], K[n],
          betaF[n], betaP[n],
          boundary_subj[n][1:Tsubj[n]], tau_subj[n][1:Tsubj[n]], urgency[n], wd[n], ws[n],
          win_indices_all, lose_indices_all, other_indices
      );
    }
    return log_lik;
  }
}

data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1> sid;
  array[N] int<lower=1> Tsubj;
  array[N] real<lower=0> minRT;
  real<lower=0> RTbound;
  array[N, T] int<lower=1, upper=4> choice;
  array[N, T] real RT;
  array[N, T] real wins;
  array[N, T] real losses;
}

transformed data {
  array[N] int subject_indices;
  for (i in 1:N) {
    subject_indices[i] = i;
  }
  int block = 20;
  
  // Pre-calculate all possible indices once
  array[4, 3] int win_indices_all;
  array[4, 9] int lose_indices_all;
  array[4, 3] int other_indices = { {2, 3, 4}, {1, 3, 4}, {1, 2, 4}, {1, 2, 3} };

  for (c in 1:4) {
    int losing_idx = 1;
    for (i in 1:3) {
      win_indices_all[c, i] = (c - 1) * 3 + i;
    }
    for (j in 1:4) {
      if (j != c) {
        for (i in 1:3) {
          lose_indices_all[c, losing_idx] = (j - 1) * 3 + i;
          losing_idx += 1;
        }
      }
    }
  }
}

parameters {
  // Individual-level parameters (non-centered, N subjects)
  vector[N] boundary1_pr;
  vector[N] boundary_pr;
  vector[N] tau1_pr;
  vector[N] tau_pr;
  vector[N] urgency_pr;
  vector[N] wd_pr;
  vector[N] ws_pr;
  vector[N] drift_con_pr;
  vector[N] Arew_pr;
  vector[N] Apun_pr;
  vector[N] K_pr;
  vector[N] betaF_pr;
  vector[N] betaP_pr;
}

transformed parameters {
  // Transformed individual-level parameters (N subjects)
  array[N] real<lower=0.001, upper=5> boundary1;
  array[N] real<lower=0.001, upper=5> boundary;
  array[N] real<lower=0> tau1;
  array[N] real<lower=0> tau;
  array[N] real<lower=0.001, upper=20> urgency;
  array[N] real<lower=0.001, upper=10> wd;
  array[N] real<lower=0.001, upper=10> ws;
  array[N] real<lower=0, upper=5> drift_con;
  array[N] real<lower=0, upper=1> Arew;
  array[N] real<lower=0, upper=1> Apun;
  array[N] real<lower=0, upper=5> K;
  array[N] real betaF;
  array[N] real betaP;
  
  // Vectorized transforms
  vector[N] minRT_vec = to_vector(minRT);
  boundary1 = to_array_1d(inv_logit(boundary1_pr) * 4.99 + 0.001);
  boundary  = to_array_1d(inv_logit(boundary_pr) * 4.99 + 0.001);
  tau1 = to_array_1d(inv_logit(tau1_pr) .* (minRT_vec - RTbound - 0.02) * 0.95 + RTbound);
  tau  = to_array_1d(inv_logit(tau_pr) .* (minRT_vec - RTbound - 0.02) * 0.95 + RTbound);

  urgency = to_array_1d(inv_logit(urgency_pr) * 19.999 + 0.001);
  wd = to_array_1d(inv_logit(wd_pr) * 9.999 + 0.001);
  ws = to_array_1d(inv_logit(ws_pr) * 9.999 + 0.001);

  drift_con = to_array_1d(inv_logit(drift_con_pr) * 5);
  Arew      = to_array_1d(inv_logit(Arew_pr));
  Apun      = to_array_1d(inv_logit(Apun_pr));
  K         = to_array_1d(inv_logit(K_pr) * 5);
  betaF     = to_array_1d(betaF_pr);
  betaP     = to_array_1d(betaP_pr);
}

model {
  // Priors on non-centered parameters (for all N subjects)
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  wd_pr ~ normal(0, 1);
  ws_pr ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  Arew_pr ~ normal(0, 1);
  Apun_pr ~ normal(0, 1);
  K_pr ~ normal(0, 1);
  betaF_pr ~ normal(0, 1);
  betaP_pr ~ normal(0, 1);

  // Subject- and trial-specific vectors for the likelihood
  // This is an array of N vectors, each of max length T
  array[N] vector[T] boundary_subj;
  array[N] vector[T] tau_subj;

  // Build boundary/tau vectors for all subjects
  for (n in 1:N) {
    // First block
    boundary_subj[n][1:block] = rep_vector(boundary1[n], block);
    tau_subj[n][1:block]      = rep_vector(tau1[n], block);
    
    // Rest of trials (if any)
    if (Tsubj[n] > block) {
      int rest_len = Tsubj[n] - block;
      boundary_subj[n][(block+1):Tsubj[n]] = rep_vector(boundary[n], rest_len);
      tau_subj[n][(block+1):Tsubj[n]]      = rep_vector(tau[n], rest_len);
    }
  }
  
  // Likelihood using parallelization
  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum, subject_indices, grainsize,
                       Tsubj, choice, wins, losses, RT,
                       win_indices_all, lose_indices_all, other_indices,
                       boundary_subj, tau_subj,
                       drift_con, Arew, Apun, K,
                       betaF, betaP,
                       urgency, wd, ws);
}
