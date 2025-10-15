// Single-Subject RL-ALBA Model for the Iowa Gambling Task
functions {
  // Optimized PDF for the Linear Ballistic Accumulator
  vector lba_pdf_vec(vector t, real boundary, vector drift_mean, real sv, real A) {
    int n = num_elements(t);
    real local_A = fmax(A, 1e-10);

    vector[n] b_minus_A = rep_vector(boundary - local_A, n);
    vector[n] b = rep_vector(boundary, n);
    vector[n] z1 = (b - t .* drift_mean) / sv;
    vector[n] z2 = (b_minus_A - t .* drift_mean) / sv;
    
    // Vectorized standard normal PDF is more efficient
    vector[n] pdf_z1 = std_normal_pdf(z1);
    vector[n] pdf_z2 = std_normal_pdf(z2);

    vector[n] term1 = drift_mean .* (std_normal_cdf(z1) - std_normal_cdf(z2));
    vector[n] term2 = sv * (pdf_z2 - pdf_z1);
    
    vector[n] pdf = (term1 + term2) / local_A;
    // Defensive clamping for stability
    for(i in 1:n) if(pdf[i] < 1e-10) pdf[i] = 1e-10;
    return pdf;
  }

  // Optimized CDF for the Linear Ballistic Accumulator
  vector lba_cdf_vec(vector t, real boundary, vector drift_mean, real sv, real A) {
    int n = num_elements(t);
    real local_A = fmax(A, 1e-10);
    
    vector[n] b_minus_A = rep_vector(boundary - local_A, n);
    vector[n] b = rep_vector(boundary, n);
    vector[n] z1 = (b - t .* drift_mean) / sv;
    vector[n] z2 = (b_minus_A - t .* drift_mean) / sv;
    
    // Vectorized standard normal CDF is more efficient
    vector[n] cdf_z1 = std_normal_cdf(z1);
    vector[n] cdf_z2 = std_normal_cdf(z2);
    
    vector[n] pdf_z1 = std_normal_pdf(z1);
    vector[n] pdf_z2 = std_normal_pdf(z2);

    vector[n] cdf = 1 + ( (b - t .* drift_mean) .* cdf_z1 - 
                          (b_minus_A - t .* drift_mean) .* cdf_z2 - 
                          sv * (pdf_z1 - pdf_z2) ) / local_A;
    // Defensive clamping for stability
    for (i in 1:n) {
      if (cdf[i] <= 0) cdf[i] = 1e-10;
      if (cdf[i] >= 1) cdf[i] = 1 - 1e-10;
    }
    return cdf;
  }
  
  // Likelihood function calling the LBA pdf/cdf for a single trial
  real alba_win_all_lpdf(real RT, int choice, real tau, real boundary, vector drift_rates,
                         real A, real sv,
                         array[,] int win_indices_all, array[,] int lose_indices_all) {
    real t = RT - tau;
    if (t <= 1e-5) return negative_infinity(); // Return -inf if RT is less than non-decision time

    array[3] int winning_indices = win_indices_all[choice];
    array[9] int losing_indices = lose_indices_all[choice];
    
    vector[3] pdf_winners = lba_pdf_vec(rep_vector(t, 3), boundary, drift_rates[winning_indices], sv, A);
    vector[9] cdf_losers = lba_cdf_vec(rep_vector(t, 9), boundary, drift_rates[losing_indices], sv, A);
    
    return sum(log(pdf_winners)) + sum(log1m(cdf_losers));
  }
  
  // Trial-level function for the ALBA model
  real igt_alba_model(array[] int choice, array[] real RT, int T, vector V_subj,
                     vector boundaries, vector taus, real urgency, real wd, real ws,
                     real A, real sv,
                     array[,] int win_indices_all, array[,] int lose_indices_all,
                     array[,] int other_indices) {
    real log_lik = 0.0;
    vector[12] drift_rates;
    int k = 1;
    // Calculate drift rates once for all trials
    for (i in 1:4) {
      drift_rates[k:k+2] = urgency + (ws + wd) * V_subj[i] + (ws - wd) * V_subj[other_indices[i]];
      k += 3;
    }

    // Sum log-likelihood over trials
    for (t in 1:T) {
      if (RT[t] != 999) { // Skip missed trials
        log_lik += alba_win_all_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates,
                                     A, sv, win_indices_all, lose_indices_all);
      }
    }
    return log_lik;
  }
}

//---

data {
  int<lower=1> sid;                         // Subject ID
  int<lower=1> T;                           // Number of trials
  real<lower=0> minRT;                      // Minimum RT for this subject
  real<lower=0> RTbound;                    // Lower bound of RT across all subjects
  array[T] int<lower=1, upper=4> choice;    // Choices made at each trial (1-4)
  array[T] real<lower=0> RT;                // Response times
}

//---

transformed data {
  // Pre-calculate all possible indices once for efficiency
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

//---

parameters {
  // Subject-level raw parameters (z-scores for non-centered parameterization)
  real<lower=-5, upper=5> boundary1_pr;
  real<lower=-5, upper=5> boundary_pr;
  real<lower=-5, upper=5> tau1_pr;
  real<lower=-5, upper=5> tau_pr;
  real<lower=-5, upper=5> urgency_pr;
  real<lower=-5, upper=5> wd_pr;
  real<lower=-5, upper=5> ws_pr;
  real<lower=-5, upper=5> A_pr;
  real<lower=-5, upper=5> sv_pr;
  
  // Static Value parameters (sampled directly)
  real V1;
  real V2;
  real V3;
  real V4;
}

//---

transformed parameters {
  // Subject-level parameters transformed to their natural scale
  real<lower=0.001, upper=5> boundary1;
  real<lower=0.001, upper=5> boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0.001> urgency;
  real<lower=0.001> wd;
  real<lower=0.001> ws;
  real<lower=0.001> A;        // Starting point of evidence accumulation
  real<lower=0.001> sv;       // Drift rate variability

  // Transformations from _pr z-scores to model parameters
  boundary1 = inv_logit(boundary1_pr) * 4.99 + 0.001;
  boundary  = inv_logit(boundary_pr) * 4.99 + 0.001;
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  urgency   = log1p_exp(urgency_pr) + 0.01;
  wd        = log1p_exp(wd_pr) + 0.01;
  ws        = log1p_exp(ws_pr) + 0.01;
  A         = log1p_exp(A_pr) + 0.01;
  sv        = log1p_exp(sv_pr) + 0.01;
}

//---

model {
  // Priors on raw z-score parameters
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  wd_pr ~ normal(0, 1);
  ws_pr ~ normal(0, 1);
  A_pr ~ normal(0, 1);
  sv_pr ~ normal(0, 1);
  
  // Priors on deck value parameters
  V1 ~ normal(0, 1);
  V2 ~ normal(0, 1);
  V3 ~ normal(0, 1);
  V4 ~ normal(0, 1);

  // --- Likelihood Calculation ---
  // 1. Initialize values for the current subject
  vector[4] V_subj = [V1, V2, V3, V4]';
  
  // 2. Create trial-varying boundary and tau vectors
  vector[T] boundaries = append_row(rep_vector(boundary1, 20), rep_vector(boundary, T - 20));
  vector[T] taus = append_row(rep_vector(tau1, 20), rep_vector(tau, T - 20));
    
  // 3. Compute log-likelihood and add to the total
  target += igt_alba_model(
        choice, RT, T,
        V_subj,
        boundaries, taus, urgency, wd, ws, A, sv,
        win_indices_all, lose_indices_all, other_indices
  );
}
