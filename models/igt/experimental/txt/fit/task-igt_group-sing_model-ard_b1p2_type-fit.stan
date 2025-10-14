// Hierarchical SSM-Only ARD Model for the Iowa Gambling Task
functions {
  // Defensive PDF for Racing Diffusion/Wald
  vector race_pdf_vec(vector t, vector boundary, vector drift) {
    int n = num_elements(t);
  
    // 1. Perform expensive math on whole vectors first
    vector[n] sqrt_t = sqrt(t);
    vector[n] denom = sqrt(2 * pi()) .* t .* sqrt_t;
    vector[n] boundary_over = boundary ./ denom;
    vector[n] drift_t_minus_boundary = drift .* t - boundary;
    vector[n] exponent = -0.5 * square(drift_t_minus_boundary) ./ t;
  
    vector[n] result;
    // 2. Now, loop ONLY for the conditional assignment
    for (i in 1:n) {
      if (exponent[i] < -30) {
        result[i] = 1e-10;
      } else {
        result[i] = boundary_over[i] * exp(exponent[i]);
      }
    }
    return result;
  }

  // Defensive CDF for Racing Diffusion/Wald
  vector race_cdf_vec(vector t, vector boundary, vector drift) {
    int n = num_elements(t);
  
    // 1. Perform expensive math on whole vectors first
    vector[n] sqrt_t = sqrt(t);
    vector[n] term1 = (drift .* t - boundary) ./ sqrt_t;
    vector[n] term2 = -(drift .* t + boundary) ./ sqrt_t;
    vector[n] expo_arg = 2.0 * drift .* boundary;
  
    vector[n] result;
    // 2. Loop for conditional logic and clamping
    for (i in 1:n) {
      if (expo_arg[i] > 30) {
        result[i] = Phi_approx(term1[i]);
      } else {
        result[i] = Phi_approx(term1[i]) + exp(expo_arg[i]) * Phi_approx(term2[i]);
      }
    
      // Clamp values
      if (result[i] <= 0) result[i] = 1e-10;
      if (result[i] >= 1) result[i] = 1 - 1e-10;
    }
    return result;
  }

  // Simplified ard_win_all using pre-calculated indices
  real ard_win_all(real RT, int choice, real tau, real boundary, vector drift_rates,
                 array[,] int win_indices_all, array[,] int lose_indices_all) {
    real t = RT - tau;
  
    array[3] int winning_indices = win_indices_all[choice];
    array[9] int losing_indices = lose_indices_all[choice];
  
    vector[3] pdf_winners = race_pdf_vec(
      rep_vector(t, 3), rep_vector(boundary, 3), drift_rates[winning_indices]
    );
    vector[9] cdf_losers = race_cdf_vec(
      rep_vector(t, 9), rep_vector(boundary, 9), drift_rates[losing_indices]
    );
  
    // Vectorized log-likelihood calculation is now safe and much faster
    return sum(log(pdf_winners)) + sum(log1m(cdf_losers));
  }

  // Optimized trial-level function
  real igt_ard_model(array[] int choice, array[] real RT, int T,
		     vector V_subj,
                     vector boundaries, vector taus, real urgency, real wd, real ws,
                     array[,] int win_indices_all, array[,] int lose_indices_all,
                     array[,] int other_indices) {
    real log_lik = 0.0;
    vector[12] drift_rates;
    int k = 1;
    
    // Calculate drift rates with minimum threshold
    for (i in 1:4) {
      drift_rates[k:k+2] = urgency + (ws + wd) * V_subj[i] + (ws - wd) * V_subj[other_indices[i]];
      k += 3;
    }

    for (t in 1:T) {
      if (RT[t] != 999) {
        real trial_lik = ard_win_all(RT[t], choice[t], taus[t], boundaries[t], drift_rates,
                                     win_indices_all, lose_indices_all);
        
        log_lik += trial_lik;
      }
    }
    return log_lik;
  }
}

//---

data {
  int<lower=1> sid;     		 // Subject ID
  int<lower=1> T;                         // Number of trials
  real<lower=0> minRT;                    // Minimum RT + small value to restrict tau
  real RTbound;                           // Lower bound of RT across all subjects
  array[T] int<lower=1, upper=4> choice;  // Choices made at each trial (1-4)
  array[T] real<lower=0> RT;              // Response times
}

//---

transformed data {
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

//---

parameters {
  // Subject-level raw parameters (z-scores for non-centered parameterization)
  real<lower=-5, upper=5> boundary1_pr;
  real<lower=-5, upper=5> boundary_pr;
  real tau1_pr;
  real tau_pr;
  real urgency_pr;
  real wd_pr;
  real ws_pr;
  // Static Value parameters
  real V1;
  real V2;
  real V3;
  real V4;
}

//---

transformed parameters {
  // Subject-level parameters
  real<lower=0, upper=5> boundary1;
  real<lower=0, upper=5> boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0> urgency;
  real<lower=0> wd;
  real<lower=0> ws;

  // ARD parameters
  boundary1 = inv_logit(boundary1_pr) * 4.99 + 0.001;
  boundary  = inv_logit(boundary_pr) * 4.99 + 0.001;
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound - 1e-6) * 0.95 + RTbound;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 1e-6) * 0.95 + RTbound;
  urgency   = log1p_exp(urgency_pr);
  wd        = log1p_exp(wd_pr);
  ws        = log1p_exp(ws_pr);
}

//---

model {
  // Priors
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  wd_pr ~ normal(0, 1);
  ws_pr ~ normal(0, 1);
  V1 ~ normal(0, 1);
  V2 ~ normal(0, 1);
  V3 ~ normal(0, 1);
  V4 ~ normal(0, 1);

  // Initialize values for the current subject
  vector[4] V_subj = [V1, V2, V3, V4]';

  // Create trial-varying boundary and tau vectors for this subject
  vector[T] boundaries = append_row(rep_vector(boundary1, 20), rep_vector(boundary, T - 20));
  vector[T] taus = append_row(rep_vector(tau1, 20), rep_vector(tau, T - 20));
    
  // Compute log-likelihood for the subject and add to the total
  target += igt_ard_model_lp(
        choice, RT, T,
        V_subj,
        boundaries, taus, urgency, wd, ws,
	win_indices_all, lose_indices_all, other_indices);

}
