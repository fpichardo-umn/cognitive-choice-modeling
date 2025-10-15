// Single-Subject RL-RD Model for the Iowa Gambling Task
functions {
  // Defensive PDF for Racing Diffusion/Wald
  vector race_pdf_vec(vector t, vector boundary, vector drift) {
    int n = num_elements(t);
    vector[n] sqrt_t = sqrt(t);
    vector[n] denom = sqrt(2 * pi()) .* t .* sqrt_t;
    vector[n] boundary_over = boundary ./ denom;
    vector[n] drift_t_minus_boundary = drift .* t - boundary;
    vector[n] exponent = -0.5 * square(drift_t_minus_boundary) ./ t;
    vector[n] result;
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
    vector[n] sqrt_t = sqrt(t);
    vector[n] term1 = (drift .* t - boundary) ./ sqrt_t;
    vector[n] term2 = -(drift .* t + boundary) ./ sqrt_t;
    vector[n] expo_arg = 2.0 * drift .* boundary;
    vector[n] result;
    for (i in 1:n) {
      if (expo_arg[i] > 30) {
        result[i] = Phi_approx(term1[i]);
      } else {
        result[i] = Phi_approx(term1[i]) + exp(expo_arg[i]) * Phi_approx(term2[i]);
      }
      if (result[i] <= 0) result[i] = 1e-10;
      if (result[i] >= 1) result[i] = 1 - 1e-10;
    }
    return result;
  }

  // Simplified likelihood for a 4-way "win-first" race
  real win_first_lpdf(real RT, int choice, real tau, real boundary, vector drift_rates) {
    real t = RT - tau;
    if (t <= 1e-5) return negative_infinity(); // Use a small threshold for safety

    array[3] int loser_indices;
    int k = 1;
    for (j in 1:4) {
      if (j != choice) {
        loser_indices[k] = j;
        k += 1;
      }
    }

    real log_pdf_winner = log(race_pdf_vec(rep_vector(t, 1), rep_vector(boundary, 1), drift_rates[choice:choice])[1]);
    vector[3] cdf_losers = race_cdf_vec(rep_vector(t, 3), rep_vector(boundary, 3), drift_rates[loser_indices]);
    
    return log_pdf_winner + sum(log1m(cdf_losers));
  }

  // Trial-level function for the RD model
  real igt_rd_model(array[] int choice, array[] real RT, int T, vector V_subj,
                    vector boundaries, vector taus, real urgency, real drift_con) {
    real log_lik = 0.0;
    // Drift rates are constant across trials for this model
    vector[4] drift_rates = urgency + drift_con * V_subj;

    for (t in 1:T) {
      if (RT[t] != 999) { // Skip missed trials
        log_lik += win_first_lpdf(RT[t] | choice[t], taus[t], boundaries[t], drift_rates);
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

parameters {
  // Subject-level raw parameters (z-scores for non-centered parameterization)
  real<lower=-5, upper=5> boundary1_pr;
  real<lower=-5, upper=5> boundary_pr;
  real<lower=-5, upper=5> tau1_pr;
  real<lower=-5, upper=5> tau_pr;
  real<lower=-5, upper=5> urgency_pr;
  real<lower=-5, upper=5> drift_con_pr; // Drift rate consistency/scaling
  
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
  real<lower=0.001> drift_con;

  // Transformations from _pr z-scores to model parameters
  boundary1 = inv_logit(boundary1_pr) * 4.99 + 0.001;
  boundary  = inv_logit(boundary_pr) * 4.99 + 0.001;
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  urgency   = log1p_exp(urgency_pr) + 0.01;
  drift_con = log1p_exp(drift_con_pr) + 0.01;
}

//---

model {
  // Priors on raw z-score parameters
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  
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
  target += igt_rd_model(
        choice, RT, T,
        V_subj,
        boundaries, taus, urgency, drift_con
  );
}
