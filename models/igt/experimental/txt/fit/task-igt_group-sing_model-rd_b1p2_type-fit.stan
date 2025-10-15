// Single-Subject RL-RD Model for the Iowa Gambling Task
functions {
  // Defensive PDF for the Wald (Inverse Gaussian) distribution
  vector race_pdf_vec(vector t, vector boundary, vector drift) {
    int n = num_elements(t);
    // 1. Perform expensive math on whole vectors first
    vector[n] sqrt_t = sqrt(t);
    vector[n] denom = sqrt(2 * pi()) .* t .* sqrt_t;
    vector[n] boundary_over = boundary ./ denom;
    vector[n] drift_t_minus_boundary = drift .* t - boundary;
    vector[n] exponent = -0.5 * square(drift_t_minus_boundary) ./ t;
  
    vector[n] result;
    // 2. Now, loop ONLY for the conditional assignment to prevent underflow
    for (i in 1:n) {
      if (exponent[i] < -30) {
        result[i] = 1e-10;
      } else {
        result[i] = boundary_over[i] * exp(exponent[i]);
      }
    }
    return result;
  }

  // Defensive CDF for the Wald (Inverse Gaussian) distribution
  vector race_cdf_vec(vector t, vector boundary, vector drift) {
    int n = num_elements(t);
    // 1. Perform expensive math on whole vectors first
    vector[n] sqrt_t = sqrt(t);
    vector[n] term1 = (drift .* t - boundary) ./ sqrt_t;
    vector[n] term2 = -(drift .* t + boundary) ./ sqrt_t;
    vector[n] expo_arg = 2.0 * drift .* boundary;
  
    vector[n] result;
    // 2. Loop for conditional logic and clamping to prevent overflow/underflow
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

  [cite_start]// Likelihood for a 4-way "win-first" race on a single trial [cite: 14]
  real win_first_lpdf(real RT, int choice, real tau, real boundary, vector drift_rates) {
    real t = RT - tau;
    if (t <= 1e-5) return negative_infinity(); [cite_start]// RT must be greater than non-decision time [cite: 15]

    array[3] int loser_indices;
    int k = 1;
    for (j in 1:4) {
      if (j != choice) {
        loser_indices[k] = j;
        k += 1;
      }
    }

    // Calculate log prob of winner and sum of log prob of losers not finishing
    [cite_start]real log_pdf_winner = log(race_pdf_vec(rep_vector(t, 1), rep_vector(boundary, 1), drift_rates[choice:choice])[1]); [cite: 17]
    [cite_start]vector[3] cdf_losers = race_cdf_vec(rep_vector(t, 3), rep_vector(boundary, 3), drift_rates[loser_indices]); [cite: 18]
    
    [cite_start]return log_pdf_winner + sum(log1m(cdf_losers)); [cite: 19]
  }

  // Trial-level function for the RD model
  real igt_rd_model(array[] int choice, array[] real RT, int T, vector V_subj,
                    vector boundaries, vector taus, real urgency, real drift_con) {
    real log_lik = 0.0;
    [cite_start]// Drift rates are constant across trials for this model [cite: 20]
    [cite_start]vector[4] drift_rates = urgency + drift_con * V_subj; [cite: 20]

    for (t in 1:T) {
      if (RT[t] != 999) { // Skip missed trials
        [cite_start]log_lik += win_first_lpdf(RT[t], choice[t], taus[t], boundaries[t], drift_rates); [cite: 21]
      }
    }
    return log_lik;
  }
}

//---

data {
  int<lower=1> sid;                         // Subject ID
  int<lower=1> T;                           // Number of trials
  real<lower=0> minRT;                      [cite_start]// Minimum RT for this subject [cite: 29]
  real<lower=0> RTbound;                    [cite_start]// Lower bound of RT across all subjects [cite: 29]
  array[T] int<lower=1, upper=4> choice;    [cite_start]// Choices made at each trial (1-4) [cite: 29]
  array[T] real<lower=0> RT;                [cite_start]// Response times [cite: 29]
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
