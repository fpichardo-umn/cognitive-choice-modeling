// Individual Static-ARD Model for the Iowa Gambling Task
// Pure SSM with static deck preferences (no learning)
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

  // Static-ARD trial-level function (no learning, static preferences)
  real igt_ard_model(
      array[] int choice, array[] real RT,
      int T,
      real sensitivity, vector deck_values,
      vector boundaries, vector taus, real urgency,
      array[,] int win_indices_all,
      array[,] int lose_indices_all) {

    real log_lik = 0.0;
    real scaled_urgency = urgency * sensitivity;

    for (t in 1:T) {
      vector[12] drift_rates;
      int k = 1;
      for (i in 1:4) {
        // All 3 accumulators for deck i have the same drift
        // based on static preference V[i]
        for (j in 1:3) {
          drift_rates[k] = scaled_urgency + deck_values[i] * sensitivity;
          k += 1;
        }
      }

      if (RT[t] != 999) {
        log_lik += ard_win_all(RT[t], choice[t], taus[t], boundaries[t], drift_rates,
                               win_indices_all, lose_indices_all);
      }
    }
    return log_lik;
  }
}

data {
  int<lower=1> T;                           // Total number of trials
  real<lower=0> minRT;                      // Minimum RT for this subject
  real<lower=0> RTbound;                    // RT bound
  array[T] int<lower=1, upper=4> choice;   // Choices
  array[T] real RT;                         // Response times
}

transformed data {
  int block = 20;
  
  // Pre-calculate all possible indices once
  array[4, 3] int win_indices_all;
  array[4, 9] int lose_indices_all;

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
  // Individual-level parameters (non-centered)
  real boundary1_pr;
  real boundary_pr;
  real tau1_pr;
  real tau_pr;
  real urgency_pr;
  real drift_con_pr;
  
  // Static deck preferences
  real V1_pr;
  real V2_pr;
  real V3_pr;
  real V4_pr;
}

transformed parameters {
  // Transformed individual-level parameters
  real<lower=0.001, upper=5> boundary1;
  real<lower=0.001, upper=5> boundary;
  real<lower=0> tau1;
  real<lower=0> tau;
  real<lower=0.001, upper=20> urgency;
  real<lower=0, upper=5> drift_con;
  
  // Static deck preferences (can be negative or positive)
  real<lower=-10, upper=10> V1;
  real<lower=-10, upper=10> V2;
  real<lower=-10, upper=10> V3;
  real<lower=-10, upper=10> V4;
  
  // Transform parameters
  boundary1 = inv_logit(boundary1_pr) * 4.99 + 0.001;
  boundary  = inv_logit(boundary_pr) * 4.99 + 0.001;
  tau1 = inv_logit(tau1_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  tau  = inv_logit(tau_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;

  urgency = inv_logit(urgency_pr) * 19.999 + 0.001;
  drift_con = inv_logit(drift_con_pr) * 5;
  
  // Transform deck preferences from (-10, 10) range
  V1 = inv_logit(V1_pr) * 20 - 10;
  V2 = inv_logit(V2_pr) * 20 - 10;
  V3 = inv_logit(V3_pr) * 20 - 10;
  V4 = inv_logit(V4_pr) * 20 - 10;
}

model {
  // Priors on non-centered parameters
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  
  // Priors on static deck preferences
  V1_pr ~ normal(0, 1);
  V2_pr ~ normal(0, 1);
  V3_pr ~ normal(0, 1);
  V4_pr ~ normal(0, 1);

  // Build boundary/tau vectors for trials
  vector[T] boundary_vec;
  vector[T] tau_vec;
  
  // First block
  boundary_vec[1:block] = rep_vector(boundary1, block);
  tau_vec[1:block]      = rep_vector(tau1, block);
  
  // Rest of trials (if T > block)
  if (T > block) {
    boundary_vec[(block+1):T] = rep_vector(boundary, T - block);
    tau_vec[(block+1):T]      = rep_vector(tau, T - block);
  }
  
  // Build deck values vector
  vector[4] deck_values;
  deck_values[1] = V1;
  deck_values[2] = V2;
  deck_values[3] = V3;
  deck_values[4] = V4;
  
  real sensitivity = pow(3, drift_con) - 1;
  
  // Likelihood
  target += igt_ard_model(choice, RT,
                       T,
                       sensitivity, deck_values,
                       boundary_vec, tau_vec, urgency,
                       win_indices_all, lose_indices_all);
}
