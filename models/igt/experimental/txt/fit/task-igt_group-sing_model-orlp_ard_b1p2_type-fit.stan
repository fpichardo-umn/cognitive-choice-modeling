// Single-Subject ORL-ARD Model for the Iowa Gambling Task
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
    return fmax(log_pdf, -50);
// fmax is differentiable
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
      real Arew, real Apun, real K,
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

    real wswd_plus = (ws + wd);
real wswd_minus = (ws - wd);

    for (t in 1:T) {
      sign_outcome[t] = wins[t] >= losses[t] ?
1.0 : -1.0;
    }

    for (t in 1:T) {
      vector[12] drift_rates;
      vector[4] value_utility; // Store value utility for all 4 decks
      int k = 1;

      // 1. Calculate value-based utility (EV + EF) for all decks first
      for (i in 1:4) {
        value_utility[i] = local_ev[i] + local_ef[i] * betaF;
      }
      
      // 2. Calculate drift rates
      for (i in 1:4) {
        // Get the value utility of the 3 other decks
        vector[3] other_utilities = value_utility[other_indices[i]];
        
        // 3. Calculate advantage-scaled drift based ONLY on value-utility
        //    (urgency + (ws+wd)*V_i + (ws-wd)*V_others)
        vector[3] advantage_drift = rep_vector(urgency, 3) +
                                    wswd_plus * value_utility[i] +
                                    wswd_minus * other_utilities;

        // 4. Add the perseverance bias *after* the advantage scaling
        drift_rates[k:k+2] = advantage_drift + pers[i] * betaP;
        
        k += 3;
      }

      // --- [ End of modification ] ---

      if (RT[t] != 999) {
        log_lik += ard_win_all(RT[t], choice[t], taus[t], boundaries[t], drift_rates,
                               win_indices_all, lose_indices_all);
}

      // --- [ ORL Learning Updates (unchanged) ] ---
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
}

data {
  int<lower=1> T;                           // Total number of trials
  real<lower=0> minRT;                      // Minimum RT for this subject
  real<lower=0> RTbound;                    // RT bound
  array[T] int<lower=1, upper=4> choice;   // Choices
  array[T] real RT;                         // Response times
  array[T] real wins;                       // Win amounts
  array[T] real losses;                     // Loss amounts
}

transformed data {
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
  // Individual-level parameters (non-centered)
  real boundary1_pr;
  real boundary_pr;
  real tau1_pr;
  real tau_pr;
  real urgency_pr;
  real wd_pr;
  real ws_pr;

  real Arew_pr;
  real Apun_pr;
  real K_pr;
  real betaF_pr;
  real betaP_pr;
}

transformed parameters {
  // Transformed individual-level parameters
  real<lower=0.001, upper=5> boundary1;
  real<lower=0.001, upper=5> boundary;
  real<lower=0> tau1;
  real<lower=0> tau;
  real<lower=0.001, upper=20> urgency;
  real<lower=0.001, upper=10> wd;
  real<lower=0.001, upper=10> ws;

  real<lower=0, upper=1> Arew;
  real<lower=0, upper=1> Apun;
  real<lower=0, upper=5> K;
  real betaF;
  real betaP;
  
  // Transform parameters
  boundary1 = inv_logit(boundary1_pr) * 4.99 + 0.001;
  boundary  = inv_logit(boundary_pr) * 4.99 + 0.001;
  tau1 = inv_logit(tau1_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;
  tau  = inv_logit(tau_pr) * (minRT - RTbound - 0.02) * 0.95 + RTbound;

  urgency = inv_logit(urgency_pr) * 19.999 + 0.001;
  wd = inv_logit(wd_pr) * 9.999 + 0.001;
  ws = inv_logit(ws_pr) * 9.999 + 0.001;

  Arew      = inv_logit(Arew_pr);
  Apun      = inv_logit(Apun_pr);
  K         = inv_logit(K_pr) * 5;
  betaF     = betaF_pr;
  betaP     = betaP_pr;
}

model {
  // Priors on non-centered parameters
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  wd_pr ~ normal(0, 1);
  ws_pr ~ normal(0, 1);

  Arew_pr ~ normal(0, 1);
  Apun_pr ~ normal(0, 1);
  K_pr ~ normal(0, 1);
  betaF_pr ~ normal(0, 1);
  betaP_pr ~ normal(0, 1);
  
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
  
  // Likelihood with ORL learning
  vector[4] ev = rep_vector(0.0, 4);
  vector[4] ef = rep_vector(0.0, 4);
  
  target += igt_ard_model(choice, wins, losses, RT,
                       ev, ef, T,
                       Arew, Apun, K,
                       betaF, betaP,
                       boundary_vec, tau_vec, urgency, wd, ws,
                       win_indices_all, lose_indices_all, other_indices);
}
