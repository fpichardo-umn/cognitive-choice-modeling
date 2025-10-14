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

  // EV-ARD trial-level function (logic remains the same)
  real igt_ard_model(
      array[] int choice, array[] real wins, array[] real losses, array[] real RT,
      vector ev_init, int T,
      real sensitivity, real update, real wgt_pun, real wgt_rew,
      vector boundaries, vector taus, real urgency, real wd, real ws,
      array[,] int win_indices_all,
      array[,] int lose_indices_all,
      array[,] int other_indices) {

    vector[4] local_ev = ev_init;
    real log_lik = 0.0;

    real scaled_urgency = urgency * sensitivity;
    real scaled_wswd_plus = (ws + wd) * sensitivity;
    real scaled_wswd_minus = (ws - wd) * sensitivity;

    for (t in 1:T) {
      vector[12] drift_rates;
      int k = 1;
      for (i in 1:4) {
        drift_rates[k:k+2] = scaled_urgency +
                             scaled_wswd_plus * local_ev[i] +
                             scaled_wswd_minus * local_ev[other_indices[i]];
        k += 3;
      }

      if (RT[t] != 999) {
        log_lik += ard_win_all(RT[t], choice[t], taus[t], boundaries[t], drift_rates,
                               win_indices_all, lose_indices_all);
      }

      real curUtil = wgt_rew * wins[t] - wgt_pun * abs(losses[t]);
      local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
    }
    return log_lik;
  }
}

data {
  int<lower=1> sid;     		 // Subject ID
  int<lower=1> T;                         // Number of trials
  real<lower=0> minRT;                    // Minimum RT + small value to restrict tau
  real RTbound;                           // Lower bound of RT across all subjects
  array[T] int<lower=1, upper=4> choice;  // Choices made at each trial (1-4)
  array[T] real<lower=0> RT;              // Response times
  array[T] real<lower=0> wins;            // Win amount at each trial
  array[T] real<lower=0> losses;          // Loss amount at each trial
}

//---

transformed data {
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
  // ARD parameters
  real<lower=-5, upper=5> boundary1_pr; // Boundary separation (first 20 trials)
  real<lower=-5, upper=5> boundary_pr;  // Boundary separation (remaining trials)
  real<lower=-5, upper=5> tau1_pr;      // Non-decision time (first 20 trials)
  real<lower=-5, upper=5> tau_pr;       // Non-decision time (remaining trials)
  real urgency_pr;   // Urgency signal (V0)
  real wd_pr;        // Advantage weight
  real ws_pr;        // Sum weight
  
  // EV learning parameters
  real drift_con_pr;                    // Consistency parameter
  real wgt_pun_pr;                      // Weight for punishments
  real wgt_rew_pr;                      // Weight for rewards
  real update_pr;                       // Updating rate
}

transformed parameters {
  // ARD parameters
  real<lower=0, upper=6> boundary1;
  real<lower=0, upper=6> boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0> urgency;
  real<lower=0> wd;
  real<lower=0> ws;
  
  // EV parameters
  real<lower=0, upper=5> drift_con;
  real<lower=0, upper=1> wgt_pun;
  real<lower=0, upper=1> wgt_rew;
  real<lower=0, upper=1> update;
  
  boundary1 = inv_logit(boundary1_pr) * 4.99 + 0.01;
  boundary  = inv_logit(boundary_pr) * 4.99 + 0.01;
  tau1 = inv_logit(tau1_pr) * (minRT - RTbound - 1e-6) * 0.95 + RTbound;
  tau = inv_logit(tau_pr) * (minRT - RTbound - 1e-6) * 0.95 + RTbound;
  urgency = log1p_exp(urgency_pr);
  wd = log1p_exp(wd_pr);
  ws = log1p_exp(ws_pr);
  drift_con = inv_logit(drift_con_pr) * 5;
  wgt_pun = inv_logit(wgt_pun_pr);
  wgt_rew = inv_logit(wgt_rew_pr);
  update = inv_logit(update_pr);
}

model {
  // Priors
  boundary1_pr ~ normal(0, 1);
  boundary_pr ~ normal(0, 1);
  tau1_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  urgency_pr ~ normal(0, 1);
  wd_pr ~ normal(0, 1);
  ws_pr ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  wgt_pun_pr ~ normal(0, 1);
  wgt_rew_pr ~ normal(0, 1);
  update_pr ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev = rep_vector(0.0, 4);
  real sensitivity = pow(3, drift_con) - 1;
  
  // Create trial-varying boundary and tau vectors
  vector[T] boundaries = append_row(rep_vector(boundary1, 20), rep_vector(boundary, T - 20));
  vector[T] taus = append_row(rep_vector(tau1, 20), rep_vector(tau, T - 20));
  
  // Combined model computation
  target += igt_ard_model(choice, wins, losses, RT, ev, T, 
                             sensitivity, update, wgt_pun, wgt_rew,
                             boundaries, taus, urgency, wd, ws,
                             win_indices_all, lose_indices_all, other_indices);
}
