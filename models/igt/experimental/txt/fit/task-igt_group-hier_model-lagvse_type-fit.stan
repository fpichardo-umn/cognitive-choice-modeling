// Optimized Hierarchical lagVSE Model for the Iowa Gambling Task
functions {
  // --- Core function to calculate log-likelihood for a SINGLE subject ---
  // This function is now self-contained and returns a 'real' value.
  real igt_subject(
        array[] int choice, array[] real wins, array[] real losses,
        int Tsub, real sensitivity, real gain, real loss, real decay, real phi
        ) {
    real log_lik = 0.0;
    vector[4] local_ev_exploit = rep_vector(0.0, 4);
    vector[4] local_choice_lag = rep_vector(0.0, 4);
    vector[4] combined_value;
    real curUtil;

    for (t in 1:Tsub) {
      // Increment lag for all options
      local_choice_lag += 1;
      
      // Combine exploitation and exploration values
      combined_value = local_ev_exploit + phi * local_choice_lag;
      
      // Calculate choice probability and add to log-likelihood
      // This now adds to a local variable instead of the global 'target'
      log_lik += categorical_logit_lpmf(choice[t] | sensitivity * combined_value);
      
      // Calculate utility
      real win_component = (wins[t] == 0) ? 0.0 : exp(gain * log(wins[t]));
      real loss_component = (losses[t] == 0) ? 0.0 : exp(gain * log(losses[t]));
      curUtil = win_component - loss * loss_component;
      
      // Exploitation: Decay all deck values
      local_ev_exploit *= (1 - decay);
      
      // Exploitation: Update chosen deck
      local_ev_exploit[choice[t]] += curUtil;
      
      // Exploration: Reset lag for chosen deck
      local_choice_lag[choice[t]] = 0;
    }
    
    return log_lik;
  }
  
  // --- Partial sum function for reduce_sum parallelization ---
  // This function loops over a "slice" of subjects.
  real partial_sum_func(array[] int slice_n, int start, int end,
                        // Data arrays
                        int N, array[,] int choice, array[,] real wins, array[,] real losses, array[] int Tsubj,
                        // Transformed Parameter arrays
                        array[] real con, array[] real gain, array[] real loss, array[] real decay, array[] real phi
                       ) {
    real log_lik = 0.0;

    for (n in start:end) {
      // Calculate sensitivity for this subject
      real sensitivity = expm1(log(3) * con[n]);
      
      // Call the single-subject function and accumulate the log-likelihood
      log_lik += igt_subject(choice[n, 1:Tsubj[n]],
                                  wins[n, 1:Tsubj[n]],
                                  losses[n, 1:Tsubj[n]],
                                  Tsubj[n],
                                  sensitivity, gain[n], loss[n], decay[n], phi[n]);
    }
    return log_lik;
  }
}

data {
  int<lower=1> N; // Number of subjects
  int<lower=1> T; // Maximum number of trials
  array[N] int<lower=1> sid; // Subject IDs
  array[N] int<lower=1> Tsubj; // Number of trials for each subject
  array[N, T] int<lower=1, upper=4> choice; // Choices made at each trial (1-4)
  array[N, T] real<lower=0> wins; // Win amount at each trial
  array[N, T] real<lower=0> losses; // Loss amount at each trial
}

// --- prepare data for parallelization ---
transformed data {
  array[N] int subject_indices;
  for (i in 1:N) {
    subject_indices[i] = i;
  }
}

parameters {
  // Group-level hyperparameters
  array[5] real mu_pr;
  array[5] real<lower=0> sigma;

  // Subject-level raw parameters
  array[N] real con_pr;
  array[N] real gain_pr;
  array[N] real loss_pr;
  array[N] real decay_pr;
  array[N] real phi_pr;
}

transformed parameters {
  // Transform subject-level raw parameters
  array[N] real<lower=0, upper=3> con;
  array[N] real<lower=0, upper=1> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> decay;
  array[N] real<lower=-1, upper=1> phi;

  // Hierarchical transformation
  con   = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(con_pr)) * 3);
  gain  = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(gain_pr)));
  loss  = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(loss_pr)) * 10);
  decay = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(decay_pr)));
  phi   = to_array_1d(-1 + inv_logit(mu_pr[5] + sigma[5] .* to_vector(phi_pr)) * 2);
}

// --- OPTIMIZED model block ---
model {
  // Hyperpriors (vectorized for efficiency)
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  // Subject-level priors (vectorized for efficiency)
  con_pr   ~ normal(0, 1);
  gain_pr  ~ normal(0, 1);
  loss_pr  ~ normal(0, 1);
  decay_pr ~ normal(0, 1);
  phi_pr   ~ normal(0, 1);

  int grainsize = max(1, N %/% 4);
  
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       // Pass data
                       N, choice, wins, losses, Tsubj,
                       // Pass transformed parameters
                       con, gain, loss, decay, phi);
}

generated quantities {
  // Group-level parameters in interpretable scale
  real<lower=0, upper=5> mu_con;
  real<lower=0, upper=1> mu_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0, upper=1> mu_decay;
  real<lower=-1, upper=1> mu_phi;
  
  mu_con   = inv_logit(mu_pr[1]) * 3;
  mu_gain  = inv_logit(mu_pr[2]);
  mu_loss  = inv_logit(mu_pr[3]) * 10;
  mu_decay = inv_logit(mu_pr[4]);
  mu_phi   = -1 + inv_logit(mu_pr[5]) * 2;
}
