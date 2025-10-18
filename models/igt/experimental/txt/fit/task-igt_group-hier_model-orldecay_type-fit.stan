functions {
  // Returns log-likelihood instead of modifying target
  real igt_subject(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, vector ef, int Tsub,
        real Drew, real Dpun, real K, real betaF, real betaP
        ) {
    real log_lik = 0.0;  // Local accumulator
    vector[4] local_ev = ev;
    vector[4] local_ef = ef;
    vector[4] pers = rep_vector(0.0, 4);
    vector[4] util;
    real PEval;
    real PEfreq;
    vector[4] PEfreq_fic;
    array[Tsub] real sign_outcome;
    real K_tr = pow(3, K) - 1;
    
    // Calculate sign for each trial
    for (t in 1:Tsub) {
      sign_outcome[t] = wins[t] >= losses[t] ? 1.0 : -1.0;
    }

    // For each trial
    for (t in 1:Tsub) {
      util = local_ev + local_ef * betaF + pers * betaP;
      
      // Accumulate to local variable instead of target
      log_lik += categorical_logit_lpmf(choice[t] | util);
      
      PEval = wins[t] - losses[t] - local_ev[choice[t]];
      PEfreq = sign_outcome[t] - local_ef[choice[t]];
      
      for (d in 1:4) {
        PEfreq_fic[d] = -sign_outcome[t]/3.0 - local_ef[d];
      }
      
      if (wins[t] >= losses[t]) {
        local_ef += PEfreq_fic * (1 - Dpun);

	// Decay all
        local_ef = local_ef * (1 - Drew);
        local_ev = local_ev * (1 - Drew);

        local_ef[choice[t]] = PEfreq;
        local_ev[choice[t]] = PEval;
      } else {
        local_ef += PEfreq_fic * (1 - Drew);

	// Decay all
        local_ef = local_ef * (1 - Dpun);
        local_ev = local_ev * (1 - Dpun);

        local_ef[choice[t]] = PEfreq;
        local_ev[choice[t]] = PEval;
      }
      
      pers[choice[t]] = 1;
      pers = pers / (1 + K_tr);
    }
    
    return log_lik;  // Return log-likelihood, not state
  }
  
  // Parallelization wrapper
  real partial_sum_func(array[] int slice_n, int start, int end,
                        // Data
                        array[,] int choice, array[,] real wins, array[,] real losses, 
                        array[] int Tsubj,
                        // Transformed parameters
                        array[] real Drew, array[] real Dpun, array[] real K,
                        array[] real betaF, array[] real betaP) {
    real log_lik = 0.0;
    
    for (n in start:end) {
      vector[4] ev = rep_vector(0., 4);
      vector[4] ef = rep_vector(0., 4);
      
      log_lik += igt_subject(choice[n, 1:Tsubj[n]],
                             wins[n, 1:Tsubj[n]], 
                             losses[n, 1:Tsubj[n]],
                             ev, ef, Tsubj[n], 
                             Drew[n], Dpun[n], K[n], 
                             betaF[n], betaP[n]);
    }
    return log_lik;
  }
}

data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1> sid;
  array[N] int<lower=1> Tsubj;
  array[N, T] int<lower=1, upper=4> choice;
  array[N, T] real<lower=0> wins;
  array[N, T] real<lower=0> losses;
}

transformed data {
  // Create subject indices with explicit loop
  array[N] int subject_indices;
  for (i in 1:N) {
    subject_indices[i] = i;
  }
}

parameters {
  array[5] real mu_pr;
  array[5] real<lower=0> sigma;

  array[N] real Drew_pr;
  array[N] real Dpun_pr;
  array[N] real K_pr;
  array[N] real betaF_pr;
  array[N] real betaP_pr;
}

transformed parameters {
  array[N] real<lower=0, upper=1> Drew;
  array[N] real<lower=0, upper=1> Dpun;
  array[N] real<lower=0, upper=5> K;
  array[N] real betaF;
  array[N] real betaP;
  
  // Vectorized transformation
  Drew  = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(Drew_pr)));
  Dpun  = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(Dpun_pr)));
  K     = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(K_pr)) * 5);
  betaF = to_array_1d(mu_pr[4] + sigma[4] .* to_vector(betaF_pr));
  betaP = to_array_1d(mu_pr[5] + sigma[5] .* to_vector(betaP_pr));
}

model {
  // Vectorized priors
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  Drew_pr  ~ normal(0, 1);
  Dpun_pr  ~ normal(0, 1);
  K_pr     ~ normal(0, 1);
  betaF_pr ~ normal(0, 1);
  betaP_pr ~ normal(0, 1);

  // Parallel processing
  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       Drew, Dpun, K, betaF, betaP);
}

generated quantities {
  real<lower=0, upper=1> mu_Drew;
  real<lower=0, upper=1> mu_Dpun;
  real<lower=0, upper=5> mu_K;
  real mu_betaF;
  real mu_betaP;
  
  mu_Drew  = inv_logit(mu_pr[1]);
  mu_Dpun  = inv_logit(mu_pr[2]);
  mu_K     = inv_logit(mu_pr[3]) * 5;
  mu_betaF = mu_pr[4];
  mu_betaP = mu_pr[5];
}
