data {
  int<lower=1> 		     sid;     // Subject ID
  int<lower=1> 		     T;       // Number of trials
  real<lower=0> 	     minRT;   // Minimum RT + small value to restrict tau
  real<lower=0>		     RTbound; // Lower bound or RT across all subjects (e.g., 0.1 second)
  array[T] real<lower=0>     RT;      // Reaction times
  array[T] int<lower=0, upper=1> choice;  // Binary choices (1=play, 0=pass)
}

parameters {
  real<lower=-5, upper=5> boundary1_pr; // Boundary separation (T1 a)
  real<lower=-5, upper=5> boundary_pr;  // Boundary separation (a)
  real<lower=-3, upper=3> tau1_pr;      // Non-decision time (T1 tau)
  real<lower=-3, upper=3> tau_pr;       // Non-decision time (tau)
  real<lower=-3, upper=3> beta_pr;      // Starting point
  real<lower=-3, upper=3> drift;        // Drift rate
}

transformed parameters {
  real<lower=0, upper=6> 	   boundary1;
  real<lower=0, upper=6> 	   boundary;
  real<lower=RTbound, upper=minRT> tau1;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0, upper=1> 	   beta;

  boundary1 = inv_logit(boundary1_pr) * 5 + 0.01;
  boundary  = inv_logit(boundary_pr) * 5 + 0.01;
  tau1      = inv_logit(tau1_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound;
  tau       = inv_logit(tau_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound;
  beta      = inv_logit(beta_pr);
}

model {
  // Priors
  boundary1_pr ~ normal(0, 1);
  boundary_pr  ~ normal(0, 1);
  tau1_pr      ~ normal(0, 1);
  tau_pr       ~ normal(0, 1);
  beta_pr      ~ normal(0, 1);
  drift	       ~ normal(0, 1);
  
  // Accumulate play/pass info with valid RTs
  array[T] int play_indices;
  array[T] int pass_indices;
  int play_count = 0;
  int pass_count = 0;
  
  // Separate by block and choice
  vector[T] boundaries;
  vector[T] nondt;
  
  for (t in 1:T) {
    // Set block-specific boundary and tau
    if (t <= 20) {
      boundaries[t] = boundary1;
      nondt[t] = tau1;
    } else {
      boundaries[t] = boundary;
      nondt[t] = tau;
    }
    
    // Store indices for play and pass - ONLY for valid RTs
    if (RT[t] != 999) {
      if (choice[t] == 1) {
        play_count += 1;
        play_indices[play_count] = t;
      } else {
        pass_count += 1;
        pass_indices[pass_count] = t;
      }
    }
  }
  
  // Likelihood
  RT[play_indices[:play_count]] ~ wiener(boundaries[play_indices[:play_count]], nondt[play_indices[:play_count]], beta, drift);
  RT[pass_indices[:pass_count]] ~ wiener(boundaries[pass_indices[:pass_count]], nondt[pass_indices[:pass_count]], 1-beta, -drift);
}
