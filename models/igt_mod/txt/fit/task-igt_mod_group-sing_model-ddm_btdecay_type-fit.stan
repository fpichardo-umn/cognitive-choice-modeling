functions {
  // Fitting DDM Decay model
  void ddm_decay_model_lp(
    array[] int choice, array[] real RT, int T,
    real boundary_init, real boundary_decay,
    real tau_init, real tau_decay, real RTbound,
    real beta, real drift
			) {

    // Accumulate info
    array[Tsub] int play_indices;
    array[Tsub] int pass_indices;
    array[Tsub] real boundary;
    array[Tsub] real tau;
    int play_count = 0;
    int pass_count = 0;

    // For each deck shown
    for (t in 1: Tsub) {
      // Store indices for play and pass
      if (choice[t] == 1) {
        play_count += 1;
        play_indices[play_count] = t;
	boundary[i] = boundary_init * exp(-boundary_decay * t);
	tau[i] = RTbound + (tau_init - RTbound) * exp(-tau_decay * t);

      } else {
        pass_count += 1;
        pass_indices[pass_count] = t;
	boundary[i] = boundary_init * exp(-boundary_decay * t);
	tau[i] = RTbound + (tau_init - RTbound) * exp(-tau_decay * t);

      }
    }

    // Compute log probability for RTs/choice
    target += wiener_lpdf(
	RT[play_indices[:play_count]] | boundary[play_indices[:play_count]], tau[play_indices[:play_count]], beta, drift_rates[play_indices[:play_count]]);
    target += wiener_lpdf(
	RT[pass_indices[:pass_count]] | boundary[pass_indices[:pass_count]], tau[pass_indices[:pass_count]], 1-beta, -drift_rates[pass_indices[:pass_count]]);
  }
}

data {
  int<lower=1> T;                     // Number of trials
  real<lower=0> minRT;                // Minimum RT + small value to restrict tau
  real RTbound;                       // Lower bound or RT across all subjects (e.g., 0.1 second)
  array[T] real<lower=0> RT;          // Reaction times
  array[T] int<lower=0, upper=1> choice;  // Binary choices made at each trial
  array[T] int<lower=0, upper=4> shown;   // Deck shown at each trial
  array[T] real outcome;              // Outcome at each trial
}

parameters {
  real<lower=-3, upper=3> boundary_init_pr;  // Initial boundary separation
  real<lower=-3, upper=3> boundary_decay_pr; // Boundary decay rate
  real<lower=-3, upper=3> tau_init_pr;       // Initial non-decision time
  real<lower=-3, upper=3> tau_decay_pr;      // Tau decay rate
  real<lower=-3, upper=3> beta_pr;           // Starting point
  real<lower=-3, upper=3> drift;             // Drift rate
}

transformed parameters {
  real<lower=0> boundary_init;
  real<lower=0, upper=1> boundary_decay;
  real<lower=RTbound, upper=minRT> tau_init;
  real<lower=0, upper=1> tau_decay;
  real<lower=0, upper=1> beta;

  boundary_init = exp(inv_logit(boundary_init_pr) * 10 - 5);
  boundary_decay = inv_logit(boundary_decay_pr);
  tau_init = inv_logit(tau_init_pr) * (minRT - RTbound) * 0.99 + RTbound;
  tau_decay = inv_logit(tau_decay_pr);
  beta = inv_logit(beta_pr);
}

model {
  // Priors
  boundary_init_pr ~ normal(0, 1);
  boundary_decay_pr ~ normal(0, 1);
  tau_init_pr ~ normal(0, 1);
  tau_decay_pr ~ normal(0, 1);
  beta_pr ~ normal(0, 1);
  drift ~ normal(0, 1);
  
  // Apply the DDM decay model
  ddm_decay_model(
    choice, RT, T,
    boundary_init, boundary_decay,
    tau_init, tau_decay, RTbound,
    beta, drift
  );
}