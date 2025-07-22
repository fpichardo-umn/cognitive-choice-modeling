functions {
  // Fitting EV DDM model with perseveration
  vector igt_model_lp(
			array[] int choice, array[] int shown, array[] real outcome,
			array[] real RT, vector ev, int Tsub, 
			real boundary, real tau, real beta,
			vector sensitivity, real gain, real loss,
			real retend, real k, real ep
			) {
    // Define values
    real         curUtil;      // Current utility
    int          curDeck;      // Current deck
    vector[Tsub] drift_rates;  // Current drift
    vector[4]    local_ev = ev;
    vector[4]    deck_pers = rep_vector(0.0, 4);
    real         gen_pers = 0.0;

    // Accumulate play/pass info
    array[Tsub] int play_indices;
    array[Tsub] int pass_indices;
    int play_count = 0;
    int pass_count = 0;

    // For each deck shown
    for (t in 1: Tsub) {
      // Deck presented to sub
      curDeck = shown[t];

      // Drift diffusion process - combine EV and perseveration
      drift_rates[t] = (local_ev[curDeck] * sensitivity[t]) + (deck_pers[curDeck] + gen_pers);

      // Compute utility
      curUtil = pow(abs(outcome[t]), gain) * (outcome[t] > 0 ? 1 : -1 * loss);
      
      // Decay perseveration values
      deck_pers = deck_pers * k;
      gen_pers = gen_pers * k;
      
      // Update general perseveration
      if (choice[t] == 1) {
        gen_pers += ep;
      } else {
        gen_pers -= ep;
      }
      
      // Update deck-specific perseveration only if played
      if (choice[t] == 1) {
        if (outcome[t] >= 0) {
          deck_pers[curDeck] += ep;
        } else {
          deck_pers[curDeck] -= ep;
        }
      }

      // Decay expected values
      local_ev *= retend;
      
      // Update selected deck
      local_ev[curDeck] += curUtil * choice[t]; // 0 update if not selected

      // Store indices for play and pass
      if (choice[t] == 1) {
        play_count += 1;
        play_indices[play_count] = t;
      } else {
        pass_count += 1;
        pass_indices[pass_count] = t;
      }
    }

    // Compute log probability for RTs/choice
    target += wiener_lpdf(
	RT[play_indices[:play_count]] | boundary, tau, beta, drift_rates[play_indices[:play_count]]);
    target += wiener_lpdf(
	RT[pass_indices[:pass_count]] | boundary, tau, 1-beta, -drift_rates[pass_indices[:pass_count]]);

    return local_ev;
  }
}

data {
  int<lower=1> 			 T; 	  // Number of trials
  real<lower=0> 		 minRT;   // Minimum RT + small value to restrict tau
  real 				 RTbound; // Lower bound or RT across all subjects (e.g., 0.1 second)
  array[T] real<lower=0> 	 RT;  	  // Reaction times
  array[T] int<lower=0, upper=1> choice;  // Binary choices made at each trial
  array[T] int<lower=0, upper=4> shown;   // Deck shown at each trial
  array[T] real 		 outcome; // Outcome at each trial
}

parameters {
  real<lower=-3, upper=3> boundary_pr;  // Boundary separation (a)
  real<lower=-3, upper=3> tau_pr;       // Non-decision time (tau)
  real<lower=-3, upper=3> beta_pr;      // Starting point
  real<lower=-3, upper=3> drift_con_pr; // Drift consistency parameter
  real<lower=-3, upper=3> gain_pr;      // Attention weight for rewards
  real<lower=-3, upper=3> loss_pr;      // Loss aversion
  real<lower=-3, upper=3> retend_pr;    // Retention rate
  real<lower=-3, upper=3> k_pr;         // Perseveration decay rate
  real<lower=-3, upper=3> ep_pr;        // Perseveration effect strength
}

transformed parameters {
  real<lower=0> 		   boundary;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0, upper=1> 	   beta;
  real<lower=-2, upper=2> 	   drift_con;
  real<lower=0, upper=1> 	   gain;
  real<lower=0, upper=1> 	   loss;
  real<lower=0, upper=1> 	   retend;
  real<lower=0, upper=1>  	   k;
  real  	   		   ep;

  boundary  = exp(boundary_pr);
  tau       = inv_logit(tau_pr) * (minRT - RTbound) * 0.99 + RTbound; // Ensures tau will always be at least 1% less than minRT
  beta      = inv_logit(beta_pr);
  drift_con = inv_logit(drift_con_pr) * 4 - 2;
  gain      = inv_logit(gain_pr);
  loss      = inv_logit(loss_pr);
  retend    = inv_logit(retend_pr);
  k         = inv_logit(k_pr);
  ep        = ep_pr;
}

model {
  // Priors
  boundary_pr  ~ normal(0, 1);
  tau_pr       ~ normal(0, 1);
  beta_pr      ~ normal(0, 1);
  drift_con_pr ~ normal(0, 1);
  gain_pr      ~ normal(0, 1);
  loss_pr      ~ normal(0, 1);
  retend_pr    ~ normal(0, 1);
  k_pr         ~ normal(0, 1);
  ep_pr        ~ normal(0, 1);

  // Initial subject-level deck expectations
  vector[4] ev = rep_vector(0., 4);

  // Initial trial data for theta
  vector[T] sensitivity = pow(to_vector(linspaced_array(T, 1, T)) / 10.0, drift_con);

  {
    ev = igt_model_lp(
			choice, shown, outcome,
			RT, ev, T,
			boundary, tau, beta,
			sensitivity, gain, loss,
			retend, k, ep
			);
  }
}
