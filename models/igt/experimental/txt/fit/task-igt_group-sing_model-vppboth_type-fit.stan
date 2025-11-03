functions {
  vector igt_model_lp(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, vector pers, int Tsub,
        real sensitivity, real update, real gain, real loss, 
        real epP, real epN, real K, real w, real decay
        ) {
    // Define values
    real curUtil;
    vector[4] local_ev = ev;
    vector[4] local_pers = pers;
    vector[4] V;

    // For each trial
    for (t in 1:Tsub) {
      // Calculate combined value
      V = w * local_ev + (1-w) * local_pers;
      
      // Choice probability using sensitivity
      target += categorical_logit_lpmf(choice[t] | sensitivity * V);
      
      // Decay perseverance
      local_pers = local_pers * K;
      
      // Calculate utility & update perseverance
      // Use separate values for wins and losses
      real win_component = (wins[t] == 0) ? 0.0 : exp(gain * log(wins[t]));
real loss_component = (losses[t] == 0) ? 0.0 : exp(gain * log(losses[t]));
curUtil = win_component - loss * loss_component;
      
      // Update perseverance based on net outcome
      if (wins[t] >= losses[t]) {
        local_pers[choice[t]] += epP;
      } else {
        local_pers[choice[t]] += epN;
      }
      
      // First decay all deck values
      local_ev = local_ev * (1 - decay);
      
      // Then update chosen deck with delta rule
      local_ev[choice[t]] += curUtil + update * (curUtil - local_ev[choice[t]]);
    }
    
    return local_ev;
  }
}

data {
  int<lower=1> sid;     		 // Subject ID
  int<lower=1> T;                        // Number of trials
  array[T] int<lower=1, upper=4> choice; // Choices made at each trial
  array[T] real<lower=0> wins;           // Win amount at each trial
  array[T] real<lower=0> losses;         // Loss amount at each trial (positive values)
}

parameters {
  real con_pr;     // Consistency parameter
  real update_pr;  // Learning rate parameter
  real gain_pr;    // Outcome sensitivity parameter
  real loss_pr;    // Loss aversion parameter
  real epP;     // Positive perseverance strength
  real epN;     // Negative perseverance strength
  real K_pr;       // Perseverance decay parameter
  real w_pr;       // Value weight parameter
  real decay_pr;   // Decay rate
}

transformed parameters {
  real<lower=0, upper=5> con;
  real<lower=0, upper=1> update;
  real<lower=0, upper=1> gain;
  real<lower=0, upper=10> loss;
  real<lower=0, upper=1> K;
  real<lower=0, upper=1> w;
  real<lower=0, upper=1> decay;
  
  con = inv_logit(con_pr) * 5;
  update = inv_logit(update_pr);
  gain = inv_logit(gain_pr);
  loss = inv_logit(loss_pr) * 10;
  K = inv_logit(K_pr);
  w = inv_logit(w_pr);
  decay = inv_logit(decay_pr);
}

model {
  // Priors
  con_pr ~ normal(0, 1);
  update_pr ~ normal(0, 1);
  gain_pr ~ normal(0, 1);
  loss_pr ~ normal(0, 1);
  epP ~ normal(0, 1);
  epN ~ normal(0, 1);
  K_pr ~ normal(0, 1);
  w_pr ~ normal(0, 1);
  decay_pr ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev = rep_vector(0., 4);
  vector[4] pers = rep_vector(0., 4);
  real sensitivity = expm1(log(3) * con);
  
  // Run model
  ev = igt_model_lp(choice, wins, losses, ev, pers, T, sensitivity, update, gain, loss, epP, epN, K, w, decay);
}
