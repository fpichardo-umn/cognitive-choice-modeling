functions {
  vector igt_model_lp(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev_exploit, vector choice_lag, int Tsub,
        real sensitivity, real gain, real loss,
	real update, real decay, real phi
        ) {
    // Define values
    real curUtil;
    vector[4] local_ev_exploit = ev_exploit;
    vector[4] local_choice_lag = choice_lag;  // Track trials since last chosen
    vector[4] combined_value;
    
    // For each trial
    for (t in 1:Tsub) {
      // Increment lag for all options
      local_choice_lag += 1;

      // Combine exploitation and exploration values
      combined_value = local_ev_exploit + phi * local_choice_lag;
      
      // Choice probability using sensitivity
      target += categorical_logit_lpmf(choice[t] | sensitivity * combined_value);
      
      // Calculate utility (using value sensitivity for wins and losses)
      real win_component = (wins[t] == 0) ? 0.0 : exp(gain * log(wins[t]));
real loss_component = (losses[t] == 0) ? 0.0 : exp(gain * log(losses[t]));
curUtil = win_component - loss * loss_component;
      
      // Exploitation: Decay all deck values
      local_ev_exploit = local_ev_exploit * (1 - decay);
      
      // Exploitation: Update chosen deck
      local_ev_exploit[choice[t]] += curUtil + update * (curUtil - local_ev_exploit[choice[t]]);;
      
      // Exploration: Reset lag for chosen deck
      local_choice_lag[choice[t]] = 0;
    }
    
    return local_ev_exploit;
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
  real con_pr;           // Consistency parameter
  real gain_pr;          // Value sensitivity parameter
  real loss_pr;          // Loss aversion
  real update_pr;        // Updating rate
  real decay_pr;         // Decay parameter for exploitation
  real phi_pr;		 // Choice lag weight
}

transformed parameters {
  real<lower=0, upper=5> con;
  real<lower=0, upper=1> gain;
  real<lower=0, upper=10> loss;
  real<lower=0, upper=1> update;
  real<lower=0, upper=1> decay;
  real<lower=-1, upper=1> phi;
  
  con = inv_logit(con_pr) * 5;
  gain = inv_logit(gain_pr);
  loss = inv_logit(con_pr) * 10;
  update = inv_logit(update_pr);
  decay = inv_logit(decay_pr);
  phi = -1 + inv_logit(phi_pr) * 2;
}

model {
  // Priors
  con_pr ~ normal(0, 1);
  gain_pr ~ normal(0, 1);
  loss_pr ~ normal(0, 1);
  update_pr ~ normal(0, 1);
  decay_pr ~ normal(0, 1);
  phi_pr ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev_exploit = rep_vector(0., 4);
  vector[4] choice_lag = rep_vector(0., 4);
  real sensitivity = expm1(log(3) * con);
  
  // Run model
  ev_exploit = igt_model_lp(choice, wins, losses,
				ev_exploit, choice_lag, T, 
                           	sensitivity, gain, loss,
				update, decay, phi);
}

