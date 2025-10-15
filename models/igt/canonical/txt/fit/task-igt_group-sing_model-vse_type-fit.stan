functions {
  vector igt_model_lp(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev_exploit, vector ev_explore, int Tsub,
        real sensitivity, real gain, real loss, real decay, 
        real explore_alpha, real explore_bonus
        ) {
    // Define values
    real curUtil;
    vector[4] local_ev_exploit = ev_exploit;
    vector[4] local_ev_explore = ev_explore;
    vector[4] combined_value;
    
    // For each trial
    for (t in 1:Tsub) {
      // Combine exploitation and exploration values
      combined_value = local_ev_exploit + local_ev_explore;
      
      // Choice probability using sensitivity
      target += categorical_logit_lpmf(choice[t] | sensitivity * combined_value);
      
      // Calculate utility (using value sensitivity for wins and losses)
      real win_component = (wins[t] == 0) ? 0.0 : exp(gain * log(wins[t]));
real loss_component = (losses[t] == 0) ? 0.0 : exp(gain * log(losses[t]));
curUtil = win_component - loss * loss_component;
      
      // Exploitation: Decay all deck values
      local_ev_exploit = local_ev_exploit * (1 - decay);
      
      // Exploitation: Update chosen deck
      local_ev_exploit[choice[t]] += curUtil;
      
      // Exploration: Reset chosen deck to zero
      local_ev_explore[choice[t]] = 0;
      
      // Exploration: Update unchosen decks (return to exploration bonus)
      for (d in 1:4) {
        if (d != choice[t]) {
          local_ev_explore[d] += explore_alpha * (explore_bonus - local_ev_explore[d]);
        }
      }
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
  real decay_pr;         // Decay parameter for exploitation
  real explore_alpha_pr; // Learning rate for exploration
  real explore_bonus_pr; // Exploration bonus parameter
}

transformed parameters {
  real<lower=0, upper=5> con;
  real<lower=0, upper=1> gain;
  real<lower=0, upper=10> loss;
  real<lower=0, upper=1> decay;
  real<lower=0, upper=1> explore_alpha;
  real<lower=-10, upper=10> explore_bonus;
  
  con = inv_logit(con_pr) * 5;
  gain = inv_logit(gain_pr);
  loss = inv_logit(con_pr) * 10;
  decay = inv_logit(decay_pr);
  explore_alpha = inv_logit(explore_alpha_pr);
  explore_bonus = -10 + inv_logit(explore_bonus_pr) * 20;
}

model {
  // Priors
  con_pr ~ normal(0, 1);
  gain_pr ~ normal(0, 1);
  loss_pr ~ normal(0, 1);
  decay_pr ~ normal(0, 1);
  explore_alpha_pr ~ normal(0, 1);
  explore_bonus_pr ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev_exploit = rep_vector(0., 4);
  vector[4] ev_explore = rep_vector(0., 4);
  real sensitivity = expm1(log(3) * con);
  
  // Run model
  ev_exploit = igt_model_lp(choice, wins, losses, 
				ev_exploit, ev_explore, T, 
                           	sensitivity, gain, loss, decay, 
				explore_alpha, explore_bonus);
}
