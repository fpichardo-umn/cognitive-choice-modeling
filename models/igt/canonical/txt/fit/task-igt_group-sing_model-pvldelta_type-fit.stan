functions {
  vector igt_model_lp(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, int Tsub,
        real sensitivity, real gain, real loss, real update
        ) {
    // Define values
    real curUtil;
    vector[4] local_ev = ev;
    
    // For each trial
    for (t in 1:Tsub) {
      // Choice probability
      target += categorical_logit_lpmf(choice[t] | sensitivity * local_ev);
      
      // Compute utility from wins and losses separately
      real win_component = (wins[t] == 0) ? 0.0 : exp(gain * log(wins[t]));
real loss_component = (losses[t] == 0) ? 0.0 : exp(gain * log(losses[t]));
curUtil = win_component - loss * loss_component;
      
      // Update expected value using delta rule
      local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
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
  real gain_pr;    // Value sensitivity
  real loss_pr;    // Loss aversion
  real update_pr;  // Updating rate
}

transformed parameters {
  real<lower=0, upper=5> con;
  real<lower=0, upper=2> gain;
  real<lower=0, upper=10> loss;
  real<lower=0, upper=1> update;
  
  con = inv_logit(con_pr) * 5;
  gain = inv_logit(gain_pr) * 2;
  loss = inv_logit(loss_pr) * 10;
  update = inv_logit(update_pr);
}

model {
  // Priors
  con_pr ~ normal(0, 1);
  gain_pr ~ normal(0, 1);
  loss_pr ~ normal(0, 1);
  update_pr ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev = rep_vector(0., 4);
  real sensitivity = expm1(log(3) * con);
  
  // Run model
  ev = igt_model_lp(choice, wins, losses, ev, T, sensitivity, gain, loss, update);
}
