functions {
  vector igt_model_lp(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev, int Tsub, 
        real sensitivity, real update, real wgt_pun, real wgt_rew
        ) {
    // Define values
    real curUtil;
    vector[4] local_ev = ev;
    
    // For each trial
    for (t in 1:Tsub) {
      // Choice probability using consistency
      target += categorical_logit_lpmf(choice[t] | sensitivity * local_ev);
      
      // Compute utility with separate weights for rewards and punishments
      curUtil = wgt_rew * wins[t] + wgt_pun * losses[t];
      
      // Update expected value of chosen deck
      local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
    }
    
    return local_ev;
  }
}

data {
  int<lower=1> sid;     		 // Subject ID
  int<lower=1> T;                        // Number of trials
  array[T] int<lower=1, upper=4> choice; // Choices made at each trial (1-4)
  array[T] real<lower=0> wins;           // Win amount at each trial
  array[T] real<lower=0> losses;         // Loss amount at each trial (positive values)
}

parameters {
  real con_pr;      // Consistency parameter
  real wgt_pun_pr;  // Weight for punishments (negative)
  real wgt_rew_pr;  // Weight for rewards (positive)
  real update_pr;   // Updating rate
}

transformed parameters {
  real<lower=0, upper=5> con;
  real<lower=0, upper=1> wgt_pun;
  real<lower=0, upper=1> wgt_rew;
  real<lower=0, upper=1> update;
  
  con     = inv_logit(con_pr) * 5;
  wgt_pun = inv_logit(wgt_pun_pr);
  wgt_rew = inv_logit(wgt_rew_pr);
  update  = inv_logit(update_pr);
}

model {
  // Priors
  con_pr     ~ normal(0, 1);
  wgt_pun_pr ~ normal(0, 1);
  wgt_rew_pr ~ normal(0, 1);
  update_pr  ~ normal(0, 1);
  
  // Initialize values
  vector[4] ev     = rep_vector(0., 4);
  real sensitivity = expm1(log(3) * con);
  
  // Run model
  ev = igt_model_lp(choice, wins, losses, ev, T, 
			sensitivity, update, wgt_pun, wgt_rew);
}
