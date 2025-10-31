
data {
  int<lower=1> 		     sid;     // Subject ID
  int<lower=1> 		     T;       // Number of trials
  int<lower=0>  	     Nplay;   // Number of play trials
  int<lower=0>  	     Npass;   // Number of pass trials
  real<lower=0> 	     minRT;   // Minimum RT + small value to restrict tau
  real<lower=0>		     RTbound; // Lower bound or RT across all subjects (e.g., 0.1 second)
  array[T] real<lower=0> RTplay;  // Reaction times for play trials
  array[T] real<lower=0> RTpass;  // Reaction times for play trials
}

parameters {
  real<lower=-3, upper=3> boundary_pr; // Boundary separation (a)
  real<lower=-3, upper=3> tau_pr;      // Non-decision time (tau)
  real<lower=-3, upper=3> beta_pr;     // Starting point
  real<lower=-3, upper=3> drift;       // Drift rate
}

transformed parameters {
  real<lower=0> 		   boundary;
  real<lower=RTbound, upper=minRT> tau;
  real<lower=0, upper=1> 	   beta;


  boundary = exp(inv_logit(boundary_pr) * 10 - 5);
  tau 	   = inv_logit(tau_pr) * (minRT - RTbound - 1e-6) * 0.99 + RTbound; // Ensures tau will always be at least 1% less than minRT
  beta     = inv_logit(beta_pr);
}

model {
  // Priors
  boundary_pr ~ normal(0, 1);
  tau_pr      ~ normal(0, 1);
  beta_pr     ~ normal(0, 1);
  drift	      ~ normal(0, 1);
  
  // Likelihood
  RTplay[1:Nplay] ~ wiener(boundary, tau, beta, drift);
  RTpass[1:Npass] ~ wiener(boundary, tau, 1-beta, -drift);
}
