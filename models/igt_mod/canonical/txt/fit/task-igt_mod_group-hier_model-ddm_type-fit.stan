// Optimized Hierarchical DDM Model for IGT_MOD
data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1> sid;
  int<lower=1> Nplay_max;
  int<lower=1> Npass_max;
  real<lower=0> RTbound;
  array[N] int<lower=0> Nplay;
  array[N] int<lower=0> Npass;
  array[N, T] real RTplay;
  array[N, T] real RTpass;
  vector[N] minRT;
}

transformed data {
  vector[N] minRTdiff = minRT - RTbound - 1e-6;
  real RTmax = max(minRT);
} 

parameters {
  array[4] real mu_pr;
  array[4] real<lower=0> sigma;

  array[N] real boundary_pr;
  array[N] real tau_pr;
  array[N] real beta_pr;
  array[N] real drift_pr;
}

transformed parameters {
  array[N] real<lower=1e-6> boundary;
  array[N] real<lower=RTbound - 1e-5, upper=RTmax> tau;
  array[N] real<lower=0, upper=1> beta;
  array[N] real drift;

  boundary = to_array_1d(exp(inv_logit(mu_pr[1] + sigma[1] .* to_vector(boundary_pr)) * 10 - 5));
  tau      = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(tau_pr)) .* (minRTdiff * 0.99) + RTbound);
  beta     = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(beta_pr)));
  drift    = to_array_1d(mu_pr[4] + sigma[4] .* to_vector(drift_pr));
}

model {
  mu_pr ~ normal(0, 5);
  sigma ~ student_t(3, 0, 1);

  boundary_pr ~ normal(0, 5);
  tau_pr      ~ normal(0, 2.5);
  beta_pr     ~ normal(0, 2.5);
  drift_pr    ~ normal(0, 5);
  
  for (n in 1:N) { 
    RTplay[n, 1:Nplay[n]] ~ wiener(boundary[n], tau[n],   beta[n],  drift[n]);
    RTpass[n, 1:Npass[n]] ~ wiener(boundary[n], tau[n], 1-beta[n], -drift[n]);
  }
}

generated quantities {
  real<lower=RTbound, upper=RTmax> mu_tau;
  real<lower=0, upper=1> mu_beta;
  real<lower=1e-6> mu_boundary;
  real mu_drift = mu_pr[4];

  {
    array[N] real mu_transformed = inv_logit(mu_pr);
    real RTlowerbound = (mean(minRT) - RTbound) + RTbound;

    mu_boundary = exp(mu_transformed[1] * 10 - 5);
    mu_tau  = mu_transformed[2] * RTlowerbound;
    mu_beta = mu_transformed[3];
  }
}
