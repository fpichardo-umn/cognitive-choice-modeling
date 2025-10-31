// Optimized Hierarchical DDM b1p2 Model for IGT_MOD
data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1> sid;
  array[N] int<lower=1> Tsubj;
  array[N] real<lower=0> minRT;
  real RTbound;
  array[N, T] real<lower=0> RT;
  array[N, T] int<lower=0, upper=1> choice;
}

transformed data {
  int block = 20; 
}

parameters {
  array[6] real mu_pr;
  array[6] real<lower=0> sigma;

  array[N] real boundary1_pr;
  array[N] real boundary_pr;
  array[N] real tau1_pr;
  array[N] real tau_pr;
  array[N] real beta_pr;
  array[N] real drift;
}

transformed parameters {
  array[N] real<lower=0, upper=6> boundary1;
  array[N] real<lower=0, upper=6> boundary;
  array[N] real<lower=RTbound, upper=max(minRT)> tau1;
  array[N] real<lower=RTbound, upper=max(minRT)> tau;
  array[N] real<lower=0, upper=1> beta;

  boundary1 = to_array_1d(inv_logit(mu_pr[1] + sigma[1] .* to_vector(boundary1_pr)) * 5 + 0.01);
  boundary  = to_array_1d(inv_logit(mu_pr[2] + sigma[2] .* to_vector(boundary_pr)) * 5 + 0.01);
  tau1      = to_array_1d(inv_logit(mu_pr[3] + sigma[3] .* to_vector(tau1_pr)) .* (to_vector(minRT) - RTbound - 1e-6) * 0.99 + RTbound);
  tau       = to_array_1d(inv_logit(mu_pr[4] + sigma[4] .* to_vector(tau_pr)) .* (to_vector(minRT) - RTbound - 1e-6) * 0.99 + RTbound);
  beta      = to_array_1d(inv_logit(mu_pr[5] + sigma[5] .* to_vector(beta_pr)));
}

model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);

  boundary1_pr ~ normal(0, 1);
  boundary_pr  ~ normal(0, 1);
  tau1_pr      ~ normal(0, 1);
  tau_pr       ~ normal(0, 1);
  beta_pr      ~ normal(0, 1);
  drift        ~ normal(mu_pr[6], sigma[6]);

  // Build per-subject boundary/tau vectors
  array[N] vector[T] boundary_subj;
  array[N] vector[T] tau_subj;
  
  for (n in 1:N) {
    int Tsubj_n = Tsubj[n];
    
    // First block
    boundary_subj[n][1:block] = rep_vector(boundary1[n], block);
    tau_subj[n][1:block]      = rep_vector(tau1[n], block);
    
    // Rest of blocks
    int rest_len = Tsubj_n - block;
    boundary_subj[n][(block+1): Tsubj_n] = rep_vector(boundary[n], rest_len);
    tau_subj[n][(block+1): Tsubj_n]      = rep_vector(tau[n], rest_len);
  }
  
  for (n in 1:N) {
    array[Tsubj[n]] int play_indices;
    array[Tsubj[n]] int pass_indices;
    int play_count = 0;
    int pass_count = 0;

    for (t in 1:Tsubj[n]) {
      if (RT[n, t] != 999) {
        if (choice[n, t] == 1) {
          play_count += 1;
          play_indices[play_count] = t;
        } else {
          pass_count += 1;
          pass_indices[pass_count] = t;
        }
      }
    }
    
    RT[n, play_indices[:play_count]] ~ wiener(boundary_subj[n][play_indices[:play_count]], tau_subj[n][play_indices[:play_count]], beta[n], drift[n]);
    RT[n, pass_indices[:pass_count]] ~ wiener(boundary_subj[n][pass_indices[:pass_count]], tau_subj[n][pass_indices[:pass_count]], 1-beta[n], -drift[n]);
  }
}

generated quantities {
  real<lower=0, upper=6> mu_boundary1;
  real<lower=0, upper=6> mu_boundary;
  real<lower=0> mu_tau1;
  real<lower=0> mu_tau;
  real<lower=0, upper=1> mu_beta;
  real mu_drift;
  
  mu_boundary1 = inv_logit(mu_pr[1]) * 5 + 0.01;
  mu_boundary  = inv_logit(mu_pr[2]) * 5 + 0.01;
  mu_tau1      = inv_logit(mu_pr[3]) * (mean(minRT) - RTbound) * 0.99 + RTbound;
  mu_tau       = inv_logit(mu_pr[4]) * (mean(minRT) - RTbound) * 0.99 + RTbound;
  mu_beta      = inv_logit(mu_pr[5]);
  mu_drift     = mu_pr[6];
}
