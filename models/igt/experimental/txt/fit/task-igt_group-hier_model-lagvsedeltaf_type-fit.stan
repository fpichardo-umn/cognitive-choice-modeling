functions {
  real igt_subject(
        array[] int choice, array[] real wins, array[] real losses,
        vector ev_exploit, vector ef, vector choice_lag, int Tsub,
        real gain, real loss,
        real update, real betaF, real phi
        ) {
    real log_lik = 0.0;
    vector[4] local_ev = ev_exploit;
    vector[4] local_ef = ef;
    vector[4] local_lag = choice_lag;
    vector[4] combined_value;
    real curUtil;
    real sign_outcome;  // Objective outcome sign
    real PEfreq;
    real efChosen;
    vector[4] PEfreq_fic;
    
    for (t in 1:Tsub) {
      local_lag += 1;
      combined_value = local_ev + betaF * local_ef + phi * local_lag;
      log_lik += categorical_logit_lpmf(choice[t] | combined_value);
      
      // Calculate subjective utility (for value learning)
      real win_component = (wins[t] == 0) ? 0.0 : pow(wins[t], gain);
      real loss_component = (losses[t] == 0) ? 0.0 : pow(losses[t], gain);
      curUtil = win_component - loss * loss_component;
      
      // Calculate objective outcome sign (for frequency learning)
      sign_outcome = wins[t] >= losses[t] ? 1.0 : -1.0;
      
      // Frequency prediction errors based on objective outcomes
      PEfreq = sign_outcome - local_ef[choice[t]];
      efChosen = local_ef[choice[t]];
      
      for (d in 1:4) {
        PEfreq_fic[d] = -sign_outcome/3.0 - local_ef[d];
      }
      
      local_ef += update * PEfreq_fic;
      local_ef[choice[t]] = efChosen + update * PEfreq;
      
      // Value learning uses subjective utility
      local_ev[choice[t]] += update * (curUtil - local_ev[choice[t]]);
      
      local_lag[choice[t]] = 0;
    }
    
    return log_lik;
  }
  
  real partial_sum_func(array[] int slice_n, int start, int end,
                        array[,] int choice, array[,] real wins, array[,] real losses,
                        array[] int Tsubj,
                        array[] real gain, array[] real loss,
                        array[] real update, array[] real betaF, array[] real phi) {
    real log_lik = 0.0;
    
    for (n in start:end) {
      vector[4] ev_exploit = rep_vector(0., 4);
      vector[4] ef = rep_vector(0., 4);
      vector[4] choice_lag = rep_vector(0., 4);
      
      log_lik += igt_subject(choice[n, 1:Tsubj[n]],
                             wins[n, 1:Tsubj[n]], losses[n, 1:Tsubj[n]],
                             ev_exploit, ef, choice_lag, Tsubj[n],
                             gain[n], loss[n],
                             update[n], betaF[n], phi[n]);
    }
    return log_lik;
  }
}

data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1> Tsubj;
  array[N, T] int<lower=1, upper=4> choice;
  array[N, T] real<lower=0> wins;
  array[N, T] real<lower=0> losses;
}

transformed data {
  array[N] int subject_indices;
  for (i in 1:N) {
    subject_indices[i] = i;
  }
}

parameters {
  array[5] real mu_pr;
  array[5] real<lower=0> sigma;
  
  array[N] real gain_pr;
  array[N] real loss_pr;
  array[N] real update_pr;
  array[N] real betaF_pr;
  array[N] real phi_pr;
}

transformed parameters {
  array[N] real<lower=0, upper=1> gain;
  array[N] real<lower=0, upper=10> loss;
  array[N] real<lower=0, upper=1> update;
  array[N] real betaF;
  array[N] real<lower=-1, upper=1> phi;
  
  gain   = to_array_1d(inv_logit(mu_pr[1] + sigma[1] * to_vector(gain_pr)));
  loss   = to_array_1d(inv_logit(mu_pr[2] + sigma[2] * to_vector(loss_pr)) * 10);
  update = to_array_1d(inv_logit(mu_pr[3] + sigma[3] * to_vector(update_pr)));
  betaF  = to_array_1d(mu_pr[4] + sigma[4] * to_vector(betaF_pr));
  phi    = to_array_1d(-1 + inv_logit(mu_pr[5] + sigma[5] * to_vector(phi_pr)) * 2);
}

model {
  mu_pr ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);
  
  gain_pr   ~ normal(0, 1);
  loss_pr   ~ normal(0, 1);
  update_pr ~ normal(0, 1);
  betaF_pr  ~ normal(0, 1);
  phi_pr    ~ normal(0, 1);
  
  int grainsize = max(1, N %/% 4);
  target += reduce_sum(partial_sum_func, subject_indices, grainsize,
                       choice, wins, losses, Tsubj,
                       gain, loss, update, betaF, phi);
}

generated quantities {
  real<lower=0, upper=1> mu_gain;
  real<lower=0, upper=10> mu_loss;
  real<lower=0, upper=1> mu_update;
  real mu_betaF;
  real<lower=-1, upper=1> mu_phi;
  
  mu_gain   = inv_logit(mu_pr[1]);
  mu_loss   = inv_logit(mu_pr[2]) * 10;
  mu_update = inv_logit(mu_pr[3]);
  mu_betaF  = mu_pr[4];
  mu_phi    = -1 + inv_logit(mu_pr[5]) * 2;
}
