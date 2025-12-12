functions {
  #include metad.stanfunctions
}

data {
  int<lower=1> N;                             // number of data points
  int<lower=1> P;                             // number of participants
  int<lower=2> K;                             // number of confidence values
  array[N] int<lower=0, upper=1> stimulus;    // stimulus per trial  (1: present, 0: absent)
  array[N] int<lower=0, upper=1> response;    // response on each trial (1: present, 0: absent)
  array[N] int<lower=1, upper=K> confidence;  // confidence on each trial (1: low, K: high)
  array[N] int<lower=1, upper=P> participant; // participant number on each trial
  int<lower=0, upper=1> prior_only;           // sample from (1) the prior or (2) the posterior

  // prior values
  real prior_sd_mu_d_prime;          // width of normal prior on mean d_prime (type-1 sensitivity)
  real prior_sd_mu_c;                // width of normal prior on mean c (type-1 threshold)
  real prior_sd_mu_log_M;            // width of lognormal prior on mean M (metacognitive efficiency)
  real prior_mean_mu_meta_c2;        // mean of lognormal prior on mean meta_c2 (type-2 thresholds)
  real prior_sd_mu_meta_c2;          // width of lognormal prior on mean meta_c2 (type-2 thresholds)
  real prior_sd_sigma_d_prime;       // width of half-normal prior on participant-level SD of d_prime (type-1 sensitivity)
  real prior_sd_sigma_c;             // width of half-normal prior on participant-level SD of c (type-1 threshold)
  real prior_sd_sigma_log_M;         // width of half-normal prior on participant-level SD of M (metacognitive efficiency)
  real prior_sd_sigma_meta_c2;       // width of half-normal prior on participant-level SD of meta_c2 (type-2 thresholds)
  // add prior values for correlation matrices
  real prior_eta;
}

transformed data {
  int k = K-1;
  
  // sum responses into signal detection table
  // C[participant, stimulus, response*confidence]
  // where C[p,s] = [0/K, ... 0/1, 1/1, ... 1/K]
  array[P, 2, 2*K] int C = rep_array(0, P, 2, 2*K);
  for (n in 1:N) {
    if (response[n])
      C[participant[n], stimulus[n]+1, K+confidence[n]] += 1;
    else
      C[participant[n], stimulus[n]+1, K+1-confidence[n]] += 1;
  }
}

parameters {
  // Mean signal detection parameters across participants
  real mu_d_prime;
  real mu_c;
  real mu_log_M;
  vector[k] mu_z_meta_c2_0;
  vector[k] mu_z_meta_c2_1;
  
  // Standard deviation of signal detection parameters across participants
  real<lower=0> sigma_d_prime;
  real<lower=0> sigma_c;
  real<lower=0> sigma_log_M;
  vector<lower=0>[k] sigma_meta_c2;
  
  // Participant(condition)-level z-scored signal detection parameters
  vector[P] z_d_prime;
  vector[P] z_c;
  vector[P] z_log_M;
  matrix[k, P] z_meta_c2_0;
  matrix[k, P] z_meta_c2_1;
  
  // correlations between thresholds across confidence levels
  cholesky_factor_corr[k] L_Omega_meta_c2;
}

transformed parameters {
  // Participant(condition)-level SDT parameters
  vector[P] d_prime = mu_d_prime + z_d_prime * sigma_d_prime;
  vector[P] c = mu_c + z_c * sigma_c;
  vector[P] M = exp(mu_log_M + z_log_M * sigma_log_M);
  vector[P] meta_d_prime = M .* d_prime;
  vector[P] meta_c = M .* c;

  // determine type 2 thresholds using matrix-normal distribution
  matrix[k, P] meta_c2_0;
  matrix[k, P] meta_c2_1;
  {
    // pre-compute covariance between confidence thresholds
    matrix[k,k] L_Sigma_meta_c2 = diag_pre_multiply(sigma_meta_c2, L_Omega_meta_c2);
    
    matrix[k, P] meta_C2_0 = exp(rep_matrix(mu_z_meta_c2_0, P) + L_Sigma_meta_c2 * z_meta_c2_0);
    matrix[k, P] meta_C2_1 = exp(rep_matrix(mu_z_meta_c2_1, P) + L_Sigma_meta_c2 * z_meta_c2_1);
    for (p in 1:P) {
      meta_c2_0[,p] = meta_c[p] - cumulative_sum(meta_C2_0[,p]);
      meta_c2_1[,p] = meta_c[p] + cumulative_sum(meta_C2_1[,p]);
    }
  }
}

model {
  // priors for mean SDT parameters
  target += normal_lpdf(mu_d_prime | 0, prior_sd_mu_d_prime);
  target += normal_lpdf(mu_c | 0, prior_sd_mu_c);
  target += normal_lpdf(mu_log_M | 0, prior_sd_mu_log_M);
  target += normal_lpdf(mu_z_meta_c2_0 | prior_mean_mu_meta_c2, prior_sd_mu_meta_c2);
  target += normal_lpdf(mu_z_meta_c2_1 | prior_mean_mu_meta_c2, prior_sd_mu_meta_c2);

  // priors for standard deviation of SDT parameters
  target += normal_lpdf(sigma_d_prime | 0, prior_sd_sigma_d_prime);
  target += normal_lpdf(sigma_c | 0, prior_sd_sigma_c);
  target += normal_lpdf(sigma_log_M | 0, prior_sd_sigma_log_M);
  target += normal_lpdf(sigma_meta_c2 | 0, prior_sd_sigma_meta_c2);
  
  // standard normal priors for z-scored participant-level SDT parameters
  target += std_normal_lpdf(z_d_prime);
  target += std_normal_lpdf(z_c);
  target += std_normal_lpdf(z_log_M);
  target += std_normal_lpdf(to_vector(z_meta_c2_0));
  target += std_normal_lpdf(to_vector(z_meta_c2_1));
  
  // LKJ priors for correlation matrices
  target += lkj_corr_cholesky_lpdf(L_Omega_meta_c2 | prior_eta);
  
  if (!prior_only) {
    for (p in 1:P) {
      for (s in 0:1) {
        target += multinomial_lpmf(C[p, s+1] |
                                   metad_joint_pmf(s, d_prime[p], c[p],
                                                   meta_d_prime[p], meta_c[p],
                                                   meta_c2_0[,p], meta_c2_1[,p]));
      }
    }
  }
}

generated quantities {
  // Mean SDT parameters across participants
  real mu_c_prime = mu_c / mu_d_prime;
  real mu_meta_d_prime = exp(mu_log_M) * mu_d_prime;
  real mu_meta_c = exp(mu_log_M) * mu_c;
  vector[k] mu_meta_c2_0 = mu_meta_c - cumulative_sum(exp(mu_z_meta_c2_0));
  vector[k] mu_meta_c2_1 = mu_meta_c + cumulative_sum(exp(mu_z_meta_c2_1));
  
  
  array[2] simplex[2] mu_theta_1;      // Average type-1 response probabilities
  array[2, 2] simplex[K] mu_theta_2;   // Average type-2 response probabilities
  for (i in 0:1) {
    for (j in 0:1) {
      mu_theta_1[i+1, j+1] = type1_pmf(i, j, mu_d_prime, mu_c);
      mu_theta_2[i+1, j+1] = type2_pmf(i, j, mu_meta_d_prime, mu_meta_c,
                                       mu_meta_c2_0, mu_meta_c2_1);
    }
  }

  array[2*K-1, 2] real mu_ROC_1 = type1_ROC(mu_theta_1, mu_theta_2);   // Average pseudo type-1 ROC
  array[2, K-1, 2] real mu_ROC_2 = type2_ROC(mu_theta_2);              // Average type-2 ROC
  
  
  
  array[P, 2] simplex[2] theta_1;      // Participant-level type-1 response probabilities
  array[P, 2, 2] simplex[K] theta_2;   // Participant-level type-1 response probabilities
  array[P, 2*K-1, 2] real ROC_1;       // Participant-level pseudo type-1 ROC
  array[P, 2, K-1, 2] real ROC_2;      // Participant-level type-2 ROC
  for (p in 1:P) {
    for (i in 0:1) {
      for (j in 0:1) {
        theta_1[p, i+1, j+1] = type1_pmf(i, j, d_prime[p], c[p]);
        theta_2[p, i+1, j+1] = type2_pmf(i, j, meta_d_prime[p], meta_c[p],
                                         meta_c2_0[,p], meta_c2_1[,p]);
      }
    }

    ROC_1[p] = type1_ROC(theta_1[p], theta_2[p]);
    ROC_2[p] = type2_ROC(theta_2[p]);
  }
  
  // Correlation matrices
  corr_matrix[k] Omega_meta_c2 = multiply_lower_tri_self_transpose(L_Omega_meta_c2);
}
