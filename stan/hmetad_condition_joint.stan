functions {
  #include metad.stanfunctions
}

data {
  int<lower=1> N;                             // number of data points
  int<lower=1> P;                             // number of participants
  int<lower=2> K;                             // number of confidence values
  int<lower=1> W;                             // number of conditions
  array[N] int<lower=0, upper=1> stimulus;    // stimulus per trial  (1: present, 0: absent)
  array[N] int<lower=0, upper=1> response;    // response on each trial (1: present, 0: absent)
  array[N] int<lower=1, upper=K> confidence;  // confidence on each trial (1: low, K: high)
  array[N] int<lower=1, upper=P> participant; // participant number on each trial
  array[N] int<lower=1, upper=W> condition;     // condition number on each trial
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
  // C[stimulus, response*confidence]
  // where C[s] = [0/K, ... 0/1, 1/1, ... 1/K]
  array[P, W, 2, 2*K] int C;
  for (p in 1:P) {
    for (w in 1:W) {
      C[p, w] = rep_array(0, 2, 2*K);
    }
  }
  
  for (n in 1:N) {
    if (response[n])
      C[participant[n], condition[n], stimulus[n]+1, K+confidence[n]] += 1;
    else
      C[participant[n], condition[n], stimulus[n]+1, K+1-confidence[n]] += 1;
  }
}

parameters {
  // Mean signal detection parameters across participants
  vector[W] mu_d_prime;
  vector[W] mu_c;
  vector[W] mu_log_M;
  matrix[k,W] mu_z_meta_c2_0;
  matrix[k,W] mu_z_meta_c2_1;
  
  // Standard deviation of signal detection parameters across participants
  vector<lower=0>[W] sigma_d_prime;
  vector<lower=0>[W] sigma_c;
  vector<lower=0>[W] sigma_log_M;
  vector<lower=0>[W] sigma_meta_c2_condition;
  vector<lower=0>[k] sigma_meta_c2_confidence;
  
  // Participant(condition)-level z-scored signal detection parameters
  matrix[W, P] z_d_prime;
  matrix[W, P] z_c;
  matrix[W, P] z_log_M;
  array[P] matrix[k,W] z_meta_c2_0;
  array[P] matrix[k,W] z_meta_c2_1;
  
  // Correlation matrices
  cholesky_factor_corr[W] L_Omega_d_prime;
  cholesky_factor_corr[W] L_Omega_c;
  cholesky_factor_corr[W] L_Omega_log_M;
  cholesky_factor_corr[W] L_Omega_meta_c2_condition;   // correlations between thresholds across conditions
  cholesky_factor_corr[k] L_Omega_meta_c2_confidence;  // correlations between thresholds across confidence levels
}

transformed parameters {
  // Participant(condition)-level SDT parameters
  matrix[W, P] d_prime = rep_matrix(mu_d_prime, P) +
       diag_pre_multiply(sigma_d_prime, L_Omega_d_prime) * z_d_prime;
  matrix[W, P] c = rep_matrix(mu_c, P) + diag_pre_multiply(sigma_c, L_Omega_c) * z_c;
  matrix[W, P] M = exp(rep_matrix(mu_log_M, P) +
                       diag_pre_multiply(sigma_log_M, L_Omega_log_M) * z_log_M);
  matrix[W, P] meta_d_prime = M .* d_prime;
  matrix[W, P] meta_c = M .* c;

  // determine type 2 thresholds using matrix-normal distribution
  array[W, P] vector[k] meta_c2_0;
  array[W, P] vector[k] meta_c2_1;
  {
    // get covariance matrices across conditions & confidence levels
    matrix[W,W] L_Sigma_meta_c2_condition =
      diag_pre_multiply(sigma_meta_c2_condition,
                        L_Omega_meta_c2_condition);
    matrix[k,k] L_Sigma_meta_c2_confidence =
      diag_pre_multiply(sigma_meta_c2_confidence,
                        L_Omega_meta_c2_confidence);
    for (p in 1:P) {
      // pre-compute threshold differences for efficiency
      matrix[k, W] c2_0 = exp(mu_z_meta_c2_0 +
                              L_Sigma_meta_c2_confidence *
                              z_meta_c2_0[p] *
                              L_Sigma_meta_c2_condition);
      matrix[k, W] c2_1 = exp(mu_z_meta_c2_1 +
                              L_Sigma_meta_c2_confidence *
                              z_meta_c2_1[p] *
                              L_Sigma_meta_c2_condition);
      for (w in 1:W) {
        meta_c2_0[w,p] = meta_c[w,p] - cumulative_sum(c2_0[,w]);
        meta_c2_1[w,p] = meta_c[w,p] + cumulative_sum(c2_1[,w]);
      }
    }
  }
}

model {
  // priors for mean SDT parameters
  target += normal_lpdf(mu_d_prime | 0, prior_sd_mu_d_prime);
  target += normal_lpdf(mu_c | 0, prior_sd_mu_c);
  target += normal_lpdf(mu_log_M | 0, prior_sd_mu_log_M);
  target += normal_lpdf(to_vector(mu_z_meta_c2_0) | prior_mean_mu_meta_c2, prior_sd_mu_meta_c2);
  target += normal_lpdf(to_vector(mu_z_meta_c2_1) | prior_mean_mu_meta_c2, prior_sd_mu_meta_c2);

  // priors for standard deviation of SDT parameters
  target += normal_lpdf(sigma_d_prime | 0, prior_sd_sigma_d_prime);
  target += normal_lpdf(sigma_c | 0, prior_sd_sigma_c);
  target += normal_lpdf(sigma_log_M | 0, prior_sd_sigma_log_M);
  target += normal_lpdf(sigma_meta_c2_condition | 0, prior_sd_sigma_meta_c2);
  target += normal_lpdf(sigma_meta_c2_confidence | 0, prior_sd_sigma_meta_c2);
  
  // standard normal priors for z-scored participant-level SDT parameters
  target += std_normal_lpdf(to_vector(z_d_prime));
  target += std_normal_lpdf(to_vector(z_c));
  target += std_normal_lpdf(to_vector(z_log_M));
  for (p in 1:P) {
    target += std_normal_lpdf(to_vector(z_meta_c2_0[p]));
    target += std_normal_lpdf(to_vector(z_meta_c2_1[p]));
  }
  
  // LKJ priors for correlation matrices
  target += lkj_corr_cholesky_lpdf(L_Omega_d_prime | prior_eta);
  target += lkj_corr_cholesky_lpdf(L_Omega_c | prior_eta);
  target += lkj_corr_cholesky_lpdf(L_Omega_log_M | prior_eta);
  target += lkj_corr_cholesky_lpdf(L_Omega_meta_c2_condition | prior_eta);
  target += lkj_corr_cholesky_lpdf(L_Omega_meta_c2_confidence | prior_eta);
  
  if (!prior_only) {
    for (p in 1:P) {
      for (w in 1:W) { // loop through conditions
        // model responses & confidence with binomial distribution 
        for (s in 0:1) {
          target += multinomial_lpmf(C[p, w, s+1] |
				     metad_joint_pmf(s, d_prime[w,p], c[w,p],
						     meta_d_prime[w,p], meta_c[w,p],
						     meta_c2_0[w,p], meta_c2_1[w,p]));
        }
      }
    }
  }
}

generated quantities {
  // Mean SDT parameters across participants
  vector[W] mu_c_prime = mu_c ./ mu_d_prime;
  vector[W] mu_meta_d_prime = exp(mu_log_M) .* mu_d_prime;
  vector[W] mu_meta_c = exp(mu_log_M) .* mu_c;
  array[W] vector[k] mu_meta_c2_0;
  array[W] vector[k] mu_meta_c2_1;
  for (w in 1:W) {
    mu_meta_c2_0[w] = mu_meta_c[w] - cumulative_sum(exp(mu_z_meta_c2_0[,w]));
    mu_meta_c2_1[w] = mu_meta_c[w] + cumulative_sum(exp(mu_z_meta_c2_1[,w]));
  }

  
  
  array[W, 2] simplex[2] mu_theta_1;      // Average type-1 response probabilities
  array[W, 2, 2] simplex[K] mu_theta_2;   // Average type-2 response probabilities
  array[W, 2*K-1, 2] real mu_ROC_1;       // Average pseudo type-1 ROC
  array[W, 2, K-1, 2] real mu_ROC_2;      // Average type-2 ROC
  for (w in 1:W) {
    for (i in 0:1) {
      for (j in 0:1) { 
        mu_theta_1[w, i+1, j+1] = type1_pmf(i, j, mu_d_prime[w], mu_c[w]);
        mu_theta_2[w, i+1, j+1] = type2_pmf(i, j, mu_meta_d_prime[w], mu_meta_c[w],
                                                  mu_meta_c2_0[w], mu_meta_c2_1[w]);
      }
    }

    mu_ROC_1[w] = type1_ROC(mu_theta_1[w], mu_theta_2[w]);
    mu_ROC_2[w] = type2_ROC(mu_theta_2[w]);
  }
  
  
  array[W, P, 2] simplex[2] theta_1;      // Participant-level type-1 response probabilities
  array[W, P, 2, 2] simplex[K] theta_2;   // Participant-level type-1 response probabilities
  array[W, P, 2*K-1, 2] real ROC_1;       // Participant-level pseudo type-1 ROC
  array[W, P, 2, K-1, 2] real ROC_2;      // Participant-level type-2 ROC
  for (p in 1:P) {
    for (w in 1:W) {
      for (i in 0:1) {
        for (j in 0:1) {
          theta_1[w, p, i+1, j+1] = type1_pmf(i, j, d_prime[w, p], c[w, p]);
          theta_2[w, p, i+1, j+1] = type2_pmf(i, j, meta_d_prime[w, p], meta_c[w, p],
                                              meta_c2_0[w,p], meta_c2_1[w,p]);
        }
      }

      ROC_1[w, p] = type1_ROC(theta_1[w, p], theta_2[w, p]);
      ROC_2[w, p] = type2_ROC(theta_2[w, p]);
    }
  }

  
  
  // Correlation matrices
  corr_matrix[W] Omega_d_prime = multiply_lower_tri_self_transpose(L_Omega_d_prime);
  corr_matrix[W] Omega_c = multiply_lower_tri_self_transpose(L_Omega_c);
  corr_matrix[W] Omega_log_M = multiply_lower_tri_self_transpose(L_Omega_log_M);
  corr_matrix[W] Omega_meta_c2_condition =
    multiply_lower_tri_self_transpose(L_Omega_meta_c2_condition);
  corr_matrix[k] Omega_meta_c2_confidence =
    multiply_lower_tri_self_transpose(L_Omega_meta_c2_confidence);
}
