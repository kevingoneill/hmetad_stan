functions {
  #include metad.stanfunctions
}

data {
  int<lower=1> N;                             // number of data points
  int<lower=2> K;                             // number of confidence values
  array[N] int<lower=0, upper=1> stimulus;    // stimulus per trial  (1: present, 0: absent)
  array[N] int<lower=0, upper=1> response;    // response on each trial (1: present, 0: absent)
  array[N] int<lower=1, upper=K> confidence;  // confidence on each trial (1: low, K: high)

  int<lower=1> N_x;           // number of predictors
  matrix[N, N_x] X;           // design matrix (used for d_prime, c, log_M, and z_meta_c2)

  // design matrix for predicted values
  int<lower=1> N_pred;
  matrix[N_pred, N_x] X_pred;
  
  int<lower=0, upper=1> prior_only;           // sample from (1) the prior or (2) the posterior

  // prior values
  real prior_sd_d_prime;         // width of normal prior on d_prime (type-1 sensitivity)
  real prior_sd_c;               // width of normal prior on c (type-1 threshold)
  real prior_sd_log_M;           // width of lognormal prior on M (metacognitive efficiency)
  real prior_mean_meta_c2;       // mean of lognormal prior on meta_c2 (type-2 thresholds)
  real prior_sd_meta_c2;         // width of lognormal prior on meta_c2 (type-2 thresholds)
}

transformed data {
  int k = K-1; // number of confidence thresholds
  
  array[N] int<lower=0, upper=2*K> joint_response;
  for (n in 1:N) {
    if (response[n])
      joint_response[n] = K + confidence[n];
    else
      joint_response[n] = K + 1 - confidence[n];
  }
}

parameters {
  // Type-1 signal detection parameters
  vector[N_x] beta_d_prime;
  vector[N_x] beta_c;
  
  // Type-2 signal detection paramters
  vector[N_x] beta_log_M; // metacognitive efficiency (meta_d_prime / d_prime)

  matrix[N_x,k] beta_z_meta_c2_0;  // type-2 threshold for response=0 (on transformed scale for efficiency)
  matrix[N_x,k] beta_z_meta_c2_1;  // type-2 threshold for response=1 (on transformed scale for efficiency)
}

transformed parameters {
  
}

model {
  // priors for SDT parameters
  target += normal_lpdf(beta_d_prime | 0, prior_sd_d_prime);
  target += normal_lpdf(beta_c | 0, prior_sd_c);
  target += normal_lpdf(beta_log_M | 0, prior_sd_log_M);
  target += normal_lpdf(to_vector(beta_z_meta_c2_0) | prior_mean_meta_c2, prior_sd_meta_c2);
  target += normal_lpdf(to_vector(beta_z_meta_c2_1) | prior_mean_meta_c2, prior_sd_meta_c2);
  
  if (!prior_only) {
    // predict SDT parameters per observation
    vector[N] d_prime = X * beta_d_prime;
    vector[N] c = X * beta_c;
    vector[N] M = exp(X * beta_log_M);
    vector[N] meta_d_prime = M .* d_prime;
    vector[N] meta_c = M .* c;
    matrix[k,N] d_meta_c2_0 = exp(X * beta_z_meta_c2_0)';
    matrix[k,N] d_meta_c2_1 = exp(X * beta_z_meta_c2_1)';
    
    vector[N] ll = rep_vector(0, N);
    for (n in 1:N) {
      ll[n] += categorical_lpmf(joint_response[n] |
				metad_joint_pmf(stimulus[n], d_prime[n], c[n],
						meta_d_prime[n], meta_c[n],
						meta_c[n]-cumulative_sum(d_meta_c2_0[,n]),
						meta_c[n]+cumulative_sum(d_meta_c2_1[,n])));
    }
    target += sum(ll);
  }
}

generated quantities {
  // predict SDT parameters per observation
  vector[N_pred] d_prime = X_pred * beta_d_prime;
  vector[N_pred] c = X_pred * beta_c;
  vector[N_pred] M = exp(X_pred * beta_log_M);
  vector[N_pred] meta_d_prime = M .* d_prime;
  vector[N_pred] meta_c = M .* c;
  matrix[k,N_pred] meta_c2_0 = exp(X_pred * beta_z_meta_c2_0)';
  matrix[k,N_pred] meta_c2_1 = exp(X_pred * beta_z_meta_c2_1)';
  for (n in 1:N_pred) {
    meta_c2_0[,n] = meta_c[n] - cumulative_sum(meta_c2_0[,n]);
    meta_c2_1[,n] = meta_c[n] + cumulative_sum(meta_c2_1[,n]);
  }  
  
  vector[N_pred] c_prime = c ./ d_prime;
  
  // type-1 and type-2 response distributions
  array[N_pred, 2] simplex[2] theta_1;
  array[N_pred, 2, 2] simplex[K] theta_2;
  
  // pseudo type-1 ROC and type-2 ROC
  array[N_pred, 2*K-1, 2] real ROC_1;
  array[N_pred, 2, K-1, 2] real ROC_2;
  for (n in 1:N_pred) {
    for (i in 0:1) {
      for (j in 0:1) {
	theta_1[n, i+1, j+1] = type1_pmf(i, j, d_prime[n], c[n]);
	theta_2[n, i+1, j+1] = type2_pmf(i, j, meta_d_prime[n], meta_c[n],
					 meta_c2_0[,n], meta_c2_1[,n]);
      }
    }
    
    ROC_1[n] = type1_ROC(theta_1[n], theta_2[n]);
    ROC_2[n] = type2_ROC(theta_2[n]);
  }
}
