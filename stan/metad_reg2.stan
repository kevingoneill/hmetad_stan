functions {
  // Convert a binary number x from {0, 1} to {-1, 1}
  int to_polar(int x) {
    return 2*x - 1;
  }
  
  // Calculate the probability of making a response given a stimulus
  real type1_pmf(int stimulus, int response, real d_prime, real c) {
    return Phi(to_polar(response) * (to_polar(stimulus)*d_prime/2 - c));
  }
  
  // Create a K-simplex of probabilities that confidence=k
  // (i.e., a K-vector [p(c=1), p(c=2), ... p(c=K)]) given:
  // int K: number of possible confidence ratings
  // int R: response
  // int A: accuracy of response
  // real meta_d_prime: type-1 sensitivity
  // real c: type-1 threshold
  // vector meta_c2_r: type-2 thresholds for response=r
  vector type2_pmf(int response, int accuracy, real meta_d_prime,
                   real meta_c, vector meta_c2_0, vector meta_c2_1) {
    int K = size(meta_c2_0) + 1;
    real mu = to_polar(accuracy) * meta_d_prime / 2;

    // set response-specific thresholds for vectorization
    vector[K+1] C2;
    if (response) {
      C2[1] = meta_c;
      C2[2:K] = meta_c2_1;
      C2[K+1] = positive_infinity();
    } else {
      C2[1] = -meta_c;
      C2[2:K] = -meta_c2_0;
      C2[K+1] = positive_infinity();
    }

    // compute response probabilities
    vector[K+1] theta = Phi(mu - C2);
    return (theta[1:K] - theta[2:(K+1)]) / theta[1];
  }

  // Compute the pseudo type-1 ROC given:
  //  theta_1: the type-1 response probabilities
  //  theta_2: the conditional type_2 response probabilities
  array[,] real type1_ROC(array[] vector theta_1, array[,] vector theta_2) {
    int K = size(theta_2[1,1]);
    array[2*K-1, 2] real ROC_1;
    ROC_1[K, 1] = theta_1[1, 2];
    ROC_1[K, 2] = theta_1[2, 2];
    
    for (i in 2:K) {
      // pseudo type-1 FA  rate for confidence i
      // stimulus == 0 (incorrect)
      ROC_1[K-i+1, 1] = 1 - theta_1[1,1]*sum(theta_2[1, 2, i:K]);
      ROC_1[K+i-1, 1] = theta_1[1,2]*sum(theta_2[2, 1, i:K]);
      
      // pseudo type-1 hit rate for confidence i
      // stimulus == 1 (correct)
      ROC_1[K-i+1, 2] = 1 - theta_1[2,1]*sum(theta_2[1, 1, i:K]);
      ROC_1[K+i-1, 2] = theta_1[2,2]*sum(theta_2[2, 2, i:K]);
    }

    return ROC_1;
  }

  // Compute the type-2 ROC given:
  //   theta_2: the conditional type-2 response probabilities
  array[,,] real type2_ROC(array[,] vector theta_2) {
    int K = size(theta_2[1,1]);
    array[2, K-1, 2] real ROC_2;
    
    for (i in 1:2) {
      for (j in 2:K) {
        ROC_2[i, j-1, 1] = sum(theta_2[i, 1, j:K]);  // type-2 FA rate
        ROC_2[i, j-1, 2] = sum(theta_2[i, 2, j:K]);  // type-2 hit rate	
      }
    }
    return ROC_2;
  }
}

data {
  int<lower=1> N;                             // number of data points
  int<lower=2> K;                             // number of confidence values
  array[N] int<lower=0, upper=1> stimulus;    // stimulus per trial  (1: present, 0: absent)
  array[N] int<lower=0, upper=1> response;    // response on each trial (1: present, 0: absent)
  array[N] int<lower=1, upper=K> confidence;  // confidence on each trial (1: low, K: high)

  int<lower=1> N_x;                           // number of predictors
  
  // design matrix for predicted values
  int<lower=1> N_pred;
  matrix[N_pred, N_x] X_pred;

  // indicator for which row of X_pred each observation belongs to
  array[N] int<lower=1, upper=N_pred> x;
  
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

  // tally up responses for binomial likelihood
  // R[stimulus, response]
  array[N_pred, 2, 2] int R = rep_array(0, N_pred, 2, 2);
  
  // tally up confidence ratings for multinomial likelihood
  // C[response, accuracy, confidence]
  array[N_pred, 2, 2, K] int C;
  for (n in 1:N_pred)
    C[n] = rep_array(0, 2, 2, K);

  // sum responses into signal detection table
  for (n in 1:N) {
    int accuracy = (response[n] == stimulus[n]);
    R[x[n], stimulus[n]+1, response[n]+1] += 1;
    C[x[n], response[n]+1, accuracy+1, confidence[n]] += 1;
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
}

model {
  // priors for SDT parameters
  target += normal_lpdf(beta_d_prime | 0, prior_sd_d_prime);
  target += normal_lpdf(beta_c | 0, prior_sd_c);
  target += normal_lpdf(beta_log_M | 0, prior_sd_log_M);
  target += normal_lpdf(to_vector(beta_z_meta_c2_0) | prior_mean_meta_c2, prior_sd_meta_c2);
  target += normal_lpdf(to_vector(beta_z_meta_c2_1) | prior_mean_meta_c2, prior_sd_meta_c2);
  
  if (!prior_only) {
    vector[N_pred] ll = rep_vector(0, N_pred);
    for (n in 1:N_pred) {
      // model responses (hits/FAs only) with binomial distribution
      for (s in 0:1) {
        ll += binomial_lpmf(R[n, s+1, 2] | R[n, s+1, 1] + R[n, s+1, 2],
			    type1_pmf(s, 1, d_prime[n], c[n]));
      }
      
      // model confidence data with multinomial distribution
      for (r in 0:1) {
        for (a in 0:1)
          ll += multinomial_lpmf(C[n, r+1, a+1] |
				 type2_pmf(r, a, meta_d_prime[n],
					   meta_c[n], meta_c2_0[,n], meta_c2_1[,n]));
      }
    }
    target += sum(ll);
  }
}

generated quantities {
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
