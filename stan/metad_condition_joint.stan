functions {
  #include metad.stanfunctions
}

data {
  int<lower=1> N;                             // number of data points
  int<lower=2> K;                             // number of confidence values
  int<lower=1> B;                             // number of (between-participant) conditions
  array[N] int<lower=0, upper=1> stimulus;    // stimulus per trial  (1: present, 0: absent)
  array[N] int<lower=0, upper=1> response;    // response on each trial (1: present, 0: absent)
  array[N] int<lower=1, upper=K> confidence;  // confidence on each trial (1: low, K: high)
  array[N] int<lower=1, upper=B> condition;   // condition number on each trial
  
  int<lower=0, upper=1> prior_only;           // sample from (1) the prior or (2) the posterior

  // prior values
  real prior_sd_d_prime;         // width of normal prior on d_prime (type-1 sensitivity)
  real prior_sd_c;               // width of normal prior on c (type-1 threshold)
  real prior_sd_log_M;           // width of lognormal prior on M (metacognitive efficiency)
  real prior_mean_meta_c2;       // mean of lognormal prior on meta_c2 (type-2 thresholds)
  real prior_sd_meta_c2;         // width of lognormal prior on meta_c2 (type-2 thresholds)
}

transformed data {
  int k = K-1;

  // sum responses into signal detection table
  // C[stimulus, response*confidence]
  // where C[s] = [0/K, ... 0/1, 1/1, ... 1/K]
  array[B, 2, 2*K] int C = rep_array(0, B, 2, 2*K);
  for (n in 1:N) {
    if (response[n])
      C[condition[n], stimulus[n]+1, K+confidence[n]] += 1;
    else
      C[condition[n], stimulus[n]+1, K+1-confidence[n]] += 1;
  }
}

parameters {
  // Type-1 signal detection parameters
  vector[B] d_prime;
  vector[B] c;

  // Type-2 signal detection paramters
  vector<lower=0>[B] M;        // metacognitive efficiency (meta_d_prime / d_prime)
  matrix[k, B] z_meta_c2_0;  // type-2 threshold for response=0 (on transformed scale for efficiency)
  matrix[k, B] z_meta_c2_1;  // type-2 threshold for response=1 (on transformed scale for efficiency)
}

transformed parameters {  
  vector[B] meta_d_prime = M .* d_prime;       // estimated type-1 sensitivity based on confidence
  vector[B] meta_c = M .* c;                   // estimated type-1 bias based on confidence such that meta_c_prime = c_prime
  
  // transform type-2 thresholds to their true scale (index = confidence level)
  matrix[k, B] meta_c2_0;
  matrix[k, B] meta_c2_1;
  for (b in 1:B) {
    meta_c2_0[,b] = meta_c[b] - cumulative_sum(exp(z_meta_c2_0[,b]));
    meta_c2_1[,b] = meta_c[b] + cumulative_sum(exp(z_meta_c2_1[,b]));
  }
}

model {
  // priors for SDT parameters
  target += normal_lpdf(d_prime | 0, prior_sd_d_prime);
  target += normal_lpdf(c | 0, prior_sd_c);
  target += lognormal_lpdf(M | 0, prior_sd_log_M);
  target += normal_lpdf(to_vector(z_meta_c2_0) | prior_mean_meta_c2, prior_sd_meta_c2);
  target += normal_lpdf(to_vector(z_meta_c2_1) | prior_mean_meta_c2, prior_sd_meta_c2);
  
  if (!prior_only) {
    for (b in 1:B) {
      // model responses (hits/FAs only) with binomial distribution
      for (s in 0:1) {
        target += multinomial_lpmf(C[b, s+1] |
				   metad_joint_pmf(s, d_prime[b], c[b],
						   meta_d_prime[b], meta_c[b],
						   meta_c2_0[,b], meta_c2_1[,b]));
      }
    }
  }
}

generated quantities {
  vector[B] c_prime = c ./ d_prime;
  
  // type-1 and type-2 response distributions
  array[B, 2] simplex[2] theta_1;
  array[B, 2, 2] simplex[K] theta_2;

  // pseudo type-1 ROC and type-2 ROC
  array[B, 2*K-1, 2] real ROC_1;
  array[B, 2, K-1, 2] real ROC_2;

  for (b in 1:B) {
    for (i in 0:1) {
      for (j in 0:1) {
        theta_1[b, i+1, j+1] = type1_pmf(i, j, d_prime[b], c[b]);
        theta_2[b, i+1, j+1] = type2_pmf(i, j, meta_d_prime[b], meta_c[b],
                                         meta_c2_0[,b], meta_c2_1[,b]);
      }
    }

    ROC_1[b] = type1_ROC(theta_1[b], theta_2[b]);
    ROC_2[b] = type2_ROC(theta_2[b]);
  }
}
