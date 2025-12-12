functions {
  #include metad.stanfunctions
}

data {
  int<lower=1> N;                             // number of data points
  int<lower=2> K;                             // number of confidence values
  array[N] int<lower=0, upper=1> stimulus;    // stimulus per trial  (1: present, 0: absent)
  array[N] int<lower=0, upper=1> response;    // response on each trial (1: present, 0: absent)
  array[N] int<lower=1, upper=K> confidence;  // confidence on each trial (1: low, K: high)
  
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
  array[2, 2*K] int C = rep_array(0, 2, 2*K);
  for (n in 1:N) {
    if (response[n])
      C[stimulus[n]+1, K+confidence[n]] += 1;
    else
      C[stimulus[n]+1, K+1-confidence[n]] += 1;
  }
}

parameters {
  // Type-1 signal detection parameters
  real d_prime;
  real c;

  // Type-2 signal detection paramters
  real<lower=0> M;        // metacognitive efficiency (meta_d_prime / d_prime)
  vector[k] z_meta_c2_0;  // type-2 threshold for response=0 (on transformed scale for efficiency)
  vector[k] z_meta_c2_1;  // type-2 threshold for response=1 (on transformed scale for efficiency)
}

transformed parameters {  
  real meta_d_prime = M * d_prime;       // estimated type-1 sensitivity based on confidence
  real meta_c = M * c;                   // estimated type-1 bias based on confidence such that meta_c_prime = c_prime
  
  // transform type-2 thresholds to their true scale (index = confidence level)
  vector[k] meta_c2_0 = meta_c - cumulative_sum(exp(z_meta_c2_0));
  vector[k] meta_c2_1 = meta_c + cumulative_sum(exp(z_meta_c2_1));
}

model {
  // priors for SDT parameters
  target += normal_lpdf(d_prime | 0, prior_sd_d_prime);
  target += normal_lpdf(c | 0, prior_sd_c);
  target += lognormal_lpdf(M | 0, prior_sd_log_M);
  target += normal_lpdf(z_meta_c2_0 | prior_mean_meta_c2, prior_sd_meta_c2);
  target += normal_lpdf(z_meta_c2_1 | prior_mean_meta_c2, prior_sd_meta_c2);
  
  if (!prior_only) {
    // model responses (hits/FAs only) with binomial distribution
    for (s in 0:1) {
      target += multinomial_lpmf(C[s+1] |
                                 metad_joint_pmf(s, d_prime, c,
                                                 meta_d_prime, meta_c,
                                                 meta_c2_0, meta_c2_1));
    }
  }
}

generated quantities {
  real c_prime = c / d_prime;
  
  // type-1 and type-2 response distributions
  array[2] simplex[2] theta_1;
  array[2, 2] simplex[K] theta_2;
  for (i in 0:1) {
    for (j in 0:1) {
      theta_1[i+1, j+1] = type1_pmf(i, j, d_prime, c);
      theta_2[i+1, j+1] = type2_pmf(i, j, meta_d_prime, meta_c, meta_c2_0, meta_c2_1);
    }
  }

  // pseudo type-1 ROC and type-2 ROC
  array[2*K-1, 2] real ROC_1 = type1_ROC(theta_1, theta_2);
  array[2, K-1, 2] real ROC_2 = type2_ROC(theta_2);
}
