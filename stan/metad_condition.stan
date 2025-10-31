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
  // tally up responses for binomial likelihood
  // R[stimulus, response]
  array[B, 2, 2] int R = rep_array(0, B, 2, 2);
  
  // tally up confidence ratings for multinomial likelihood
  // C[response, accuracy, confidence]
  array[B, 2, 2, K] int C;
  for (b in 1:B)
    C[b] = rep_array(0, 2, 2, K);

  // sum responses into signal detection table
  int k = K-1;
  for (n in 1:N) {
    int accuracy = (response[n] == stimulus[n]);
    R[condition[n], stimulus[n]+1, response[n]+1] += 1;
    C[condition[n], response[n]+1, accuracy+1, confidence[n]] += 1;
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
        target += binomial_lpmf(R[b, s+1, 2] | R[b, s+1, 1] + R[b, s+1, 2],
                                type1_pmf(s, 1, d_prime[b], c[b]));
      }
      
      // model confidence data with multinomial distribution
      for (r in 0:1) {
        for (a in 0:1)
          target += multinomial_lpmf(C[b, r+1, a+1] |
                                     type2_pmf(r, a, meta_d_prime[b],
                                               meta_c[b], meta_c2_0[,b], meta_c2_1[,b]));
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
