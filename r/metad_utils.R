#' to_polar(x):
#'   convert binary x from {0, 1} to {-1, 1}
to_polar <- function(x) 2*x-1

#' type1_pmf(stimulus, response, d_prime, c):
#'   Determine the probability of `response` given `stimulus` for
#'   a signal detection agent with sensitivity `d_prime` and bias `c`.
sdt_type1_pmf <- function(stimulus, response, d_prime, c) {
  pnorm(to_polar(response) * (to_polar(stimulus)*d_prime/2 - c))
}

#' type2_pmf(response, accuracy, meta_d_prime, meta_c, meta_c2_0, meta_c2_1):
#'   Determine a simplex over confidence ratings given a `response`
#'   and its `accuracy` for a signal detection agent with
#'   metacognitive sensitivity `meta_d_prime`, metacognitive bias
#'   `meta_c`, and response-specific thresholds `meta_c2_0` (for
#'   `response == 0`) and `meta_c2_1` (for `response == 1`).
#' 
sdt_type2_pmf <- function(response, accuracy, meta_d_prime, meta_c, meta_c2_0, meta_c2_1) {
  ## mean of signal distribution
  mu <- to_polar(accuracy) * meta_d_prime / 2

  ## thresholds
  K <- length(meta_c2_0) + 1
  C2 <- c(-meta_c, -meta_c2_0, Inf)
  if (response)
    C2 <- c(meta_c, meta_c2_1, Inf)

  (pnorm(mu - C2[1:K]) - pnorm(mu - C2[2:(K+1)])) / pnorm(mu - C2[1])
}

#' joint_response(response, confidence, K):
#'   convert a type 1 response (0 or 1) and a corresponding confidence level
#'   (between 1 and K) into a joint response between 1 and 2*K, with outer
#'   values reflecting confident responses and intermediate values reflecting
#'   uncertainty.
joint_response <- function(response, confidence, K) {
  ifelse(response, confidence+K, K+1-confidence)
}

response_probabilities <- function(x) {
  L <- length(x) / 2
  
  c(x[1:L] / sum(x[1:L]),
    x[(L+1):(2*L)] / sum(x[(L+1):(2*L)]))
}

sdt_joint_pmf <- function(stimulus, d_prime, c,
                          meta_d_prime, meta_c,
                          meta_c2_0, meta_c2_1, log=FALSE) {
  # number of confidence levels
  K <- length(meta_c2_0) + 1
  
  # type-1 response probabilities
  lp_1 <- pnorm(to_polar(stimulus)*d_prime/2 - c, log.p=TRUE)
  lp_0 <- pnorm(to_polar(stimulus)*d_prime/2 - c, lower.tail=FALSE, log.p=TRUE)
  
  # calculate normal cdfs (log scale)
  lp2_1 <- pnorm(to_polar(stimulus) * meta_d_prime/2 - c(meta_c, meta_c2_1), log.p=TRUE)
  lp2_0 <- pnorm(-to_polar(stimulus) * meta_d_prime/2 + c(meta_c, meta_c2_0), log.p=TRUE)
  
  # response probabilities
  log_theta <- rep(0, 2*K)
  for (k in 1:(K-1)) {
    log_theta[K-k+1] <- log(exp(lp2_0[k]) - exp(lp2_0[k+1]))
    log_theta[K+k] <- log(exp(lp2_1[k]) - exp(lp2_1[k+1]))
  }
  log_theta[1] <- lp2_0[K]
  log_theta[2*K] <- lp2_1[K]
  
  # weight by P(response|stimulus) and normalize
  log_theta[1:K] <- log_theta[1:K] + lp_0 - lp2_0[1]
  log_theta[(K+1):(2*K)] <- log_theta[(K+1):(2*K)] + lp_1 - lp2_1[1]

  if (log)
    log_theta
  else
    exp(log_theta)
}



#' sdt_type1_ROC(d_prime=1, c=seq(-10, 10, by=.05)):
#'   Get a tibble of P(Hit) and P(FA) for varying levels of d' and c
sdt_type1_ROC <- function(d_prime=1, c=0,
                          log_M=0, meta_c2=seq(-10, 10, by=.01)) {
  expand_grid(d_prime=d_prime, c=c, log_M=log_M, meta_c2=meta_c2) |>
    mutate(meta_d_prime=exp(log_M)*d_prime,
           meta_c=exp(log_M)*c,
           p_hit=pmap_dbl(list(d_prime, c, meta_d_prime, meta_c, meta_c2),
                          function(d, c, md, mc, mc2)
                            ifelse(mc2 > mc,
                                   sdt_type1_pmf(stimulus=1, response=1, d, c) *
                                     last(sdt_type2_pmf(response=1, accuracy=1,
                                                        md, mc, mc2, mc2)),
                                   1 - sdt_type1_pmf(stimulus=1, response=0, d, c) *
                                     last(sdt_type2_pmf(response=0, accuracy=0,
                                                        md, mc, mc2, mc2)))),
           p_fa=pmap_dbl(list(d_prime, c, meta_d_prime, meta_c, meta_c2),
                         function(d, c, md, mc, mc2)
                           ifelse(mc2 > mc,
                                  sdt_type1_pmf(stimulus=0, response=1, d, c) *
                                    last(sdt_type2_pmf(response=1, accuracy=0,
                                                       md, mc, mc2, mc2)),
                                  1-sdt_type1_pmf(stimulus=0, response=0, d, c) *
                                    last(sdt_type2_pmf(response=0, accuracy=1,
                                                       md, mc, mc2, mc2)))))
}


#' sim_metad(N_trials, d_prime, c, log_M, c2_0_diff, c2_1_diff, metac_absolute, summarize):
#'   Simulate from the meta-d' model with sensitivity d_prime,
#'   response bias c, metacognitive efficiency log_M, and 
#'   confidence threshold distances c2_0_diff and c2_1_diff (for the two responses).
#'
#'   Confidence thresholds are defined relative to meta-c, such that
#'   meta_c2_0 = meta_c - cumsum(c2_0_diff) and meta_c2_1 = meta_c + cumsum(c2_1_diff).
#'
#'   If metac_absolute=TRUE, meta_c = c. Otherwise, meta_c = M * c.
#' 
#'   If summarize=FALSE, returns a dataset with one row per observation.
#'   If summarize=TRUE, returns an aggregated dataset where `n` is the
#'   number of observations per response, accuracy, and confidence level.
sim_metad <- function(N_trials=100, d_prime=1, c=0, log_M=0,
                      c2_0_diff=rep(.5, 3), c2_1_diff=rep(.5, 3),
                      metac_absolute=TRUE, summarize=FALSE) {
  if (N_trials <= 0)
    stop("Error: `N_trials` must be greater than 0.")
  if (!all(length(d_prime)==1, length(c)==1, length(log_M)==1,
           is.numeric(d_prime), is.numeric(c), is.numeric(log_M)))
    stop("Error: `d_prime`, `c`, and `log_M` must be single numbers.")
  if (!is.numeric(c2_0_diff) || !is.numeric(c2_1_diff) ||
        length(c2_0_diff) != length(c2_1_diff) ||
        !all(c2_0_diff > 0) || !all(c2_1_diff > 0))
    stop("Error: c2_0_diff and c2_1_diff must be positive vectors of the same length")
  
  M <- exp(log_M)
  meta_d_prime <- M * d_prime
  meta_c <- NULL
  if (metac_absolute) {
    meta_c <- c
  } else {
    meta_c <- M * c
  }
  meta_c2_0 <- meta_c - cumsum(c2_0_diff)
  meta_c2_1 <- meta_c + cumsum(c2_1_diff)
  
  d <- expand_grid(stimulus=0:1, response=0:1, confidence=1:4) |>
    mutate(correct=as.integer(stimulus==response),
           joint_response=joint_response(response, confidence, K)) |>
    arrange(stimulus, joint_response) |>
    group_by(stimulus) |>
    mutate(theta=sdt_joint_pmf(first(stimulus), d_prime, c, meta_d_prime,
                               meta_c, meta_c2_0, meta_c2_1),
           theta_1=sdt_type1_pmf(first(stimulus), response=response, d_prime, c),
           theta_2=theta/theta_1,
           n=as.vector(rmultinom(1, N_trials/2, theta))) |>
    mutate(d_prime=d_prime, c=c, meta_d_prime=meta_d_prime,
           M=M, meta_c2_0=list(meta_c2_0), meta_c2_1=list(meta_c2_1)) |>
    select(stimulus, response, correct, confidence, n, d_prime, c, meta_d_prime,
           M, meta_c2_0, meta_c2_1, theta, theta_1, theta_2) |>
    arrange(stimulus, response, confidence)
  
  if (summarize) {
    d
  } else {
    d |>
      uncount(n) |>
      mutate(trial=row_number()) |>
      relocate(trial)
  }
}






#' cov_matrix(S, OMEGA):
#'   Generate a covariance matrix from:
#'      S: a vector of standard deviations
#'  OMEGA: a correlation matrix
cov_matrix <- function(S, OMEGA) {
  diag(S) %*% OMEGA %*% diag(S)
}

#' corr_matrix(r, nrow=2):
#'   Generate an [nrow x nrow] correlation matrix
#'   with all off-diagonal values equal to r
corr_matrix <- function(r, nrow=2) {
  diag(1-r, nrow) + r
}

#' rmatrixnorm(mu, L_sigma_1, L_sigma_2):
#'   Sample from a matrix-normal distribution with mean mu (a matrix),
#'   row-wise covariances sigma_1, and column-wise covariances sigma_2,
#'   where L_sigma_1 and L_sigma_2 are the Cholesky decomposed
#'   covariance matrices
rmatrixnorm <- function(mu, L_sigma_1, L_sigma_2) {
  mu +
    L_sigma_1 %*%
    matrix(rmvnorm(1, mean=rep(0, nrow(L_sigma_1)*nrow(L_sigma_2))),
           nrow=nrow(L_sigma_1)) %*%
    L_sigma_2
}

#' ordered_transform(x):
#'   for an unconstrained vector x, transform x
#'   such that the elements of x are positive ordered
#'   (i.e., the cumulative sum of exp(x)).
ordered_transform <- function(x) {
  if (is.vector(x) || nrow(x)==1)
    return(cumsum(exp(x)))
  apply(exp(x), MARGIN=2, cumsum)
}

#' sim_metad_condition(N_trials, d_prime, c, log_M, c2_0, c2_1, summarize):
#'   Simulate a set of conditions from the meta-d' model with sensitivity d_prime,
#'   response bias c, metacognitive efficiency log_M, and positive ordered
#'   confidence thresholds c2_0 and c2_1 (for the two responses).
#'
#'   Confidence thresholds are defined relative to meta-c, such that
#'   meta_c2_0 = meta_c - c2_0 and meta_c2_1 = meta_c + c2_1.
#'   Each parameter is a vector, where each index is the parameter value for
#'   a separate condition.
#'
#'   If summarize=FALSE, returns a dataset with one row per observation.
#'   If summarize=TRUE, returns an aggregated dataset where `n` is the
#'   number of observations per condition, response, accuracy, and confidence level.
sim_metad_condition <- function(N_trials=100, d_prime=rep(1, 2), c=rep(0, 2), log_M=rep(0, 2),
                              c2_0_diff=list(rep(.5, 3), rep(.5, 3)),
                              c2_1_diff=list(rep(.5, 3), rep(.5, 3)),
                              metac_absolute=TRUE, summarize=FALSE) {
  tibble(condition=seq_along(d_prime),
         d_prime=d_prime, c=c, log_M=log_M,
         c2_0_diff=c2_0_diff, c2_1_diff=c2_1_diff) |>
    mutate(data=pmap(list(d_prime, c, log_M, c2_0_diff, c2_1_diff), sim_metad,
                     N=N_trials, summarize=summarize, metac_absolute=metac_absolute)) |>
    select(condition, data) |>
    unnest(data)
}


#' sim_metad_participant(N_participants, N_trials, mu_d_prime, sd_d_prime, mu_c, sd_c,
#'                     mu_log_M, sd_log_M, mu_z_c2, sd_z_c2, r_z_c2, summarize):
#'   Simulate a set of participants from the meta-d' model with sensitivity d_prime,
#'   response bias c, metacognitive efficiency log_M, and positive ordered
#'   confidence thresholds c2_0 and c2_1 (for the two responses).
#'
#'   Participant-level parameters are sampled hierarchically according
#'   to normal distributions. The confidence criteria z_c2_0 and z_c2_1
#'   are unconstrained, and the true confidence criteria are generated
#'   using the ordered transform c2_0 = meta_c - cumulative_sum(exp(z_c2_0))
#'   and c2_1 = meta_c + cumulative_sum(exp(z_c2_1)).
#'
#'   If summarize=FALSE, returns a dataset with one row per observation.
#'   If summarize=TRUE, returns an aggregated dataset where `n` is the
#'   number of observations per participant, response, accuracy, and confidence level.
sim_metad_participant <- function(N_participants=100, N_trials=100,
                                  mu_d_prime=1, sd_d_prime=.5, mu_c=0, sd_c=.5,
                                  mu_log_M=0, sd_log_M=.5,
                                  mu_z_c2=rep(-1, 3), sd_z_c2=rep(.1, 3), r_z_c2=diag(3),
                                  metac_absolute=TRUE, summarize=FALSE) {
  expand_grid(participant=1:N_participants) |>
    mutate(d_prime=rnorm(n(), mu_d_prime, sd_d_prime),
           c=rnorm(n(), mu_c, sd_c),
           log_M=rnorm(n(), mu_log_M, sd_log_M),
           c2_0_diff=map(participant,
                         ~ exp(rmvnorm(1, mean=mu_z_c2,
                                       sigma=cov_matrix(sd_z_c2, r_z_c2)))),
           c2_1_diff=map(participant,
                         ~ exp(rmvnorm(1, mean=mu_z_c2,
                                       sigma=cov_matrix(sd_z_c2, r_z_c2)))),
           data=pmap(list(d_prime, c, log_M, c2_0_diff, c2_1_diff), sim_metad,
                     N=N_trials, metac_absolute=metac_absolute, summarize=summarize)) |>
    select(participant, data) |>
    unnest(data)
}

#' sim_metad_participant_condition(N_participants, N_trials, mu_d_prime, sd_d_prime, r_d_prime,
#'                     mu_c, sd_c, r_c, mu_log_M, sd_log_M, r_log_M,
#'                     mu_z_c2, sd_z_c2, r_z_c2_condition, r_z_c2_confidence, summarize):
#'   Simulate a set of participants over multiple within-participant conditions
#'   from the meta-d' model with sensitivity d_prime,
#'   response bias c, metacognitive efficiency log_M, and positive ordered
#'   confidence thresholds c2_0 and c2_1 (for the two responses).
#'
#'   Participant-level parameters are sampled hierarchically according
#'   to normal distributions, except the confidence thresholds which are sampled
#'   according to matrix-normal distributions (i.e., it uses separate correlation matrices
#'   for between-condition and between-confidence-level correlations).
#'   The sampled confidence criteria z_c2_0 and z_c2_1
#'   are unconstrained, and the true confidence criteria are generated
#'   using the ordered transform c2_0 = meta_c - cumulative_sum(exp(z_c2_0))
#'   and c2_1 = meta_c + cumulative_sum(exp(z_c2_1)).
#'
#'   If summarize=FALSE, returns a dataset with one row per observation.
#'   If summarize=TRUE, returns an aggregated dataset where `n` is the
#'   number of observations per participant, response, accuracy, and confidence level.
sim_metad_participant_condition <- function(N_participants=100, N_trials=100,
                                            mu_d_prime=rep(1, 2), sd_d_prime=rep(.5, 2), r_d_prime=diag(2),
                                            mu_c=rep(0, 2), sd_c=rep(.5, 2), r_c=diag(2),
                                            mu_log_M=rep(0, 2), sd_log_M=rep(.5, 2), r_log_M=diag(2),
                                            mu_z_c2=matrix(rep(-1, 6), nrow=3, ncol=2),
                                            sd_z_c2_condition=rep(.1, 2), r_z_c2_condition=diag(2),
                                            sd_z_c2_confidence=rep(.1, 3), r_z_c2_confidence=diag(3),
                                            metac_absolute=TRUE, summarize=FALSE) {
  ## calculate covariance matrices
  sigma_d_prime <- cov_matrix(sd_d_prime, r_d_prime)
  sigma_c <- cov_matrix(sd_c, r_c)
  sigma_log_M <- cov_matrix(sd_log_M, r_log_M)
  L_sigma_z_c2_condition <- chol(cov_matrix(sd_z_c2_condition, r_z_c2_condition))
  L_sigma_z_c2_confidence <- chol(cov_matrix(sd_z_c2_confidence, r_z_c2_confidence))

  expand_grid(participant=1:N_participants,
              condition=seq_along(mu_d_prime)) |>
    group_by(participant) |>
    mutate(d_prime=map_dbl(condition, function(condition, d) d[,condition],
                           rmvnorm(1, mu_d_prime, sigma_d_prime)),
           c=map_dbl(condition, function(condition, c) c[,condition],
                     rmvnorm(1, mu_c, sigma_c)),
           log_M=map_dbl(condition, function(condition, log_m) log_m[,condition],
                         rmvnorm(1, mu_log_M, sigma_log_M)),
           c2_0_diff=map(condition, function(condition, c2) c2[,condition],
                         exp(rmatrixnorm(mu_z_c2,
                                         L_sigma_z_c2_confidence,
                                         L_sigma_z_c2_condition))),
           c2_1_diff=map(condition, function(condition, c2) c2[,condition],
                         exp(rmatrixnorm(mu_z_c2,
                                         L_sigma_z_c2_confidence,
                                         L_sigma_z_c2_condition)))) |>
    ungroup() |>
    mutate(data=pmap(list(d_prime, c, log_M, c2_0_diff, c2_1_diff),
                     sim_metad, N=N_trials, metac_absolute=metac_absolute, summarize=summarize)) |>
    select(participant, condition, data) |>
    unnest(data)
}

