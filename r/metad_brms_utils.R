## generate a custom brms family for K confidence levels
metad <- function(K) {
  k <- K-1
  custom_family(
    name='metad', 
    dpars=c('mu', 'dprime', 'c', paste0('metac2zero', 1:k, 'diff'),
            paste0('metac2one', 1:k, 'diff')),
    links=c('log', 'identity', 'identity',
            rep('log', 2*k)),
    lb=c(0, NA,  NA, rep(0, 2*k)),
    type='int',
    specials=c('multinomial'),
    log_lik=log_lik_metad,
    posterior_predict=posterior_predict_metad,
    posterior_epred=posterior_epred_metad)
}

## write out the stan code for the log likelihood
metad_lpdf <- function(K, distribution='std_normal') {
  k <- K-1

  paste0("	// Convert a binary int x from {0, 1} to {-1, 1}
	int to_signed(int x) {
	  return 2*x - 1;
	}

	// P(response, confidence | stimulus) given as simplex
	// [P(resp=0, conf=K), .... P(resp=0, conf=1), P(resp=1, conf=1), ... P(resp=1, conf=K)]
	vector metad_", distribution, "_pmf(int stimulus, real d_prime, real c, real meta_d_prime, real meta_c, vector meta_c2_0, vector meta_c2_1) {
		// number of confidence levels
		int K = size(meta_c2_0)+1;
    
  	// type-1 response probabilities
	  real lp_1 = ", distribution, "_lcdf(to_signed(stimulus)*d_prime/2 - c);
  	real lp_0 = ", distribution, "_lccdf(to_signed(stimulus)*d_prime/2 - c);

  	// means of type-2 distributions
  	real meta_mu_1 = to_signed(stimulus) * meta_d_prime/2;
  	real meta_mu_0 = -meta_mu_1;

	  vector[K] lp2_1;         // CDFs (response == 1)
  	vector[K] lp2_0;         // CDFs (response == 0)
		vector[2*K] log_theta;   // joint (type-1 x type-2) response probabilities

	  lp2_1[1] = ", distribution, "_lcdf(meta_mu_1 - meta_c);
  	lp2_0[1] = ", distribution, "_lcdf(meta_mu_0 + meta_c);
  	for (k in 2:K) {
    	lp2_1[k] = ", distribution, "_lcdf(meta_mu_1 - meta_c2_1[k-1]);
    	lp2_0[k] = ", distribution, "_lcdf(meta_mu_0 + meta_c2_0[k-1]);

			log_theta[K-k+2] = log_diff_exp(lp2_0[k-1], lp2_0[k]);
    	log_theta[K+k-1] = log_diff_exp(lp2_1[k-1], lp2_1[k]);
  	}
  	log_theta[1] = lp2_0[K];
  	log_theta[2*K] = lp2_1[K];

	  // weight by P(response|stimulus) and normalize
  	log_theta[1:K] += lp_0 - lp2_0[1];
  	log_theta[(K+1):(2*K)] += lp_1 - lp2_1[1];

	  return exp(log_theta);
	}

	real metad_lpmf(array[] int Y, real M, real d_prime, real c, ",
  paste0("real z_meta_c2_0_", 1:k, collapse=", "), ", ",
  paste0("real z_meta_c2_1_", 1:k, collapse=", "),
  ") {
		int K = size(Y) %/% 4; // number of confidence levels

		real meta_d_prime = M * d_prime;
		real meta_c = M * c;
		vector[K-1] meta_c2_0 = meta_c - cumulative_sum([",
  paste0("z_meta_c2_0_", 1:k, collapse=", "),
  "]');
		vector[K-1] meta_c2_1 = meta_c + cumulative_sum([",
  paste0("z_meta_c2_1_", 1:k, collapse=", "),
  "]');

		// use multinomial likelihood
		return multinomial_lpmf(Y[1:(2*K)] | metad_", distribution, "_pmf(0, d_prime, c,
														meta_d_prime, meta_c, meta_c2_0, meta_c2_1)) +
  		multinomial_lpmf(Y[(2*K+1):(4*K)] |  metad_", distribution, "_pmf(1, d_prime, c,
											 meta_d_prime, meta_c, meta_c2_0, meta_c2_1));
	}")
}

posterior_epred_metad <- function(prep) {
  post_M <- brms::get_dpar(prep, 'mu')
  post_d_prime <- brms::get_dpar(prep, 'dprime')
  post_c <- brms::get_dpar(prep, 'c')
  
  # align dimensions
  n_obs <- dim(post_M)[2]
  if (is.vector(post_d_prime))
    post_d_prime <- replicate(n_obs, post_d_prime)
  if (is.vector(post_c))
    post_c <- replicate(n_obs, post_c)
  post_meta_d_prime <- post_M * post_d_prime
  post_meta_c <- post_M * post_c
  
  # determine confidence thresholds
  dpars <- names(prep$dpars)
  post_meta_c2_0 <- NULL
  post_meta_c2_1 <- NULL

  if (is.vector(brms::get_dpar(prep, 'metac2zero1diff'))) {
    post_meta_c2_0 <- dpars[str_detect(dpars, 'metac2zero')] |>
      sapply(function(s) brms::get_dpar(prep, s)) |>
      apply(1, cumsum) |>
      t() |>
      replicate(last(dim(post_meta_c)), expr=_) |>
      aperm(c(1, 3, 2))
  } else {
    post_meta_c2_0 <- dpars[str_detect(dpars, 'metac2zero')] |>
      lapply(function(s) brms::get_dpar(prep, s)) |>
      abind(along=3) |>
      apply(1:2, cumsum) |>
      aperm(c(2, 3, 1))
  }

  if (is.vector(brms::get_dpar(prep, 'metac2one1diff'))) {
    post_meta_c2_1 <- dpars[str_detect(dpars, 'metac2one')] |>
      sapply(function(s) brms::get_dpar(prep, s)) |>
      apply(1, cumsum) |>
      t() |>
      replicate(last(dim(post_meta_c)), expr=_) |>
      aperm(c(1, 3, 2))
  } else {    
    post_meta_c2_1 <- dpars[str_detect(dpars, 'metac2one')] |>
      lapply(function(s) brms::get_dpar(prep, s)) |>
      abind(along=3) |>
      apply(1:2, cumsum) |>
      aperm(c(2, 3, 1))
  }
  
  # calculate number of confidence thresholds
  k <- last(dim(post_meta_c2_0))
  K <- k+1
  
  # calculate confidence threhsolds
  post_meta_c2_0 <- replicate(k, post_meta_c) - post_meta_c2_0
  post_meta_c2_1 <- replicate(k, post_meta_c) + post_meta_c2_1
  
  # calculate joint response & confidence probabilities
  p <- array(dim=c(dim(post_d_prime), 4*K))
  for (s in 1:first(dim(post_d_prime))) {
    for (i in 1:last(dim(post_d_prime))) {
      p[s,i,1:(2*K)] <-
        sdt_joint_pmf(0, post_d_prime[s,i], post_c[s,i], post_meta_d_prime[s,i],
                      post_meta_c[s,i], post_meta_c2_0[s,i,], post_meta_c2_1[s,i,])
      p[s,i,(2*K+1):(4*K)] <-
        sdt_joint_pmf(1, post_d_prime[s,i], post_c[s,i], post_meta_d_prime[s,i],
                      post_meta_c[s,i], post_meta_c2_0[s,i,], post_meta_c2_1[s,i,])
    }
  }

  p
}

lp_metad <- function(i, prep) {
  post_M <- brms::get_dpar(prep, 'mu', i = i)
  post_d_prime <- brms::get_dpar(prep, 'dprime', i = i)
  post_c <- brms::get_dpar(prep, 'c', i = i)
  post_meta_d_prime <- post_M * post_d_prime
  post_meta_c <- post_M * post_c

  # determine confidence thresholds
  dpars <- names(prep$dpars)
  post_meta_c2_0 <- dpars[str_detect(dpars, 'metac2zero')] |>
    sapply(function(s) brms::get_dpar(prep, s, i=i)) |>
    apply(1, cumsum) |>
    t()
  post_meta_c2_1 <- dpars[str_detect(dpars, 'metac2one')] |>
    sapply(function(s) brms::get_dpar(prep, s, i=i)) |>
    apply(1, cumsum) |>
    t()
  post_meta_c2_0 <- post_meta_c - post_meta_c2_0
  post_meta_c2_1 <- post_meta_c + post_meta_c2_1
  post_meta_c2_0 <- split(post_meta_c2_0, row(post_meta_c2_0))
  post_meta_c2_1 <- split(post_meta_c2_1, row(post_meta_c2_1))
  
  # calculate joint response & confidence probabilities
  PMF <- Vectorize(sdt_joint_pmf)
  lp_0 <- PMF(0, post_d_prime, post_c, post_meta_d_prime, post_meta_c,
              post_meta_c2_0, post_meta_c2_1, log=TRUE)
  lp_1 <- PMF(1, post_d_prime, post_c, post_meta_d_prime, post_meta_c,
              post_meta_c2_0, post_meta_c2_1, log=TRUE)

  t(rbind(lp_0, lp_1))
}

log_lik_metad <- function(i, prep) {
  p <- exp(lp_metad(i, prep))

  if (any(is.na(prep$data$Y))) {
    stop('Error: please provide sample data y with trial counts')
  }
  
  y <- prep$data$Y[i,]
  N_0 <- sum(y[1:(length(y)/2)])
  N_1 <- sum(y[(length(y)/2+1):length(y)])

  # calculate multinomial response probabilities
  apply(p[,1:(ncol(p)/2)], 1,
        function(prob) dmultinom(y[1:(length(y)/2)],
                                 size=N_0, prob=prob, log=TRUE)) +
    apply(p[,(ncol(p)/2+1):ncol(p)], 1,
          function(prob) dmultinom(y[(length(y)/2+1):length(y)],
                                   size=N_1, prob=prob, log=TRUE))
}

posterior_predict_metad <- function(i, prep, ...) {
  p <- exp(lp_metad(i, prep))

  if (any(is.na(prep$data$Y))) {
    stop('Error: please provide sample data y with trial counts')
  }
  
  y <- prep$data$Y[i,]
  N_0 <- sum(y[1:(length(y)/2)])
  N_1 <- sum(y[(length(y)/2+1):length(y)])
  
  # simulate from a multinomial distribution
  rbind(apply(p[,1:(ncol(p)/2)], 1, rmultinom, n=1, size=N_0),
        apply(p[,(ncol(p)/2+1):ncol(p)], 1, rmultinom, n=1, size=N_1)) |>
    t()
}

#' metad_aggregate(data, ..., .response='y'):
#'   Aggregate `data` by columns `response`, `confidence`,
#'   and any other variables in `...`.
#'
#'   Arguments:
#'     `data`: the tibble to aggregate
#'     `...`: grouping columns in the tibble
#'     `.response`: the name of the column containing trial counts
#'
#'   Value:
#'     A tibble with one row per combination of the variables in `...`,
#'     and another column named by the value of `.response` containing trial counts.
#'     For K confidence levels, this will be an N x K*4 matrix, such that the
#'     columns represent:
#'    [N(stimulus==0, confidence==K), ..., N(stimulus==0, confidence==1),
#'     N(stimulus==0, confidence==1), ..., N(stimulus==0, confidence==K),
#'     N(stimulus==1, confidence==K), ..., N(stimulus==1, confidence==1),
#'     N(stimulus==1, confidence==1), ..., N(stimulus==1, confidence==K)]
metad_aggregate <- function(data, ..., .response='N') {
  # number of confidence levels
  K <- n_distinct(data$confidence)

  data <- data |>
    mutate(joint_response=factor(joint_response(response, confidence, K)),
           stimulus=factor(stimulus),
           across(c(...), factor)) |>
    group_by(...) |>
    count(stimulus, joint_response, .drop=FALSE) |>
    pivot_wider(names_from=c(stimulus, joint_response),
                values_from=n, names_prefix=glue::glue('{.response}_')) %>%
    mutate('{.response}_0' := sum(c_across(starts_with(glue::glue('{.response}_0_'),
                                                       ignore.case=FALSE))),
           '{.response}_1' := sum(c_across(starts_with(glue::glue('{.response}_1_'),
                                                       ignore.case=FALSE))))
  # convert counts into a matrix-column
  tibble(select(data, ...,
                all_of(c(glue::glue('{.response}_0'),
                         glue::glue('{.response}_1')))),
         '{.response}' := data |> ungroup() |>
           select(matches(glue::glue('{.response}_([[:digit:]]+)_([[:digit:]]+)'),
                          ignore.case=FALSE)) |>
           as.matrix())
}


fit_metad <- function(formula, data, ..., aggregate=TRUE,
                      distribution='std_normal', stanvars=NULL) {
  K <- NULL
  data.aggregated <- NULL

  # determine response variable
  .response <- all.vars(formula$formula)[attr(terms(formula$formula), 'response')]
  
  # aggregate data by formula terms
  if (aggregate) {
    K <- n_distinct(data$confidence)
    
    # get a list of variables by which to aggregate
    terms <- all.vars(brmsterms(bf(formula, family=metad(K)))$allvars)
    terms <- syms(terms[!(terms %in% c(.response, 'Intercept'))])
    data.aggregated <- metad_aggregate(data, !!!terms, .response=.response)
  } else {
    K <- ncol(pull(data, .response)) / 4
    data.aggregated <- data
  }

  # add metad stanvars to any user-defined stanvars
  sv <- stanvar(scode=metad_lpdf(K, distribution=distribution), block='functions')
  if (!is.null(stanvars))
    sv <- sv + stanvars
  
  brm(formula, data.aggregated, family=metad(K), stanvars=sv, ...)
}




mean_confidence_draws <- function(object, newdata, ...) {
  draws <- epred_draws(object, newdata, ...)

  ## number of confidence levels
  K <- as.integer(n_distinct(draws$.category) / 4)

  ## grouping columns
  .cols <- names(newdata)
  .cols <- .cols[!(.cols %in% c('.row', 'stimulus', '.draw'))]
  
  
  draws |>
    mutate(.category=as.integer(.category),
           stimulus=as.integer(.category > 2*K),
           joint_response=ifelse(stimulus, .category - 2*K, .category),
           response=as.integer(joint_response > K),
           confidence=ifelse(joint_response > K,
                             joint_response - K, K + 1 - joint_response)) |>
    group_by(.row, stimulus, !!!syms(.cols), .draw) |>
    summarize(.epred=sum(.epred*confidence))
}

add_mean_confidence_draws <- function(newdata, object, ...) {
  mean_confidence_draws(object, newdata, ...)
}




roc1_draws <- function(object, newdata, ..., bounds=FALSE) {
  draws <- epred_draws(object=object, newdata=newdata, ...)

  ## number of confidence levels
  K <- as.integer(n_distinct(draws$.category) / 4)

  ## grouping columns
  .cols <- names(newdata)
  .cols <- .cols[!(.cols %in% c('.row', 'joint_response', 'response', 'confidence'))]
  
  draws <- draws |>
    mutate(.category=as.integer(.category),
           stimulus=as.integer(.category > 2*K),
           joint_response=ifelse(stimulus, .category - 2*K, .category),
           response=as.integer(joint_response > K),
           confidence=ifelse(joint_response > K,
                             joint_response - K, K + 1 - joint_response)) |>
    filter(joint_response < 2*K) |>
    group_by(.row, stimulus, .draw) |>
    mutate(.epred=1-cumsum(.epred)) |>
    select(-.category) |>
    pivot_wider(names_from=stimulus, values_from=.epred, names_prefix='p_') |>
    rename(p_hit=p_1, p_fa=p_0) |>
    group_by(.row, !!!syms(.cols), joint_response, response, confidence)
  
  if (bounds) {
    ## add (0, 0) and (1, 1) points to ROC
    draws <- draws |>
      bind_rows(draws |>
                  ungroup() |>
                  distinct(.row, !!!syms(.cols), .chain, .iteration, .draw) |>
                  expand_grid(tibble(joint_response=c(0, K*2),
                                     response=c(0, 1),
                                     confidence=c(K+1, K),
                                     p_fa=c(1, 0),
                                     p_hit=c(1, 0)))) |>
      arrange(joint_response)
  }

  draws
}

add_roc1_draws <- function(newdata, object, ..., bounds=FALSE) {
  roc1_draws(object, newdata, ..., bounds=bounds)
}

roc2_draws <- function(object, newdata, ..., bounds=FALSE) {
  draws <- epred_draws(object=object, newdata=newdata, ...)

  ## number of confidence levels
  K <- as.integer(n_distinct(draws$.category) / 4)
  
  ## grouping columns
  .cols <- names(newdata)
  .cols <- .cols[!(.cols %in% c('.row', 'response', 'confidence'))]
  
  draws <- draws |>
    mutate(.category=as.integer(.category),
           stimulus=as.integer(.category > 2*K),
           joint_response=ifelse(stimulus, .category - 2*K, .category),
           response=as.integer(joint_response > K),
           accuracy=as.integer(stimulus == response),
           confidence=ifelse(joint_response > K,
                             joint_response - K, K + 1 - joint_response)) |>
    group_by(.row, .draw, accuracy, response) |>
    mutate(.epred=cumsum(.epred) / sum(.epred),
           .epred=ifelse(response, 1-.epred, .epred)) |>
    filter(!(response==0 & confidence == 1),
           !(response==1 & confidence == K)) |>
    select(-.category, -joint_response, -stimulus) |>
    ungroup() |>
    pivot_wider(names_from=accuracy, values_from=.epred, names_prefix='p_') |>
    rename(p_hit2=p_1, p_fa2=p_0) |>
    group_by(.row, response, confidence, !!!syms(.cols))

  if (bounds) {
    ## add (0, 0) and (1, 1) points to ROC
    draws <- draws |>
      bind_rows(draws |>
                  ungroup() |>
                  distinct(.row, !!!syms(.cols), .chain, .iteration, .draw) |>
                  expand_grid(tibble(response=c(0, 0, 1, 1),
                                     confidence=c(1, K+1, 0, K),
                                     p_fa2=c(1, 0, 1, 0),
                                     p_hit2=c(1, 0, 1, 0))))
  }

  draws
}

add_roc2_draws <- function(newdata, object, ..., bounds=FALSE) {
  roc2_draws(object, newdata, ..., bounds=bounds)
}
