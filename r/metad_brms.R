library(tidyverse)
library(brms)
library(tidybayes)
library(rlang)
library(mvtnorm)
library(abind)

source('metad_utils.R')

## generate a custom brms family for K confidence levels
metad <- function(K) {
  k <- K-1
  custom_family(
    name='metad', 
    dpars=c('mu', 'c', 'M', paste0('metac2zero', 1:k, 'diff'),
            paste0('metac2one', 1:k, 'diff')),
    links=c('identity', 'identity', 'log',
            rep('log', 2*k)),
    lb=c(NA,  NA, 0, rep(0, 2*k)),
    type='int',
    specials=c('multinomial'),
    log_lik=log_lik_metad,
    posterior_predict=posterior_predict_metad,
    posterior_epred=posterior_epred_metad)
}

## write out the stan code for the log likelihood
metad_lpdf <- function(K) {
  k <- K-1
  paste0('real metad_lpmf(array[] int Y, real d_prime, real c, real M, ',
         paste0('real z_meta_c2_0_', 1:k, collapse=', '), ', ',
         paste0('real z_meta_c2_1_', 1:k, collapse=', '),
         ') {
int K = size(Y) %/% 4; // number of confidence levels
	
real meta_d_prime = M * d_prime;
real meta_c = M * c;
vector[K-1] meta_c2_0 = meta_c - cumulative_sum([',
paste0('z_meta_c2_0_', 1:k, collapse=', '),
"]');
vector[K-1] meta_c2_1 = meta_c + cumulative_sum([",
paste0('z_meta_c2_1_', 1:k, collapse=', '),
"]');

// use multinomial likelihood
return multinomial_lpmf( Y[1:(2*K)] | metad_joint_pmf(0, d_prime, c,
																							 meta_d_prime, meta_c, meta_c2_0, meta_c2_1)) +
  multinomial_lpmf(Y[(2*K+1):(4*K)] |  metad_joint_pmf(1, d_prime, c,
																					 meta_d_prime, meta_c, meta_c2_0, meta_c2_1));
}")
}


posterior_epred_metad <- function(prep) {
  post_d_prime <- brms::get_dpar(prep, 'mu')
  post_c <- brms::get_dpar(prep, 'c')
  post_M <- brms::get_dpar(prep, 'M')
  
  # align dimensions
  n_obs <- dim(post_d_prime)[2]
  if (is.vector(post_c))
    post_c <- replicate(n_obs, post_c)
  if (is.vector(post_M))
    post_M <- replicate(n_obs, post_M)
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
  post_d_prime <- brms::get_dpar(prep, 'mu', i = i)
  post_c <- brms::get_dpar(prep, 'c', i = i)
  post_M <- brms::get_dpar(prep, 'M', i = i)
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




fit_metad <- function(formula, data, aggregate=TRUE, .response='N', ...) {
  K <- n_distinct(data$confidence)
  
  # get a list of variables by which to aggregate
  terms <- all.vars(brmsterms(bf(formula, family=metad(K)))$allvars)
  terms <- syms(terms[!(terms %in% c(.response, 'Intercept'))])
  
  # aggregate data by formula terms
  data.aggregated <- data
  if (aggregate) {
    data.aggregated <- metad_aggregate(data, !!!terms, .response=.response)
  }
  
  brm(formula, data.aggregated,
      family=metad(K),
      stanvars=stanvar(scode=read_file('../stan/metad.stanfunctions'),
                       block='functions') +
        stanvar(scode=metad_lpdf(K), block='functions'), ...)
}


################################################################################
#                          Basic metad' model
################################################################################
d <- sim_sdt(N_trials=100000, c=0, log_M=0, c2_0=c(.5, 1), c2_1=c(.5, 1))

m <- fit_metad(bf(N ~ 1, center=FALSE), d,
               init=0, cores=4, backend='cmdstanr',
               prior=set_prior('normal(0, 1)') +
                 set_prior('normal(0, 1)', class='c') +
                 set_prior('lognormal(0, 1)', class='M') +
                 set_prior('lognormal(0, 1)', class='metac2zero1diff') +
                 set_prior('lognormal(0, 1)', class='metac2zero2diff') +
                 set_prior('lognormal(0, 1)', class='metac2one1diff') +
                 set_prior('lognormal(0, 1)', class='metac2one2diff'))
summary(m, prior=TRUE)

predicted_draws(m, newdata=m$data) |>
  group_by(.row, .category) |>
  median_qi(.prediction) |>
  mutate(y=as.integer(m$data$y)) |>
  separate(.category, into=c('var', 'stimulus', 'joint_response'),
           sep='_', convert=TRUE) |>
  mutate(response=factor(as.integer(joint_response > max(joint_response)/2)),
         confidence=factor(ifelse(joint_response > max(joint_response)/2,
                                  joint_response-max(joint_response)/2,
                                  max(joint_response)/2 - joint_response + 1)))|>
  ggplot(aes(x=joint_response)) +
  geom_col(aes(y=y, fill=response, alpha=confidence), ) +
  geom_pointrange(aes(y=.prediction, ymin=.lower, ymax=.upper)) +
  scale_alpha_discrete(range=c(.25, 1)) +
  facet_wrap(~ stimulus, labeller=label_both) +
  theme_classic(18)


epred_draws(m, newdata=m$data) |>
  group_by(.row, .category) |>
  median_qi(.epred) |>
  mutate(.true=response_probabilities(m$data$y[1,])) |>
  separate(.category, into=c('var', 'stimulus', 'joint_response'), sep='_', convert=TRUE) |>
  mutate(response=factor(as.integer(joint_response > max(joint_response)/2)),
         confidence=factor(ifelse(joint_response > max(joint_response)/2,
                                  joint_response-max(joint_response)/2,
                                  max(joint_response)/2 - joint_response + 1))) |>
  #ungroup() |>
  #select(-.row, -var, -.lower, -.upper, -.width, -.point, -.interval)
  ggplot(aes(x=joint_response)) +
  geom_col(aes(y=.true, fill=response, alpha=confidence), ) +
  geom_pointrange(aes(y=.epred, ymin=.lower, ymax=.upper)) +
  scale_alpha_discrete(range=c(.25, 1)) +
  facet_wrap(~ stimulus, labeller=label_both) +
  theme_classic(18)


################################################################################
#                          Condition-level regression
################################################################################
d <- sim_sdt_condition(N_trials=100000,
                       d_prime=c(1, 2), c=c(-1, 1),
                       log_M=c(-3/4, -1/3))

m <- fit_metad(bf(y ~ 0 + Intercept + condition,
                  c ~ 0 + Intercept + condition,
                  M ~ 0 + Intercept + condition,
                  center=FALSE), d,
               ##sample_prior='yes',
               init=0, cores=4, backend='cmdstanr',
               prior=set_prior('normal(0, 1)') +
                 set_prior('normal(0, 1)', dpar='c') +
                 set_prior('normal(0, 1)', dpar='M') +
                 set_prior('lognormal(0, 1)', class='metac2zero1diff') +
                 set_prior('lognormal(0, 1)', class='metac2zero2diff') +
                 set_prior('lognormal(0, 1)', class='metac2zero3diff') +
                 set_prior('lognormal(0, 1)', class='metac2one1diff') +
                 set_prior('lognormal(0, 1)', class='metac2one2diff') +
                 set_prior('lognormal(0, 1)', class='metac2one3diff'))
summary(m, prior=TRUE)

d |>
  distinct(condition) |>
  add_epred_draws(m) |>
  median_qi() |>
  group_by(condition) |>
  mutate(.true=response_probabilities(m$data$y[first(condition),]),
         .category=as.integer(.category),
         condition=factor(condition),
         joint_response=factor((.category-1) %% (n_distinct(d$confidence)*2) + 1),
         stimulus=factor(as.integer(.category > n_distinct(d$confidence)*2))) |>
  ggplot(aes(x=joint_response, fill=condition)) +
  geom_col(aes(y=.true)) +
  geom_pointrange(aes(y=.epred, ymin=.lower, ymax=.upper)) +
  scale_alpha_discrete(range=c(.25, 1)) +
  facet_wrap(condition ~ stimulus, labeller=label_both) +
  theme_classic(18)


linpred_draws(m, newdata=tibble(condition=1:2), dpar=TRUE, transform=TRUE) |>
  pivot_longer(.linpred:metac2one3diff, names_to='.variable', values_to='.value') |>
  group_by(.variable, .add=TRUE) |>
  median_qi()

predicted_draws(m, newdata=select(m$data, condition, y)) |>
  group_by(condition, .row, .category) |>
  median_qi(.prediction) |>
  mutate(y=as.integer(t(m$data$y)),
         stimulus=rep(0:1, each=8, times=2),
         joint_response=rep(1:8, 4)) |>
  ggplot(aes(x=factor(joint_response))) +
  geom_bar(aes(y=y, fill=condition), position=position_dodge(1), stat='identity') +
  geom_pointrange(aes(y=.prediction, ymin=.lower, ymax=.upper, group=condition),
                  position=position_dodge(1)) +
  facet_grid( ~ stimulus, labeller=label_both) +
  theme_classic(18)




################################################################################
#       Condition-level regression (variable confidence threhsolds)
################################################################################
d <- sim_sdt_condition(N_trials=100000,
                       d_prime=c(1, 2), c=c(-1, 1),
                       log_M=c(-3/4, -1/3),
                       c2_0=list(c(.5, 1, 1.5),
                                 c(.25, .5, .75)),
                       c2_1=list(c(1/3, 2/3, 1),
                                 c(1/3, 2/3, 1)))

m <- fit_metad(bf(y ~ 0 + Intercept + condition,
                  c ~ 0 + Intercept + condition,
                  M ~ 0 + Intercept + condition,
                  metac2zero1diff ~ 0 + Intercept + condition,
                  metac2zero2diff ~ 0 + Intercept + condition,
                  metac2zero3diff ~ 0 + Intercept + condition,
                  metac2one1diff ~ 0 + Intercept + condition,
                  metac2one2diff ~ 0 + Intercept + condition,
                  metac2one3diff ~ 0 + Intercept + condition,
                  center=FALSE), d,
               ##sample_prior='yes',
               init=0, cores=4, backend='cmdstanr',
               prior=set_prior('normal(0, 1)') +
                 set_prior('normal(0, 1)', dpar='c') +
                 set_prior('normal(0, 1)', dpar='M') +
                 set_prior('normal(0, 1)', dpar='metac2zero1diff') +
                 set_prior('normal(0, 1)', dpar='metac2zero2diff') +
                 set_prior('normal(0, 1)', dpar='metac2zero3diff') +
                 set_prior('normal(0, 1)', dpar='metac2one1diff') +
                 set_prior('normal(0, 1)', dpar='metac2one2diff') +
                 set_prior('normal(0, 1)', dpar='metac2one3diff'))
summary(m, prior=TRUE)

d |>
  distinct(condition) |>
  add_epred_draws(m) |>
  median_qi() |>
  group_by(condition) |>
  mutate(.true=response_probabilities(m$data$y[first(condition),]),
         .category=as.integer(.category),
         condition=factor(condition),
         joint_response=factor((.category-1) %% (n_distinct(d$confidence)*2) + 1),
         stimulus=factor(as.integer(.category > n_distinct(d$confidence)*2))) |>
  ggplot(aes(x=joint_response, fill=condition)) +
  geom_col(aes(y=.true)) +
  geom_pointrange(aes(y=.epred, ymin=.lower, ymax=.upper)) +
  scale_alpha_discrete(range=c(.25, 1)) +
  facet_wrap(condition ~ stimulus, labeller=label_both) +
  theme_classic(18)

linpred_draws(m, newdata=tibble(condition=1:2), dpar=TRUE, transform=TRUE) |>
  pivot_longer(.linpred:metac2one3diff, names_to='.variable', values_to='.value') |>
  group_by(.variable, .add=TRUE) |>
  median_qi()

predicted_draws(m, newdata=select(m$data, condition, y)) |>
  group_by(condition, .row, .category) |>
  median_qi(.prediction) |>
  mutate(y=as.integer(t(m$data$y)),
         stimulus=rep(0:1, each=8, times=2),
         joint_response=rep(1:8, 4)) |>
  ggplot(aes(x=factor(joint_response))) +
  geom_bar(aes(y=y, fill=condition), position=position_dodge(1), stat='identity') +
  geom_pointrange(aes(y=.prediction, ymin=.lower, ymax=.upper, group=condition),
                  position=position_dodge(1)) +
  facet_grid( ~ stimulus) +
  theme_classic(18)


## fit a model without fixed effects to test main effect
d.summary <- d |> metad_aggregate(condition)
m2 <- fit_metad(bf(y ~ 0 + Intercept,
                   c ~ 0 + Intercept,
                   M ~ 0 + Intercept,
                   metac2zero1diff ~ 0 + Intercept,
                   metac2zero2diff ~ 0 + Intercept,
                   metac2zero3diff ~ 0 + Intercept,
                   metac2one1diff ~ 0 + Intercept,
                   metac2one2diff ~ 0 + Intercept,
                   metac2one3diff ~ 0 + Intercept,
                   center=FALSE), d.summary, aggregate=FALSE,
                ##sample_prior='yes',
                init=0, cores=4, backend='cmdstanr',
                prior=set_prior('normal(0, 1)') +
                  set_prior('normal(0, 1)', dpar='c') +
                  set_prior('normal(0, 1)', dpar='M') +
                  set_prior('normal(0, 1)', dpar='metac2zero1diff') +
                  set_prior('normal(0, 1)', dpar='metac2zero2diff') +
                  set_prior('normal(0, 1)', dpar='metac2zero3diff') +
                  set_prior('normal(0, 1)', dpar='metac2one1diff') +
                  set_prior('normal(0, 1)', dpar='metac2one2diff') +
                  set_prior('normal(0, 1)', dpar='metac2one3diff'))
summary(m2, prior=TRUE)

d |>
  distinct(condition) |>
  add_epred_draws(m2) |>
  median_qi() |>
  group_by(condition) |>
  mutate(.true=response_probabilities(m$data$y[first(condition),]),
         .category=as.integer(.category),
         condition=factor(condition),
         joint_response=factor((.category-1) %% (n_distinct(d$confidence)*2) + 1),
         stimulus=factor(as.integer(.category > n_distinct(d$confidence)*2))) |>
  ggplot(aes(x=joint_response, fill=condition)) +
  geom_col(aes(y=.true)) +
  geom_pointrange(aes(y=.epred, ymin=.lower, ymax=.upper)) +
  scale_alpha_discrete(range=c(.25, 1)) +
  facet_wrap(condition ~ stimulus, labeller=label_both) +
  theme_classic(18)

loo(m, m2)


################################################################################
#                          Participant-level regression
################################################################################
d <- sim_sdt_participant(N_participants=500, N_trials=250,
                         mu_d_prime=1, sd_d_prime=.5,
                         mu_c=0, sd_c=.5,
                         mu_log_M=0, sd_log_M=.5,
                         mu_z_c2=rep(-1, 3), sd_z_c2=rep(.1, 3),
                         r_z_c2=corr_matrix(.5, nrow=3))

m <- fit_metad(bf(y ~ 0 + Intercept + (1 | participant),
                  c ~ 0 + Intercept + (1 | participant),
                  M ~ 0 + Intercept + (1 | participant),
                  metac2zero1diff ~  0 + Intercept + (1 |p0| participant),
                  metac2zero2diff ~ 0 + Intercept + (1 |p0| participant),
                  metac2zero3diff ~ 0 + Intercept + (1 |p0| participant),
                  metac2one1diff ~ 0 + Intercept + (1 |p1| participant),
                  metac2one2diff ~ 0 + Intercept + (1 |p1| participant),
                  metac2one3diff ~ 0 + Intercept + (1 |p1| participant),
                  center=FALSE), d,
               init=0, cores=4, backend='cmdstanr',
               prior=set_prior('normal(0, 1)') +
                 set_prior('normal(0, 1)', dpar='c') +
                 set_prior('normal(0, 1)', dpar='M') +
                 set_prior('normal(0, 1)', dpar='metac2zero1diff') +
                 set_prior('normal(0, 1)', dpar='metac2zero2diff') +
                 set_prior('normal(0, 1)', dpar='metac2zero3diff') +
                 set_prior('normal(0, 1)', dpar='metac2one1diff') +
                 set_prior('normal(0, 1)', dpar='metac2one2diff') +
                 set_prior('normal(0, 1)', dpar='metac2one3diff'))
summary(m, prior=TRUE)



draws <- d |>
  distinct(participant) |>
  add_linpred_draws(m, dpar=TRUE, transform=TRUE)

draws |>
  mutate(meta_c=M*c,
         meta_c2_0_1=meta_c - metac2zero1diff,
         meta_c2_0_2=meta_c2_0_1 - metac2zero2diff,
         meta_c2_0_3=meta_c2_0_2 - metac2zero3diff,
         meta_c2_1_1=meta_c + metac2one1diff,
         meta_c2_1_2=meta_c2_1_1 + metac2one2diff,
         meta_c2_1_3=meta_c2_1_2 + metac2one3diff) |>
  select(-ends_with('diff'), -mu, -meta_c) |>
  rename(d_prime=.linpred) |>
  median_qi(.simple_names=FALSE) |>
  rename_with(~ paste0(., '.median'),
              .cols=c(d_prime, c, M,
                      matches('meta_c2_[[:digit:]]_[[:digit:]]$'))) |>
  left_join(d |> group_by(participant) |>
              reframe(across(d_prime:meta_c2_1, first)) |>
              select(-meta_d_prime) |>
              group_by(participant) |>
              mutate(conf=row_number()) |>
              pivot_wider(names_from=conf, values_from=meta_c2_0:meta_c2_1) |>
              rename_with(~ paste0(., '.true'), d_prime:meta_c2_1_3)) |>
  select(-.width, -.point, -.interval) |>
  pivot_longer(d_prime.median:meta_c2_1_3.true, names_sep='\\.',
               names_to=c('.variable', '.estimate')) |>
  arrange(participant, .variable, .estimate) |>
  pivot_wider(names_from=.estimate) |>
  ggplot(aes(y=median, x=true)) +
  geom_abline(linetype='dashed', slope=1, intercept=0) +
  geom_pointrange(aes(ymin=lower, ymax=upper)) +
  facet_wrap(~ .variable) +
  ##coord_fixed() +
  theme_bw()
