library(tidyverse)
library(brms)
library(tidybayes)
library(rlang)
library(mvtnorm)
library(abind)

source('metad_utils.R')
source('metad_brms_utils.R')


################################################################################
#                          Basic metad' model
################################################################################
d <- sim_metad(N_trials=100000, d_prime=1, c=1, log_M=-1,
               c2_0=c(.5, 1), c2_1=c(.5, 1))

m <- fit_metad(bf(N ~ 0 + Intercept), data=d,
               cores=4, backend='cmdstanr',
               prior=prior(normal(0, 1)) +
                 prior(normal(0, 1), class=dprime) +
                 prior(normal(0, 1), class=c) +
                 prior(lognormal(0, 1), class=metac2zero1diff) +
                 prior(lognormal(0, 1), class=metac2zero2diff) +
                 prior(lognormal(0, 1), class=metac2one1diff) +
                 prior(lognormal(0, 1), class=metac2one2diff))
summary(m, prior=TRUE)

## posterior predictions (counts of type 1/type 2 responses)
predicted_draws(m, metad_aggregate(d)) |>
  group_by(.row, .category) |>
  median_qi(.prediction) |>
  mutate(N=as.integer(m$data$N)) |>
  separate(.category, into=c('var', 'stimulus', 'joint_response'),
           sep='_', convert=TRUE) |>
  mutate(response=factor(as.integer(joint_response > max(joint_response)/2)),
         confidence=factor(ifelse(joint_response > max(joint_response)/2,
                                  joint_response-max(joint_response)/2,
                                  max(joint_response)/2 - joint_response + 1)))|>
  ggplot(aes(x=joint_response)) +
  geom_col(aes(y=N, fill=response, alpha=confidence), ) +
  geom_pointrange(aes(y=.prediction, ymin=.lower, ymax=.upper)) +
  scale_alpha_discrete(range=c(.25, 1)) +
  facet_wrap(~ stimulus, labeller=label_both) +
  theme_classic(18)

## posterior predictions (probabilities of type 1/type 2 responses)
epred_draws(m, newdata=tibble(.row=1)) |>
  group_by(.row, .category) |>
  median_qi(.epred) |>
  mutate(.true=response_probabilities(m$data$N[1,])) |>
  separate(.category, into=c('var', 'stimulus', 'joint_response'), sep='_', convert=TRUE) |>
  mutate(response=factor(as.integer(joint_response > max(joint_response)/2)),
         confidence=factor(ifelse(joint_response > max(joint_response)/2,
                                  joint_response-max(joint_response)/2,
                                  max(joint_response)/2 - joint_response + 1))) |>
  ggplot(aes(x=joint_response)) +
  geom_col(aes(y=.true, fill=response, alpha=confidence), ) +
  geom_pointrange(aes(y=.epred, ymin=.lower, ymax=.upper)) +
  scale_alpha_discrete(range=c(.25, 1)) +
  facet_wrap(~ stimulus, labeller=label_both) +
  theme_classic(18)

## calculate mean confidence per stimulus
tibble(.row=1) |>
  add_mean_confidence_draws(m) |>
  median_qi(.epred) |>
  left_join(d |>
              group_by(stimulus, response, confidence) |>
              summarize(across(theta_1:theta_2, first)) |>
              group_by(stimulus) |>
              summarize(.true=sum(theta_1 * theta_2 * confidence)))

# psuedo type-1 ROC
tibble(.row=1) |>
  add_roc1_draws(m, bounds=FALSE) |>
  median_qi(p_fa, p_hit) |>
  ggplot(aes(x=p_fa, xmin=p_fa.lower, xmax=p_fa.upper,
             y=p_hit, ymin=p_hit.lower, ymax=p_hit.upper)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_errorbar(orientation='y', width=.01) +
  geom_errorbar(orientation='x', width=.01) +
  geom_line() +
  coord_fixed(xlim=0:1, ylim=0:1, expand=FALSE) +
  xlab('P(False Alarm)') + ylab('P(Hit)') +
  theme_bw(18)

# type 2 ROC
tibble(.row=1) |>
  add_roc2_draws(m, bounds=TRUE) |>
  median_qi(p_hit2, p_fa2) |>
  mutate(response=factor(response)) |>
  ggplot(aes(x=p_fa2, xmin=p_fa2.lower, xmax=p_fa2.upper,
             y=p_hit2, ymin=p_hit2.lower, ymax=p_hit2.upper,
             color=response)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_errorbar(orientation='y', width=.01) +
  geom_errorbar(orientation='x', width=.01) +
  geom_line() +
  coord_fixed(xlim=0:1, ylim=0:1, expand=FALSE) +
  xlab('P(Type 2 False Alarm)') + ylab('P(Type 2 Hit)') +
  theme_bw(18)

# model parameters
tibble(.row=1) |>
  add_linpred_draws(m, dpar=TRUE, transform=TRUE) |>
  mutate(meta_dprime=mu * dprime,
         across(starts_with('metac2'), ~ . / meta_dprime))

################################################################################
#                    Basic meta-d' model w/ Gumbel-min noise
################################################################################
gumbel_min <- stanvar(scode='
real gumbel_min_lccdf(real y) {
  return -exp(y);
}
real gumbel_min_lcdf(real y) {
  return log1m_exp(-exp(y));
}
', block='functions')

g <- fit_metad(bf(M ~ 0 + Intercept), data=d,
               cores=4, backend='cmdstanr',
               prior=prior(normal(0, 1)) +
                 prior(normal(0, 1), class=dprime) +
                 prior(normal(0, 1), class=c) +
                 prior(lognormal(0, 1), class=metac2zero1diff) +
                 prior(lognormal(0, 1), class=metac2zero2diff) +
                 prior(lognormal(0, 1), class=metac2one1diff) +
                 prior(lognormal(0, 1), class=metac2one2diff),
               stanvars=gumbel_min, distribution='gumbel_min')
summary(g, prior=TRUE)

epred_draws(g, newdata=tibble(.row=1)) |>
  group_by(.row, .category) |>
  median_qi(.epred) |>
  mutate(.true=response_probabilities(g$data$M[1,])) |>
  separate(.category, into=c('var', 'stimulus', 'joint_response'),
           sep='_', convert=TRUE) |>
  mutate(response=factor(as.integer(joint_response > max(joint_response)/2)),
         confidence=factor(ifelse(joint_response > max(joint_response)/2,
                                  joint_response-max(joint_response)/2,
                                  max(joint_response)/2 - joint_response + 1))) |>
  ggplot(aes(x=joint_response)) +
  geom_col(aes(y=.true, fill=response, alpha=confidence), ) +
  geom_pointrange(aes(y=.epred, ymin=.lower, ymax=.upper)) +
  scale_alpha_discrete(range=c(.25, 1)) +
  facet_wrap(~ stimulus, labeller=label_both) +
  theme_classic(18)

## calculate mean confidence per stimulus
epred_draws(g, tibble(.row=1)) |>
  group_by(.draw) |>
  mutate(stimulus=rep(0:1, each=n_distinct(d$confidence)*2)) |>
  group_by(stimulus, .draw) |>
  summarize(.epred=list(.epred)) |>
  mutate(.mean=map_dbl(.epred, `%*%`,
                       c(n_distinct(d$confidence):1,
                         1:(n_distinct(d$confidence))))) |>
  median_qi(.mean) |>
  left_join(d |>
              group_by(stimulus, response, confidence) |>
              summarize(across(theta_1:theta_2, first)) |>
              group_by(stimulus) |>
              summarize(.true=sum(theta_1 * theta_2 * confidence)))


## calculate mean confidence per stimulus
tibble(.row=1) |>
  add_mean_confidence_draws(g) |>
  median_qi(.epred) |>
  left_join(d |>
              group_by(stimulus, response, confidence) |>
              summarize(across(theta_1:theta_2, first)) |>
              group_by(stimulus) |>
              summarize(.true=sum(theta_1 * theta_2 * confidence)))


# psuedo type-1 ROC
tibble(.row=1) |>
  add_roc1_draws(g, bounds=TRUE) |>
  median_qi(p_fa, p_hit) |>
  ggplot(aes(x=p_fa, xmin=p_fa.lower, xmax=p_fa.upper,
             y=p_hit, ymin=p_hit.lower, ymax=p_hit.upper)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_errorbar(orientation='y', width=.01) +
  geom_errorbar(orientation='x', width=.01) +
  geom_line() +
  coord_fixed(xlim=0:1, ylim=0:1, expand=FALSE) +
  xlab('P(False Alarm)') + ylab('P(Hit)') +
  theme_bw(18)

# type 2 ROC
roc2_draws(g, tibble(.row=1), bounds=TRUE) |>
  median_qi(p_hit2, p_fa2) |>
  mutate(response=factor(response)) |>
  ggplot(aes(x=p_fa2, xmin=p_fa2.lower, xmax=p_fa2.upper,
             y=p_hit2, ymin=p_hit2.lower, ymax=p_hit2.upper,
             color=response)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_errorbar(orientation='y', width=.01) +
  geom_errorbar(orientation='x', width=.01) +
  geom_line() +
  coord_fixed(xlim=0:1, ylim=0:1, expand=FALSE) +
  xlab('P(Type 2 False Alarm)') + ylab('P(Type 2 Hit)') +
  theme_bw(18)


################################################################################
#                          Condition-level regression
################################################################################
d <- sim_metad_condition(N_trials=100000,
                         d_prime=c(1, 2), c=c(-1, 1),
                         log_M=c(-3/4, -1/3))

m <- fit_metad(bf(N ~ 0 + Intercept + condition,
                  dprime+c ~ 0 + Intercept + condition),
               data=d, cores=4, backend='cmdstanr',
               prior=prior(normal(0, 1)) +
                 prior(normal(0, 1), dpar=dprime) +
                 prior(normal(0, 1), dpar=c) +
                 prior(lognormal(0, 1), class=metac2zero1diff) +
                 prior(lognormal(0, 1), class=metac2zero2diff) +
                 prior(lognormal(0, 1), class=metac2zero3diff) +
                 prior(lognormal(0, 1), class=metac2one1diff) +
                 prior(lognormal(0, 1), class=metac2one2diff) +
                 prior(lognormal(0, 1), class=metac2one3diff))
summary(m, prior=TRUE)

d |>
  distinct(condition) |>
  add_epred_draws(m) |>
  median_qi() |>
  group_by(condition) |>
  mutate(.true=response_probabilities(m$data$N[first(condition),]),
         .category=as.integer(.category),
         condition=factor(condition),
         joint_response=factor((.category-1) %% (n_distinct(d$confidence)*2) + 1),
         stimulus=factor(as.integer(.category > n_distinct(d$confidence)*2))) |>
  ggplot(aes(x=joint_response, fill=condition)) +
  geom_col(aes(y=.true)) +
  geom_pointrange(aes(y=.epred, ymin=.lower, ymax=.upper)) +
  facet_wrap(condition ~ stimulus, labeller=label_both) +
  theme_classic(18)


linpred_draws(m, newdata=tibble(condition=1:2), dpar=TRUE, transform=TRUE) |>
  pivot_longer(.linpred:metac2one3diff, names_to='.variable', values_to='.value') |>
  group_by(.variable, .add=TRUE) |>
  median_qi()

predicted_draws(m, newdata=select(m$data, condition, N)) |>
  group_by(condition, .row, .category) |>
  median_qi(.prediction) |>
  mutate(y=as.integer(t(m$data$N)),
         stimulus=rep(0:1, each=8, times=2),
         joint_response=rep(1:8, 4)) |>
  ggplot(aes(x=factor(joint_response))) +
  geom_bar(aes(y=y, fill=condition), position=position_dodge(1), stat='identity') +
  geom_pointrange(aes(y=.prediction, ymin=.lower, ymax=.upper, group=condition),
                  position=position_dodge(1)) +
  facet_grid( ~ stimulus, labeller=label_both) +
  theme_classic(18)


## calculate mean confidence per stimulus
d |>
  distinct(condition) |>
  add_mean_confidence_draws(m) |>
  median_qi(.epred) |>
  left_join(d |>
              group_by(condition, stimulus, response, confidence) |>
              summarize(across(theta_1:theta_2, first)) |>
              group_by(condition, stimulus) |>
              summarize(.true=sum(theta_1 * theta_2 * confidence)))


# psuedo type-1 ROC
d |>
  distinct(condition) |>
  mutate(condition=factor(condition)) |>
  add_roc1_draws(m, bounds=TRUE) |>
  group_by(condition, .add=TRUE) |>
  median_qi(p_fa, p_hit) |>
  ggplot(aes(x=p_fa, xmin=p_fa.lower, xmax=p_fa.upper,
             y=p_hit, ymin=p_hit.lower, ymax=p_hit.upper,
             group=condition, color=condition)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_errorbar(orientation='y', width=.01) +
  geom_errorbar(orientation='x', width=.01) +
  geom_line() +
  coord_fixed(xlim=0:1, ylim=0:1, expand=FALSE) +
  xlab('P(False Alarm)') + ylab('P(Hit)') +
  theme_bw(18)

# type 2 ROC
d |>
  distinct(condition) |>
  mutate(condition=factor(condition)) |>
  add_roc2_draws(m, bounds=TRUE) |>
  median_qi(p_hit2, p_fa2) |>
  mutate(response=factor(response)) |>
  ggplot(aes(x=p_fa2, xmin=p_fa2.lower, xmax=p_fa2.upper,
             y=p_hit2, ymin=p_hit2.lower, ymax=p_hit2.upper,
             color=response, linetype=condition)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_errorbar(orientation='y', width=.01) +
  geom_errorbar(orientation='x', width=.01) +
  geom_line() +
  coord_fixed(xlim=0:1, ylim=0:1, expand=FALSE) +
  xlab('P(Type 2 False Alarm)') + ylab('P(Type 2 Hit)') +
  theme_bw(18)


################################################################################
#       Condition-level regression (variable confidence threhsolds)
################################################################################
d <- sim_metad_condition(N_trials=100000,
                         d_prime=c(1, 2), c=c(-1, 1),
                         log_M=c(-3/4, -1/3),
                         c2_0=list(c(.5, 1, 1.5),
                                   c(.25, .5, .75)),
                         c2_1=list(c(1/3, 2/3, 1),
                                   c(1/3, 2/3, 1)))

m <- fit_metad(bf(N ~ 0 + Intercept + condition,
                  dprime + c + metac2zero1diff + metac2zero2diff + metac2zero3diff +
                    metac2one1diff + metac2one2diff + metac2one3diff ~
                      0 + Intercept + condition),
               data=d, cores=4, backend='cmdstanr',
               prior=prior(normal(0, 1)) +
                 prior(normal(0, 1), dpar=dprime) +
                 prior(normal(0, 1), dpar=c) +
                 prior(normal(0, 1), dpar=metac2zero1diff) +
                 prior(normal(0, 1), dpar=metac2zero2diff) +
                 prior(normal(0, 1), dpar=metac2zero3diff) +
                 prior(normal(0, 1), dpar=metac2one1diff) +
                 prior(normal(0, 1), dpar=metac2one2diff) +
                 prior(normal(0, 1), dpar=metac2one3diff))
summary(m, prior=TRUE)

d |>
  distinct(condition) |>
  add_epred_draws(m) |>
  median_qi() |>
  group_by(condition) |>
  mutate(.true=response_probabilities(m$data$N[first(condition),]),
         .category=as.integer(.category),
         condition=factor(condition),
         joint_response=factor((.category-1) %% (n_distinct(d$confidence)*2) + 1),
         stimulus=factor(as.integer(.category > n_distinct(d$confidence)*2))) |>
  ggplot(aes(x=joint_response, fill=condition)) +
  geom_col(aes(y=.true)) +
  geom_pointrange(aes(y=.epred, ymin=.lower, ymax=.upper)) +
  facet_wrap(condition ~ stimulus, labeller=label_both) +
  theme_classic(18)

linpred_draws(m, newdata=tibble(condition=1:2), dpar=TRUE, transform=TRUE) |>
  pivot_longer(.linpred:metac2one3diff, names_to='.variable', values_to='.value') |>
  group_by(.variable, .add=TRUE) |>
  median_qi() |>
  arrange(.variable)

predicted_draws(m, newdata=select(m$data, condition, N)) |>
  group_by(condition, .row, .category) |>
  median_qi(.prediction) |>
  mutate(y=as.integer(t(m$data$N)),
         stimulus=rep(0:1, each=8, times=2),
         joint_response=rep(1:8, 4)) |>
  ggplot(aes(x=factor(joint_response))) +
  geom_bar(aes(y=y, fill=condition), position=position_dodge(1), stat='identity') +
  geom_pointrange(aes(y=.prediction, ymin=.lower, ymax=.upper, group=condition),
                  position=position_dodge(1)) +
  facet_grid( ~ stimulus) +
  theme_classic(18)


## calculate mean confidence per stimulus
d |>
  distinct(condition) |>
  add_mean_confidence_draws(m) |>
  median_qi(.epred) |>
  left_join(d |>
              group_by(condition, stimulus, response, confidence) |>
              summarize(across(theta_1:theta_2, first)) |>
              group_by(condition, stimulus) |>
              summarize(.true=sum(theta_1 * theta_2 * confidence)))


# psuedo type-1 ROC
d |>
  distinct(condition) |>
  mutate(condition=factor(condition)) |>
  add_roc1_draws(m, bounds=FALSE) |>
  group_by(condition, .add=TRUE) |>
  median_qi(p_fa, p_hit) |>
  ggplot(aes(x=p_fa, xmin=p_fa.lower, xmax=p_fa.upper,
             y=p_hit, ymin=p_hit.lower, ymax=p_hit.upper,
             group=condition, color=condition)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_errorbar(orientation='y', width=.01) +
  geom_errorbar(orientation='x', width=.01) +
  geom_line() +
  geom_point() +
  coord_fixed(xlim=0:1, ylim=0:1, expand=FALSE) +
  xlab('P(False Alarm)') + ylab('P(Hit)') +
  theme_bw()

# type 2 ROC
d |>
  distinct(condition) |>
  mutate(condition=factor(condition)) |>
  add_roc2_draws(m) |>
  median_qi(p_hit2, p_fa2) |>
  mutate(response=factor(response)) |>
  ggplot(aes(x=p_fa2, xmin=p_fa2.lower, xmax=p_fa2.upper,
             y=p_hit2, ymin=p_hit2.lower, ymax=p_hit2.upper,
             color=response, linetype=condition)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_errorbar(orientation='y', width=.01) +
  geom_errorbar(orientation='x', width=.01) +
  geom_line() +
  geom_point() +
  coord_fixed(xlim=0:1, ylim=0:1, expand=FALSE) +
  xlab('P(Type 2 False Alarm)') + ylab('P(Type 2 Hit)') +
  theme_bw(18)



## fit a model without fixed effects to test main effect
d.summary <- d |> metad_aggregate(condition)
m2 <- fit_metad(bf(N ~ 0 + Intercept + condition,
                   dprime + c ~ 0 + Intercept + condition,
                   metac2zero1diff + metac2zero2diff + metac2zero3diff +
                     metac2one1diff + metac2one2diff + metac2one3diff ~
                       0 + Intercept),
                data=d.summary, aggregate=FALSE,
                cores=4, backend='cmdstanr',
                prior=prior(normal(0, 1)) +
                  prior(normal(0, 1), dpar=dprime) +
                  prior(normal(0, 1), dpar=c) +
                  prior(normal(0, 1), dpar=metac2zero1diff) +
                  prior(normal(0, 1), dpar=metac2zero2diff) +
                  prior(normal(0, 1), dpar=metac2zero3diff) +
                  prior(normal(0, 1), dpar=metac2one1diff) +
                  prior(normal(0, 1), dpar=metac2one2diff) +
                  prior(normal(0, 1), dpar=metac2one3diff))
summary(m2, prior=TRUE)

d |>
  distinct(condition) |>
  add_epred_draws(m2) |>
  median_qi() |>
  group_by(condition) |>
  mutate(.true=response_probabilities(m$data$N[first(condition),]),
         .category=as.integer(.category),
         condition=factor(condition),
         joint_response=factor((.category-1) %% (n_distinct(d$confidence)*2) + 1),
         stimulus=factor(as.integer(.category > n_distinct(d$confidence)*2))) |>
  ggplot(aes(x=joint_response, fill=condition)) +
  geom_col(aes(y=.true)) +
  geom_pointrange(aes(y=.epred, ymin=.lower, ymax=.upper)) +
  facet_wrap(condition ~ stimulus, labeller=label_both) +
  theme_classic(18)

loo(m, m2)


################################################################################
#                          Participant-level regression
#                         (currently takes ~1.5mins to fit)
################################################################################
if (file.exists('../data/data_participant.csv')) {
  d <- read_csv('../data/data_participant.csv') |>
    mutate(meta_c2_0=map(meta_c2_0, ~ as.numeric(unlist(str_split(., '\\|')))),
           meta_c2_1=map(meta_c2_1, ~ as.numeric(unlist(str_split(., '\\|')))))
} else {
  d <- sim_metad_participant(N_participants=100, N_trials=500,
                             mu_d_prime=1, sd_d_prime=.5,
                             mu_c=0, sd_c=.5,
                             mu_log_M=0, sd_log_M=.5,
                             mu_z_c2=rep(-1, 3), sd_z_c2=rep(.1, 3),
                             r_z_c2=corr_matrix(.5, nrow=3))
  ## collapse list_cols to strings
  d |>
    mutate(meta_c2_0=map_chr(meta_c2_0, str_flatten, collapse='|'),
           meta_c2_1=map_chr(meta_c2_1, str_flatten, collapse='|')) |>
    write_csv('../data/data_participant.csv')
}


m <- fit_metad(bf(N ~ 0 + Intercept + (1 | participant),
                  dprime + c ~ 0 + Intercept + (1 | participant),
                  metac2zero1diff + metac2zero2diff + metac2zero3diff ~
                    0 + Intercept + (1 |p0| participant),
                  metac2one1diff + metac2one2diff + metac2one3diff
                  ~ 0 + Intercept + (1 |p1| participant)),
               data=d, init=0, cores=4, backend='cmdstanr',
               file='brms_participant', file_refit='on_change',
               prior=prior(normal(0, 1)) +
                 prior(normal(0, 1), dpar=dprime) +
                 prior(normal(0, 1), dpar=c) +
                 prior(normal(0, 1), dpar=metac2zero1diff) +
                 prior(normal(0, 1), dpar=metac2zero2diff) +
                 prior(normal(0, 1), dpar=metac2zero3diff) +
                 prior(normal(0, 1), dpar=metac2one1diff) +
                 prior(normal(0, 1), dpar=metac2one2diff) +
                 prior(normal(0, 1), dpar=metac2one3diff))
summary(m, prior=TRUE)


d |>
  distinct(participant) |>
  add_linpred_draws(m, dpar=TRUE, transform=TRUE) |>
  rename(M=mu, d_prime=dprime) |>
  mutate(meta_c2_0_1=M*c - metac2zero1diff,
         meta_c2_0_2=meta_c2_0_1 - metac2zero2diff,
         meta_c2_0_3=meta_c2_0_2 - metac2zero3diff,
         meta_c2_1_1=M*c + metac2one1diff,
         meta_c2_1_2=meta_c2_1_1 + metac2one2diff,
         meta_c2_1_3=meta_c2_1_2 + metac2one3diff) |>
  select(-ends_with('diff'), -.linpred) |>
  pivot_longer(M:meta_c2_1_3, names_to='.variable', values_to='.value') |>
  group_by(participant, .variable) |>
  median_qi() |>
  left_join(
    d |>
      group_by(participant) |>
      reframe(across(d_prime:meta_c2_1, first)) |>
      group_by(participant) |>
      mutate(conf=row_number()) |>
      pivot_wider(names_from=conf, values_from=meta_c2_0:meta_c2_1) |>
      pivot_longer(d_prime:meta_c2_1_3, names_to='.variable', values_to='.true')) |>
  ggplot(aes(y=.value, x=.true)) +
  geom_abline(linetype='dashed', slope=1, intercept=0) +
  geom_pointrange(aes(ymin=.lower, ymax=.upper)) +
  facet_wrap(~ .variable, scales='free') +
  theme_bw()


## calculate mean confidence per stimulus
d |>
  distinct(participant) |>
  add_mean_confidence_draws(m) |>
  median_qi(.epred) |>
  left_join(d |>
              group_by(participant, stimulus, response, confidence) |>
              summarize(across(theta_1:theta_2, first)) |>
              group_by(participant, stimulus) |>
              summarize(.true=sum(theta_1 * theta_2 * confidence))) |>
  mutate(stimulus=factor(stimulus)) |>
  ggplot(aes(x=.true, y=.epred, ymin=.lower, ymax=.upper, color=stimulus)) +
  geom_abline(linetype='dashed', slope=1, intercept=0) +
  geom_pointrange() +
  facet_wrap(stimulus ~ ., labeller=label_both) +
  coord_fixed(xlim=c(1, 4), ylim=c(1, 4)) +
  theme_bw()


# psuedo type-1 ROC
tibble(.row=1) |>
  add_roc1_draws(m, re_formula=NA) |>
  median_qi(p_fa, p_hit) |>
  ggplot(aes(x=p_fa, xmin=p_fa.lower, xmax=p_fa.upper,
             y=p_hit, ymin=p_hit.lower, ymax=p_hit.upper)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_line(aes(x=p_fa, y=p_hit, group=participant),
            alpha=.1, inherit.aes=FALSE,
            data=d |>
              distinct(participant) |>
              mutate(participant=factor(participant)) |>
              add_roc1_draws(m) |>
              group_by(participant, .add=TRUE) |>
              median_qi(p_fa, p_hit)) +
  geom_errorbar(orientation='y', width=.01) +
  geom_errorbar(orientation='x', width=.01) +
  geom_line() +
  coord_fixed(xlim=0:1, ylim=0:1, expand=FALSE) +
  xlab('P(False Alarm)') + ylab('P(Hit)') +
  theme_bw()


# type 2 ROC
tibble(.row=1) |>
  add_roc2_draws(m, re_formula=NA) |>
  median_qi(p_hit2, p_fa2) |>
  mutate(response=factor(response)) |>
  ggplot(aes(x=p_fa2, xmin=p_fa2.lower, xmax=p_fa2.upper,
             y=p_hit2, ymin=p_hit2.lower, ymax=p_hit2.upper)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_line(aes(x=p_fa2, y=p_hit2, group=participant),
            alpha=.1, inherit.aes=FALSE,
            data=d |>
              distinct(participant) |>
              mutate(participant=factor(participant)) |>
              add_roc2_draws(m) |>
              group_by(participant, .add=TRUE) |>
              median_qi(p_fa2, p_hit2)) +
  geom_errorbar(orientation='y', width=.01) +
  geom_errorbar(orientation='x', width=.01) +
  geom_line() +
  geom_point() +
  coord_fixed(xlim=0:1, ylim=0:1, expand=FALSE) +
  xlab('P(Type 2 False Alarm)') + ylab('P(Type 2 Hit)') +
  facet_wrap(~ response, labeller=label_both) +
  theme_bw(18)



################################################################################
#                     Participant-level regression by condition
#                         (currently takes ~3.5mins to fit)
################################################################################
if (file.exists('../data/data_participant_condition.csv')) {
  d <- read_csv('../data/data_participant_condition.csv') |>
    mutate(meta_c2_0=map(meta_c2_0, ~ as.numeric(unlist(str_split(., '\\|')))),
           meta_c2_1=map(meta_c2_1, ~ as.numeric(unlist(str_split(., '\\|')))))
} else {
  d <- sim_metad_participant_condition(
    N_participants=100, N_trials=500,
    mu_d_prime=rep(1, 2), sd_d_prime=rep(.5, 2), r_d_prime=corr_matrix(.5),
    mu_c=rep(0, 2), sd_c=rep(.5, 2), r_c=corr_matrix(.25),
    mu_log_M=rep(0, 2), sd_log_M=rep(.75, 2), r_log_M=corr_matrix(.5),
    mu_z_c2=cbind(c(-1,-1,-1), c(-1.5,-1.5,-1.5)),
    r_z_c2_condition=corr_matrix(.9))
  ## collapse list_cols to strings
  d |>
    mutate(meta_c2_0=map_chr(meta_c2_0, str_flatten, collapse='|'),
           meta_c2_1=map_chr(meta_c2_1, str_flatten, collapse='|')) |>
    write_csv('../data/data_participant_condition.csv')
}


m <- fit_metad(bf(N ~ 0 + Intercept + condition + (condition | participant),
                  dprime + c ~ 0 + Intercept + condition + (condition | participant),
                  metac2zero1diff + metac2zero2diff + metac2zero3diff ~
                    0 + Intercept + condition + (condition |p0| participant),
                  metac2one1diff + metac2one2diff + metac2one3diff ~
                    0 + Intercept + condition + (condition |p1| participant)),
               data=d, init=0, cores=4, backend='cmdstanr',
               file='brms_participant_condition', file_refit='on_change',
               prior=prior(normal(0, 1)) +
                 prior(normal(0, 1), dpar=dprime) +
                 prior(normal(0, 1), dpar=c) +
                 prior(normal(0, 1), dpar=metac2zero1diff) +
                 prior(normal(0, 1), dpar=metac2zero2diff) +
                 prior(normal(0, 1), dpar=metac2zero3diff) +
                 prior(normal(0, 1), dpar=metac2one1diff) +
                 prior(normal(0, 1), dpar=metac2one2diff) +
                 prior(normal(0, 1), dpar=metac2one3diff))
summary(m, prior=TRUE)


d |>
  distinct(condition) |>
  add_linpred_draws(m, dpar=TRUE, re_formula=NA) |>
  select(-.linpred) |>
  rename(M=mu) |>
  pivot_longer(M:metac2one3diff, names_to='.variable', values_to='.value') |>
  group_by(.variable, condition) |>
  median_qi()


draws.participant <- d |>
  distinct(participant, condition) |>
  add_linpred_draws(m, dpar=TRUE, transform=TRUE) |>
  rename(M=mu, d_prime=dprime) |>
  mutate(meta_c2_0_1=M*c - metac2zero1diff,
         meta_c2_0_2=meta_c2_0_1 - metac2zero2diff,
         meta_c2_0_3=meta_c2_0_2 - metac2zero3diff,
         meta_c2_1_1=M*c + metac2one1diff,
         meta_c2_1_2=meta_c2_1_1 + metac2one2diff,
         meta_c2_1_3=meta_c2_1_2 + metac2one3diff) |>
  select(-ends_with('diff'), -.linpred) |>
  pivot_longer(M:meta_c2_1_3, names_to='.variable', values_to='.value') |>
  group_by(participant, condition, .variable) |>
  median_qi() |>
  left_join(
    d |>
      group_by(participant, condition) |>
      reframe(across(d_prime:meta_c2_1, first)) |>
      group_by(participant, condition) |>
      mutate(conf=row_number()) |>
      pivot_wider(names_from=conf, values_from=meta_c2_0:meta_c2_1) |>
      pivot_longer(d_prime:meta_c2_1_3, names_to='.variable', values_to='.true'))

ggplot(draws.participant, aes(y=.value, x=.true)) +
  geom_abline(linetype='dashed', slope=1, intercept=0) +
  geom_pointrange(aes(ymin=.lower, ymax=.upper, color=factor(condition))) +
  facet_wrap(~ .variable, scales='free') +
  theme_bw()



## calculate mean confidence per stimulus
draws.mean <- d |>
  distinct(participant, condition) |>
  add_epred_draws(m, ndraws=100) |>
  group_by(participant, condition, .draw) |>
  mutate(stimulus=rep(0:1, each=n_distinct(d$confidence)*2)) |>
  group_by(participant, condition, stimulus, .draw) |>
  summarize(.epred=list(.epred)) |>
  mutate(.mean=map_dbl(.epred, `%*%`,
                       c(n_distinct(d$confidence):1,
                         1:(n_distinct(d$confidence))))) |>
  median_qi(.mean) |>
  left_join(d |>
              group_by(participant, condition, stimulus, response, confidence) |>
              summarize(across(theta_1:theta_2, first)) |>
              group_by(participant, condition, stimulus) |>
              summarize(.true=sum(theta_1 * theta_2 * confidence)))

draws.mean |>
  mutate(stimulus=factor(stimulus)) |>
  ggplot(aes(x=.true, y=.mean, ymin=.lower, ymax=.upper, color=stimulus)) +
  geom_abline(linetype='dashed', slope=1, intercept=0) +
  geom_pointrange() +
  facet_grid(stimulus ~ condition, labeller=label_both) +
  coord_fixed(xlim=c(1, 4), ylim=c(1, 4)) +
  theme_bw()

