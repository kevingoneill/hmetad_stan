library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(mvtnorm)

source('metad_utils.R')

## simulate an SDT agent for N trials
d <- sim_metad_participant_condition(
  N_participants=100, N_trials=100,
  mu_d_prime=rep(1, 2), sd_d_prime=rep(.5, 2), r_d_prime=corr_matrix(.5),
  mu_c=rep(0, 2), sd_c=rep(.5, 2), r_c=corr_matrix(.25),
  mu_log_M=rep(0, 2), sd_log_M=rep(.5, 2), r_log_M=corr_matrix(.5))


## format data for stan
data.simulated.stan <- list(N=nrow(d),
                            K=n_distinct(d$confidence),
                            P=n_distinct(d$participant),
                            W=n_distinct(d$condition),
                            stimulus=d$stimulus,
                            response=d$response,
                            confidence=d$confidence,
                            participant=d$participant,
                            condition=d$condition,
                            prior_sd_mu_d_prime=1,
                            prior_sd_mu_c=1, 
                            prior_sd_mu_log_M=1,
                            prior_mean_mu_meta_c2=-1,
                            prior_sd_mu_meta_c2=.5,
                            prior_sd_sigma_d_prime=2,
                            prior_sd_sigma_c=2,
                            prior_sd_sigma_log_M=2,
                            prior_sd_sigma_meta_c2=2,
                            prior_eta=1)

## fit stan model
m <- cmdstan_model('../stan/hmetad_condition.stan')

prior <- m$sample(c(data.simulated.stan, prior_only=TRUE), chains=4, parallel_chains=4, init=0)

fit <- m$sample(c(data.simulated.stan, prior_only=FALSE), chains=4, parallel_chains=4, init=0)

fit$summary(c('mu_d_prime', 'mu_c', 'mu_log_M',
              'sigma_d_prime', 'sigma_c', 'sigma_log_M',
              'mu_z_meta_c2_0', 'mu_z_meta_c2_1',
              'sigma_meta_c2_condition', 'sigma_meta_c2_confidence',
              'Omega_d_prime[2,1]', 'Omega_c[2,1]', 'Omega_log_M[2,1]',
              'Omega_meta_c2_condition[2,1]', 'Omega_meta_c2_confidence[2,1]')) %>%
  print(., n=nrow(.))



################################################################################
##                      Plot participant-level recovery
################################################################################
d.average <- d %>%
  group_by(participant, condition) %>%
  filter(trial==1) %>%
  mutate(log_M=mean(log(M)),
         meta_c2_0_1=map_dbl(meta_c2_0, ~ .[1]),
         meta_c2_0_2=map_dbl(meta_c2_0, ~ .[2]),
         meta_c2_0_3=map_dbl(meta_c2_0, ~ .[3]),
         meta_c2_1_1=map_dbl(meta_c2_1, ~ .[1]),
         meta_c2_1_2=map_dbl(meta_c2_1, ~ .[2]),
         meta_c2_1_3=map_dbl(meta_c2_1, ~ .[3])) %>%
  select(participant, condition, d_prime, c, log_M, meta_c2_0_1:meta_c2_1_3) %>%
  pivot_longer(d_prime:meta_c2_1_3, names_to='.variable', values_to='.true_value') %>%
  left_join(fit %>%
              gather_draws(d_prime[condition, participant],
                           c[condition, participant],
                           M[condition, participant],
                           meta_c2_0[condition, participant, k],
                           meta_c2_1[condition, participant, k]) %>%
              mutate(.value=ifelse(.variable=='M', log(.value), .value),
                     .variable=ifelse(.variable=='M', 'log_M', .variable),
                     .variable=ifelse(is.na(k), .variable, paste0(.variable, '_', k))) %>%
              median_qi())

d.average %>%
  filter(.variable=='d_prime') %>%
  ggplot(aes(x=.true_value, y=.value, ymin=.lower, ymax=.upper)) +
  geom_abline(intercept=0, slope=1, linetype='dashed') +
  geom_pointrange(size=.25) +
  facet_wrap(~ condition, labeller=label_both) +
  xlab('True d\'') + ylab('Estimated d\'') +
  coord_fixed(xlim=round(range(d$d_prime)) + c(-.5, .5),
              ylim=round(range(d$d_prime)) + c(-.5, .5)) +
  theme_classic(18) +
  theme(panel.grid.major=element_line(linewidth=.5, color='grey80'),
        panel.border=element_rect(linewidth=1.5, fill=NA),
        axis.line=element_blank())
ggsave('../plots/hmetad_condition/recovery_d_prime.png', width=8, height=5)

d.average %>%
  filter(.variable=='c') %>%
  ggplot(aes(x=.true_value, y=.value, ymin=.lower, ymax=.upper)) +
  geom_abline(intercept=0, slope=1, linetype='dashed') +
  geom_pointrange(size=.25) +
  facet_wrap(~ condition, labeller=label_both) +
  xlab('True c') + ylab('Estimated c') +
  coord_fixed(xlim=round(range(d$c)) + c(-.5, .5),
              ylim=round(range(d$c)) + c(-.5, .5)) +
  theme_classic(18) +
  theme(panel.grid.major=element_line(linewidth=.5, color='grey80'),
        panel.border=element_rect(linewidth=1.5, fill=NA),
        axis.line=element_blank())
ggsave('../plots/hmetad_condition/recovery_c.png', width=8, height=5)

d.average %>%
  filter(.variable=='log_M') %>%
  ggplot(aes(x=.true_value, y=.value, ymin=.lower, ymax=.upper)) +
  geom_abline(intercept=0, slope=1, linetype='dashed') +
  geom_pointrange(size=.25) +
  facet_wrap(~ condition, labeller=label_both) +
  xlab('True log(M)') + ylab('Estimated log(M)') +
  coord_fixed(xlim=round(range(log(d$M))) + c(-1, 1),
              ylim=round(range(log(d$M))) + c(-1, 1)) +
  theme_classic(18) +
  theme(panel.grid.major=element_line(linewidth=.5, color='grey80'),
        panel.border=element_rect(linewidth=1.5, fill=NA),
        axis.line=element_blank())
ggsave('../plots/hmetad_condition/recovery_log_M.png', width=8, height=5)

d.average %>%
  filter(str_starts(.variable, 'meta_c2')) %>%
  separate_wider_regex(.variable, patterns=c(.variable='meta_c2', '_',
                                             response='.', '_', confidence='.')) %>%
  ggplot(aes(x=.true_value, y=.value, ymin=.lower, ymax=.upper)) +
  geom_abline(intercept=0, slope=1, linetype='dashed') +
  geom_pointrange(size=.25) +
  facet_grid(response+condition ~ confidence, labeller=label_both) +
  xlab('True meta-c2') + ylab('Estimated meta-c2') +
  ##coord_fixed() +
  theme_classic(18) +
  theme(panel.grid.major=element_line(linewidth=.5, color='grey80'),
        panel.border=element_rect(linewidth=1.5, fill=NA),
        axis.line=element_blank())

################################################################################
##                      Plot estimated correlations
################################################################################
fit %>%
  gather_draws(Omega_d_prime[i,j], Omega_c[i,j], Omega_log_M[i,j],
               Omega_meta_c2_condition[i,j], Omega_meta_c2_confidence[i,j]) %>%
  filter(i==1, j==2) %>%
  mutate(.variable=str_replace_all(.variable, 'Omega', 'r')) %>%
  ggplot(aes(x=.value)) +
  stat_histinterval(.width=.95) +
  geom_vline(aes(xintercept=.value),
             data=tibble(.variable=c('r_d_prime', 'r_c', 'r_log_M',
                                     'r_meta_c2_condition',
                                     'r_meta_c2_confidence'),
                         .value=c(.5, .25, .5, 0, 0))) +
  facet_wrap(~ .variable) +
  theme_classic(18) +
  theme(panel.grid.major.x=element_line(linewidth=.5, color='grey80'),
        panel.border=element_rect(linewidth=1.5, fill=NA),
        axis.line=element_blank(),
        axis.title=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())
ggsave('../plots/hmetad_condition/recovery_corr.png', width=8, height=4)
