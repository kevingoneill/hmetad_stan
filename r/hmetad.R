library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(mvtnorm)

source('metad_utils.R')

## simulate an SDT agent for N trials
d <- sim_sdt_participant(
    N_participants=100, N_trials=100,
    mu_d_prime=1, sd_d_prime=.5,
    mu_c=0, sd_c=.5,
    mu_log_M=0, sd_log_M=.5)


## format data for stan
data.simulated.stan <- list(N=nrow(d),
                            K=n_distinct(d$confidence),
                            P=n_distinct(d$participant),
                            stimulus=d$stimulus,
                            response=d$response,
                            confidence=d$confidence,
                            participant=d$participant,
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
m <- cmdstan_model('../stan/hmetad.stan')

prior <- m$sample(c(data.simulated.stan, prior_only=TRUE), chains=4, parallel_chains=4, init=0)

fit <- m$sample(c(data.simulated.stan, prior_only=FALSE), chains=4, parallel_chains=4, init=0)

fit$summary(c('mu_d_prime', 'mu_c', 'mu_log_M',
              'sigma_d_prime', 'sigma_c', 'sigma_log_M',
              'mu_z_meta_c2_0', 'mu_z_meta_c2_1',
              'sigma_meta_c2'##, 'Omega_meta_c2[2,1]'
              )) %>%
    print(., n=nrow(.))



################################################################################
##                      Plot participant-level recovery
################################################################################
d.average <- d %>%
    group_by(participant) %>%
    summarize(d_prime=mean(d_prime),
              c=mean(c),
              log_M=mean(log(M))) %>%
    pivot_longer(d_prime:log_M, names_to='.variable', values_to='.true_value') %>%
    left_join(fit %>%
              gather_draws(d_prime[participant],
                           c[participant],
                           M[participant]) %>%
              mutate(.value=ifelse(.variable=='M', log(.value), .value),
                     .variable=ifelse(.variable=='M', 'log_M', .variable)) %>%
              median_qi())

d.average %>%
    filter(.variable=='d_prime') %>%
    ggplot(aes(x=.true_value, y=.value, ymin=.lower, ymax=.upper)) +
    geom_abline(intercept=0, slope=1, linetype='dashed') +
    geom_pointrange(size=.25) +
    xlab('True d\'') + ylab('Estimated d\'') +
    coord_fixed(xlim=round(range(d$d_prime)) + c(-.5, .5),
                ylim=round(range(d$d_prime)) + c(-.5, .5)) +
    theme_classic(18) +
    theme(panel.grid.major=element_line(linewidth=.5, color='grey80'),
          panel.border=element_rect(linewidth=1.5),
          axis.line=element_blank())
ggsave('../plots/hmetad/recovery_d_prime.png', width=8, height=5)

d.average %>%
    filter(.variable=='c') %>%
    ggplot(aes(x=.true_value, y=.value, ymin=.lower, ymax=.upper)) +
    geom_abline(intercept=0, slope=1, linetype='dashed') +
    geom_pointrange(size=.25) +
    xlab('True c') + ylab('Estimated c') +
    coord_fixed(xlim=round(range(d$c)) + c(-.5, .5),
                ylim=round(range(d$c)) + c(-.5, .5)) +
    theme_classic(18) +
    theme(panel.grid.major=element_line(linewidth=.5, color='grey80'),
          panel.border=element_rect(linewidth=1.5),
          axis.line=element_blank())
ggsave('../plots/hmetad/recovery_c.png', width=8, height=5)

d.average %>%
    filter(.variable=='log_M') %>%
    ggplot(aes(x=.true_value, y=.value, ymin=.lower, ymax=.upper)) +
    geom_abline(intercept=0, slope=1, linetype='dashed') +
    geom_pointrange(size=.25) +
    xlab('True log(M)') + ylab('Estimated log(M)') +
    coord_fixed(xlim=round(range(log(d$M))) + c(-1, 1),
                ylim=round(range(log(d$M))) + c(-1, 1)) +
    theme_classic(18) +
    theme(panel.grid.major=element_line(linewidth=.5, color='grey80'),
          panel.border=element_rect(linewidth=1.5),
          axis.line=element_blank())
ggsave('../plots/hmetad/recovery_log_M.png', width=8, height=5)

################################################################################
##                      Plot estimated correlations
################################################################################
fit %>%
    spread_draws(Omega_meta_c2[i,j]) %>%
    filter(i==1, j==2) %>%
    ggplot(aes(x=Omega_meta_c2)) +
    stat_histinterval(.width=.95) +
    geom_vline(aes(xintercept=Omega_meta_c2),
               data=tibble(Omega_meta_c2=0)) +
    theme_classic(18) +
    theme(panel.grid.major.x=element_line(linewidth=.5, color='grey80'),
          panel.border=element_rect(linewidth=1.5),
          axis.line=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank())
ggsave('../plots/hmetad/recovery_corr.png', width=8, height=4)
