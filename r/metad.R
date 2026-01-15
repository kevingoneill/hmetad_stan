library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(mvtnorm)

source('metad_utils.R')

## simulate an SDT agent for N trials
d <- sim_sdt(N_trials=100000, d_prime=1, c=0, log_M=0)


## format data for stan
data.simulated.stan <- list(N=nrow(d),
                            K=n_distinct(d$confidence),
                            stimulus=d$stimulus,
                            response=d$response,
                            confidence=d$confidence,
                            prior_sd_d_prime=1,
                            prior_sd_c=1, 
                            prior_sd_log_M=1,
                            prior_mean_meta_c2=-1,
                            prior_sd_meta_c2=.5)

## fit stan model
m <- cmdstan_model('../stan/metad.stan')

prior <- m$sample(c(data.simulated.stan, prior_only=TRUE), chains=4, parallel_chains=4, init=0)

fit <- m$sample(c(data.simulated.stan, prior_only=FALSE),
                chains=4, parallel_chains=4, init=0)

fit$summary(c('d_prime', 'c', 'M', 'meta_c2_0', 'meta_c2_1')) %>%
    print(., n=nrow(.))



################################################################################
##                      Plot group-level recovery
################################################################################
draws <- d %>%
    summarize(d_prime=first(d_prime),
              c=first(c),
              log_M=first(log(M)),
              meta_c2_0_1=first(meta_c2_0)[1],
              meta_c2_1_1=first(meta_c2_1)[1],
              meta_c2_0_2=first(meta_c2_0)[2],
              meta_c2_1_2=first(meta_c2_1)[2],
              meta_c2_0_3=first(meta_c2_0)[3],
              meta_c2_1_3=first(meta_c2_1)[3]) %>%
  pivot_longer(d_prime:meta_c2_1_3, names_to='.variable',
               values_to='.true_value') %>%
    right_join(fit %>%
              gather_draws(d_prime, c, M, meta_c2_0[k], meta_c2_1[k]) %>%
              mutate(.value=ifelse(.variable=='M', log(.value), .value),
                     .variable=ifelse(.variable=='M', 'log_M', .variable)) %>%
              mutate(.variable=ifelse(!is.na(k), paste0(.variable, '_', k), .variable)) %>%
              group_by(.variable) %>%
              select(-k))


ggplot(draws, aes(x=.value)) +
    stat_halfeye(.width=.95) +
    geom_vline(aes(xintercept=.true_value)) +
    facet_wrap(~ .variable, scales='free_x') +
    theme_classic(18) +
    theme(panel.grid.major.x=element_line(linewidth=.5, color='grey80'),
          panel.border=element_rect(linewidth=1.5, fill=NA),
          axis.line=element_blank(),
          axis.title=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank())
ggsave('../plots/metad/recovery.png', width=8, height=8)


# pseudo-type 1 ROC
fit |>
  spread_draws(ROC_1[k, stimulus]) |>
  mutate(stimulus=stimulus-1) |>
  pivot_wider(names_from=stimulus, values_from=ROC_1, names_prefix='p_') |>
  median_qi(p_0, p_1) |>
  ggplot(aes(x=p_0, xmin=p_0.lower, xmax=p_0.upper,
             y=p_1, ymin=p_1.lower, ymax=p_1.upper)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_errorbar(orientation='y', width=.01) +
  geom_errorbar(orientation='x', width=.01) +
  geom_line() +
  geom_point() +
  coord_fixed(xlim=0:1, ylim=0:1, expand=FALSE) +
  xlab('P(False Alarm)') + ylab('P(Hit)') +
  theme_bw()

# pseudo-type 2 ROC
fit |>
  spread_draws(ROC_2[response, k, accuracy]) |>
  mutate(response=response-1,
         accuracy=accuracy-1) |>
  pivot_wider(names_from=accuracy, values_from=ROC_2,
              names_prefix='p_') |>
  median_qi(p_0, p_1) |>
  ggplot(aes(x=p_0, xmin=p_0.lower, xmax=p_0.upper,
             y=p_1, ymin=p_1.lower, ymax=p_1.upper)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_errorbar(orientation='y', width=.01) +
  geom_errorbar(orientation='x', width=.01) +
  geom_line() +
  geom_point() +
  coord_fixed(xlim=0:1, ylim=0:1, expand=FALSE) +
  xlab('P(False Alarm 2)') + ylab('P(Hit 2)') +
  theme_bw()

