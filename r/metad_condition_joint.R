library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(mvtnorm)

source('metad_utils.R')

## simulate an SDT agent for N trials
d <- sim_sdt_condition(N_trials=10000,
                       d_prime=c(1, 2), c=c(-1, 1),
                       log_M=c(-3/4, -1/3))


## format data for stan
data.simulated.stan <- list(N=nrow(d),
                            K=n_distinct(d$confidence),
                            B=n_distinct(d$condition),
                            stimulus=d$stimulus,
                            response=d$response,
                            confidence=d$confidence,
                            condition=d$condition,
                            prior_sd_d_prime=1,
                            prior_sd_c=1, 
                            prior_sd_log_M=1,
                            prior_mean_meta_c2=-1,
                            prior_sd_meta_c2=.5)

## fit stan model
m <- cmdstan_model('../stan/metad_condition_joint.stan')

prior <- m$sample(c(data.simulated.stan, prior_only=TRUE), chains=4, parallel_chains=4, init=0)

fit <- m$sample(c(data.simulated.stan, prior_only=FALSE), chains=4, parallel_chains=4, init=0)

fit$summary(c('d_prime', 'c', 'M', 'meta_c2_0', 'meta_c2_1')) %>%
    print(., n=nrow(.))



################################################################################
##                      Plot participant-level recovery
################################################################################
draws <- d %>%
    group_by(condition) %>%
    summarize(d_prime=first(d_prime),
              c=first(c),
              log_M=first(log(M)),
              meta_c2_0_1=first(meta_c2_0)[1],
              meta_c2_1_1=first(meta_c2_1)[1],
              meta_c2_0_2=first(meta_c2_0)[2],
              meta_c2_1_2=first(meta_c2_1)[2],
              meta_c2_0_3=first(meta_c2_0)[3],
              meta_c2_1_3=first(meta_c2_1)[3]) %>%
    pivot_longer(d_prime:meta_c2_1_3, names_to='.variable', values_to='.true_value') %>%
    right_join(fit %>%
               gather_draws(d_prime[condition], c[condition], M[condition],
                            meta_c2_0[k,condition], meta_c2_1[k,condition]) %>%
               mutate(.value=ifelse(.variable=='M', log(.value), .value),
                      .variable=ifelse(.variable=='M', 'log_M', .variable)) %>%
               mutate(.variable=ifelse(!is.na(k), paste0(.variable, '_', k), .variable)) %>%
               group_by(.variable) %>%
               select(-k))


ggplot(draws, aes(x=.value)) +
    geom_vline(aes(xintercept=.true_value, color=factor(condition))) +
    stat_halfeye(aes(fill=factor(condition), group=condition), .width=.95, slab_alpha=.5) +
    facet_wrap( ~ .variable, scales='free_x') +
    scale_fill_discrete(name='Condition') +
    scale_color_discrete(name='Condition') +
    theme_classic(18) +
    theme(panel.grid.major.x=element_line(linewidth=.5, color='grey80'),
          panel.border=element_rect(linewidth=1.5, fill=NA),
          axis.line=element_blank(),
          axis.title=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank())
ggsave('../plots/metad_condition_joint/recovery.png', width=8, height=6)
