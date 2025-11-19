library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(mvtnorm)
library(bayesplot)

source('metad_utils.R')

## simulate an SDT agent for N trials
d <- sim_sdt(N_trials=10000, d_prime=1, c=0, log_M=0)


## format data for stan
data.simulated.stan <- list(N=nrow(d),
                            K=n_distinct(d$confidence),
                            stimulus=d$stimulus,
                            response=d$response,
                            confidence=d$confidence,
                            prior_sd_d_prime=1,
                            prior_sd_c=1, 
                            prior_sd_log_M=1,
                            prior_alpha_meta_c2=2)

m <- cmdstan_model('../stan/metad_dirichlet.stan')


prior <- m$sample(c(data.simulated.stan, prior_only=TRUE),
                  iter_warmup=3000,
                  chains=4, parallel_chains=4, adapt_delta=.99)

fit <- m$sample(c(data.simulated.stan, prior_only=FALSE),
                  chains=4, parallel_chains=4)


mcmc_pairs(prior$draws('meta_c2_1'))
mcmc_pairs(fit$draws('meta_c2_1'))

prior %>%
    spread_draws(theta_2[response,correct,confidence]) %>%
    mutate(response=factor(response-1),
           correct=factor(correct-1)) %>%
    ggplot(aes(x=confidence, y=theta_2,
               group=correct,
               fill=correct)) +
    stat_ccdfinterval(aes(slab_alpha = after_stat(f)),
                      thickness=1, .width=.95, position='dodge',
                      show.legend=FALSE) +
    facet_grid(response ~ ., labeller=label_both) +
    theme_classic(18)

fit %>%
    spread_draws(theta_2[response,correct,confidence]) %>%
    mutate(response=factor(response-1),
           correct=factor(correct-1)) %>%
    ggplot(aes(x=confidence, y=theta_2,
               group=correct,
               fill=correct)) +
    stat_ccdfinterval(aes(slab_alpha = after_stat(f)),
                      thickness=1, .width=.95, position='dodge',
                      show.legend=FALSE) +
    geom_errorbar(aes(y=p, ymin=p, ymax=p), position=position_dodge(1),
                 data=d %>% count(response, correct, confidence) %>%
                     group_by(response, correct) %>%
                     mutate(p=n/sum(n),
                            response=factor(response),
                            correct=factor(correct))) +
    facet_grid(response ~ ., labeller=label_both) +
    theme_classic(18)

fit %>%
    spread_draws(meta_c2_0[k], meta_c2_1[k]) %>%
    group_by(.variable, .draw) %>%
    arrange(.draw) %>%
    left_join(spread_draws(fit, meta_c)) %>%
    mutate()


## inspect prior over z_meta_c2 (on log difference scale)
draws.zc2 <- prior %>%
    spread_draws(meta_c, meta_c2_0[k], meta_c2_1[k]) %>%
    group_by(.draw) %>%
    arrange(.draw) %>%
    mutate(d_meta_c2_0=-diff(c(meta_c[1], meta_c2_0)),
           d_meta_c2_1=diff(c(meta_c[1], meta_c2_1)),
           z_meta_c2_0=log(d_meta_c2_0),
           z_meta_c2_1=log(d_meta_c2_1)) %>%
    select(.draw, k, z_meta_c2_0, z_meta_c2_1) %>%
    pivot_wider(names_from=k, values_from=z_meta_c2_0:z_meta_c2_1)

draws.zc2 %>%
    ungroup %>%
    select(starts_with('z_meta_c2_0')) %>%
    as.matrix() %>%
    mcmc_pairs()

draws.zc2 %>%
    ungroup %>%
    select(starts_with('z_meta_c2_1')) %>%
    as.matrix() %>%
    mcmc_pairs()
