
library("tidyverse")
library("brms")
library("tidybayes")
library("cmdstanr")


source("sim_fun.R")

#### compile model:

### data needed for compliation

comp_dat <- gen_one_data(N_cond = 500, g = 0.85, e = 0.2, mu = 0.4, 
             diff = 0.1, phi = exp(1.11))
str(comp_dat)

comp_dat %>% 
  ggplot(aes(x = response)) +
  geom_histogram(binwidth = 0.02, boundary = 0) +
  facet_wrap(vars(condition))

head(comp_dat)

### model
zoib2 <- custom_family(
  "zoib2", dpars = c("mu", "phi", "zipp", "coi"),
  links = c("logit", "log", "logit", "logit"), lb = c(NA, 0, NA, NA),
  type = "real"
)

stan_funs <- "
/* zero-one-inflated beta log-PDF of a single response 
   * Args: 
   *   y: response value 
   *   mu: mean parameter of the beta part
   *   phi: precision parameter of the beta part
   *   zipp: zero-inflation probability parameter
   *   coi: conditional one-inflation probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zoib2_lpdf(real y, real mu, real phi,
                                    real zipp, real coi) {
     row_vector[2] shape = [mu * phi, (1 - mu) * phi]; 
     if (y == 0) { 
       return bernoulli_lpmf(0 | zipp); 
     } else if (y == 1) {
       return bernoulli_lpmf(1 | zipp) + bernoulli_lpmf(1 | coi);
     } else { 
       return bernoulli_lpmf(1 | zipp) + bernoulli_lpmf(0 | coi) + beta_lpdf(y | shape[1], shape[2]);
     } 
   }
"
stanvars <- stanvar(scode = stan_funs, block = "functions")

zoib_model <- bf(
  response ~ 0 + Intercept + condition,
  phi ~ 1,
  zipp ~ 0 + Intercept + condition,
  coi ~ 0 + Intercept + condition, 
  family = zoib2, center = FALSE, cmc = FALSE, nl = FALSE
)

# np <- prior_string("target += logistic_lpdf(Intercept_zipp | 0, 1);", 
#                    check = FALSE)

make_stancode(zoib_model, data = comp_dat, stanvars = stanvars, #prior = np, 
              save_model = "zoibr_model.stan")

compiled_model <- cmdstan_model("zoibr_model.stan")
compiled_model$print()

#### run one simulation:

d1 <- gen_one_data(N_cond = 500, g = 0.85, e = 0.2, mu = 0.4, 
                   diff = 0.1, phi = exp(1.11))
fit_one_dat(data = d1, compiled_model = compiled_model, 
            formula = zoib_model, stanvars = stanvars)



### multicore simulations: 

library("parallel")

## testing:
nsim <- 20
sim_df_05 <- make_sim_df(nsim = nsim, diff = 0.05)

run_one_sim(index = 1, sim_df = sim_df_05, 
            compiled_model = compiled_model, 
            formula = zoib_model, stanvars = stanvars)
sim_out_05 <- mclapply(seq_len(nsim), FUN = run_one_sim, sim_df = sim_df_05, 
                       compiled_model = compiled_model, 
                       formula = zoib_model, stanvars = stanvars, 
                       mc.cores = 16, mc.preschedule = FALSE, 
                       mc.allow.recursive = FALSE)
bind_rows(sim_out_05) %>% 
  group_by(parameter) %>% 
  summarise(mean(include_0))

mean(bind_rows(sim_out_05)$obs_diff)

#### main simulations:

NSIM <- 5000
differences <- c(0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15)
sim_list <- vector("list", length(differences))

for (i in seq_along(differences)) {
  ptm <- proc.time()
  print(i)
  print(differences[i])
  sim_df <-  make_sim_df(nsim = NSIM, diff = differences[i])
  sim_out <- mclapply(seq_len(NSIM), FUN = run_one_sim, sim_df = sim_df,
                      compiled_model = compiled_model,
                      formula = zoib_model, stanvars = stanvars,
                      mc.cores = 16, mc.preschedule = FALSE,
                      mc.allow.recursive = FALSE)
  save(sim_df, sim_out, 
       file = print(paste0("zoibr_sim_diff_", differences[i], "_v2.rda")))
  sim_list[[i]] <- sim_out
  print(proc.time() - ptm)
}

sim_res <- bind_rows(sim_list)

summary_sim <- sim_res %>% 
  filter(parameter == "prop_bet") %>% 
  group_by(diff) %>% 
  summarise(sig = mean(exclude_0), 
            sig_neg = mean(exclude_0 & .mean < 0),
            sig_pos = mean(exclude_0 & .mean > 0), 
            n = n(), 
            success = if_else(diff[1] == 0, sum(exclude_0), 
                              sum(exclude_0 & .mean < 0))) 
theme_set(theme_bw(base_size = 14))
library("binom")
summary_sim <- summary_sim %>% 
  mutate(
    upper = binom.profile(x = success, n = n)$upper,
    lower = binom.profile(x = success, n = n)$lower)


ggplot(summary_sim, aes(x = diff, y = sig)) +
  geom_hline(yintercept = 0.05, colour = "darkgrey") +
  geom_line() +
  #geom_point() +
  geom_pointrange(aes(ymax = upper, ymin = lower)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Difference in proportion bet", 
       y = "Probability of 95%-CI excl. 0")

ggsave("sensitivity.pdf", width = 5, height = 3.5)

#### second run of main simulation

#### main simulations:

NSIM <- 10000
differences <- seq(0, 0.1, by = 0.01)
sim_list <- vector("list", length(differences))

for (i in seq_along(differences)) {
  ptm <- proc.time()
  print(i)
  print(differences[i])
  sim_df <-  make_sim_df(nsim = NSIM, diff = differences[i])
  sim_out <- mclapply(seq_len(NSIM), FUN = run_one_sim, sim_df = sim_df,
                      compiled_model = compiled_model,
                      formula = zoib_model, stanvars = stanvars,
                      mc.cores = 16, mc.preschedule = FALSE,
                      mc.allow.recursive = FALSE)
  save(sim_df, sim_out, 
       file = print(paste0("zoibr_sim_diff_", differences[i], "_v3.rda")))
  sim_list[[i]] <- sim_out
  print(proc.time() - ptm)
}
sim_res <- bind_rows(sim_list)

#### load simulation data ###
sim_res <- map_dfr(list.files(pattern = "zoibr_sim_diff_.+_v3.rda"), ~{
  load(.)
  sim_out
})

summary_sim <- sim_res %>% 
  filter(parameter == "prop_bet") %>% 
  group_by(diff) %>% 
  summarise(sig = mean(exclude_0), 
            sig_neg = mean(exclude_0 & .mean < 0),
            sig_pos = mean(exclude_0 & .mean > 0), 
            n = n(), 
            success = if_else(diff[1] == 0, sum(exclude_0), 
                              sum(exclude_0 & .mean < 0))) %>% 
  mutate(sig = if_else(diff == 0, sig, sig_neg))
theme_set(theme_bw(base_size = 14))
library("binom")
summary_sim <- summary_sim %>% 
  mutate(
    upper = binom.profile(x = success, n = n)$upper,
    lower = binom.profile(x = success, n = n)$lower)


ggplot(summary_sim, aes(x = diff, y = sig)) +
  geom_hline(yintercept = 0.05, colour = "darkgrey") +
  geom_line() +
  #geom_point() +
  geom_pointrange(aes(ymax = upper, ymin = lower)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Difference in proportion bet", 
       y = "Probability of 95%-CI excl. 0") +
  scale_x_continuous(breaks = seq(0, 0.1, 0.02)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2))

ggsave("sensitivity_2.pdf", width = 5, height = 3.5)
