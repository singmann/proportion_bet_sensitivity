
calc_proportion_bet <- function(invec) {
  g <- invec[1]
  e <- invec[2]
  mu <- invec[3]
  g*e + g*(1-e)*mu
}

calc_proportion_bet_2 <- function(g, e, mu) {
  g*e + g*(1-e)*mu
}

g_e_from_weights <- function(weights) {
  g <- 1-weights[1]
  c(g = g, e = weights[2]/g)
}

weights_from_par <- function(par) {
  c((1-par[1]), par[1]*par[2], par[1]*(1-par[2]))
}

check_par <- function(par, case, diff) {
  if (any(par < 0) | any(par > 1)) {
    stop("illegal par for case ", case, " diff = ", diff, call. = FALSE)
  }
}

gen_one_data <- function(N_cond, g, e, mu, diff, phi) {
  
  which_case <- sample(1:3, 1, prob = c(0.25, 0.25, 0.5))
  #which_case <- 4
  
  prop_bet <- g*e + g*(1-e)*mu
  weights <- c((1-g), g*e, g*(1-e))
  
  if (which_case == 1) { ## achieve diff by moving g and e around
    diff_scale <- diff / (1 - weights[3]) 
    
    ## need to reparametrise g and e:
    ## x1 =  weights[1] / (1 - weights[3]) 
    ## x2 = weights[2] / (1 - weights[3])
    x1_a <- weights[1] / (1 - weights[3]) - diff_scale/2
    x1_b <- weights[1] / (1 - weights[3]) + diff_scale/2
    if ((x1_a < 0) | (x1_b > 1)) {
      stop("no case == 1 possible for diff = ", diff, call. = FALSE)
    }
    
    weights_a <- c(x1_a*(1 - weights[3]), (1-x1_a)*(1 - weights[3]), weights[3])
    weights_b <- c(x1_b*(1 - weights[3]), (1-x1_b)*(1 - weights[3]), weights[3])
    
    new_par_a <- c(g_e_from_weights(weights = weights_a), mu)
    new_par_b <- c(g_e_from_weights(weights = weights_b), mu)
    if (!isTRUE(all.equal(calc_proportion_bet(new_par_a) - 
                          calc_proportion_bet(new_par_b), 
                          diff, check.attributes = FALSE))) {
      stop("case 1 difference false for diff = ", diff, call. = FALSE)
    }
  } else if (which_case == 2) {
    diff_scale <- diff / (weights[3]) 
    
    new_par_a <- c(g, e, mu + diff_scale/2)
    check_par(new_par_a, case = 2, diff = diff)
    new_par_b <- c(g, e, mu - diff_scale/2)
    check_par(new_par_b, case = 2, diff = diff)
    
    if (!isTRUE(all.equal(calc_proportion_bet(new_par_a) - 
                          calc_proportion_bet(new_par_b), 
                          diff, check.attributes = FALSE))) {
      stop("case 2 difference false for diff = ", diff, call. = FALSE)
    }
  } else if (which_case == 3) {
    diff_scale_prop <- diff / (1 - weights[3])
    diff_scale_mu <- diff / weights[3]
    case_3_weight <- runif(1)
    
    x1_a <- weights[1] / (1 - weights[3]) - (diff_scale_prop * case_3_weight)/2
    x1_b <- weights[1] / (1 - weights[3]) + (diff_scale_prop * case_3_weight)/2
    if ((x1_a < 0) | (x1_b > 1)) {
      stop("no case == 3 possible for diff = ", diff, call. = FALSE)
    }
    
    weights_a <- c(x1_a*(1 - weights[3]), (1-x1_a)*(1 - weights[3]), 
                   weights[3])
    weights_b <- c(x1_b*(1 - weights[3]), (1-x1_b)*(1 - weights[3]), 
                   weights[3])
    
    new_par_a <- c(g_e_from_weights(weights = weights_a), 
                   mu + (diff_scale_mu * (1 - case_3_weight))/2)
    new_par_b <- c(g_e_from_weights(weights = weights_b), 
                   mu - (diff_scale_mu * (1- case_3_weight))/2)
    check_par(new_par_a, case = 3, diff = diff)
    check_par(new_par_b, case = 3, diff = diff)
    
    
    if (!isTRUE(all.equal(calc_proportion_bet(new_par_a) - 
                          calc_proportion_bet(new_par_b), 
                          diff, check.attributes = FALSE))) {
      stop("case 3 difference false for diff = ", diff, call. = FALSE)
    }
    
  } else if (which_case == 4) {
    ## currently not working, needs to update weights better
    case_3_weight <- runif(1)
    within_3_weight <- runif(1)
    
    
    diff_scale_mu <- diff / weights[3]
    
    w3_diff <- (diff_scale_mu * (1 - case_3_weight) * within_3_weight) / 2
    
    weights_a <- vector("numeric", 3)
    weights_b <- vector("numeric", 3)
    weights_a[3] <- weights[3] - w3_diff
    weights_b[3] <- weights[3] + w3_diff
    
    diff_scale_prop_a <- diff / (1 - weights_a[3])
    diff_scale_prop_b <- diff / (1 - weights_b[3])
    #prop_diff <- (diff_scale_prop * case_3_weight)/2
    
    x1_a <- weights[1] / (1 - weights_a[3]) - 
      (diff_scale_prop_a * case_3_weight)/2
    x1_b <- weights[1] / (1 - weights_b[3]) + 
      (diff_scale_prop_b * case_3_weight)/2
    if ((x1_a < 0) | (x1_b > 1)) {
      stop("no case == 4 possible for diff = ", diff, call. = FALSE)
    }
    
    weights_a[1:2] <- c(x1_a*(1 - weights_a[3]), (1-x1_a)*(1 - weights_a[3]))
    weights_b[1:2] <- c(x1_b*(1 - weights_b[3]), (1-x1_b)*(1 - weights_b[3]))
    
    new_par_a <- c(g_e_from_weights(weights = weights_a), 
                   mu + ((diff / weights_a[3]) * 
                           (1 - case_3_weight) * 
                           (1 - within_3_weight))/2)
    new_par_b <- c(g_e_from_weights(weights = weights_b), 
                   mu - ((diff / weights_b[3]) * 
                           (1 - case_3_weight) * 
                           (1 - within_3_weight))/2)
    check_par(new_par_a, case = 3, diff = diff)
    check_par(new_par_b, case = 3, diff = diff)
    
    if (!isTRUE(all.equal(calc_proportion_bet(new_par_a) - 
                          calc_proportion_bet(new_par_b), 
                          diff, check.attributes = FALSE))) {
      stop("case 3 difference false for diff = ", diff, call. = FALSE)
    }
  }
  
  new_par_a
  new_par_b
  data_out <- bind_rows(
    gen_from_par(N = N_cond, par = new_par_a, phi = phi, cond = "A"),
    gen_from_par(N = N_cond, par = new_par_b, phi = phi, cond = "B")
  )
  data_out$condition <- factor(data_out$condition)
  data_out
}

gen_from_par <- function(N, par, phi, condition) {
  
  weights <- weights_from_par(par)
  dmnom <- rmultinom(1, size = N, prob = weights)
  data_out <- data.frame(
    id = factor(seq_len(N)),
    response = c(rep(0, dmnom[1,1]), 
                 rep(1, dmnom[2,1]), 
                 rep(NA_real_, dmnom[3,1])),
    condition = condition
  )
  a <- par[3] * phi
  b <- (1 - par[3]) * phi
  
  dbeta <- rbeta(dmnom[3,1], shape1 = a, shape2 = b)
  data_out$response[is.na(data_out$response)] <- dbeta
  data_out
}

draws_to_summary <- function(draws_df) {
  res_mu <- draws_df %>% 
    mean_qi(`b[2]`) %>% 
    rename(.mean = `b[2]`) %>% 
    mutate(parameter = "mu")
  res_g <- draws_df %>% 
    mean_qi(`b_zipp[2]`) %>% 
    rename(.mean = `b_zipp[2]`) %>% 
    mutate(parameter = "g")
  res_e <- draws_df %>% 
    mean_qi(`b_coi[2]`) %>% 
    rename(.mean = `b_coi[2]`) %>% 
    mutate(parameter = "e")
  
  draws2 <- draws_df %>% 
    mutate(mu_a = plogis(`b[1]`), 
           mu_b = plogis(`b[1]` + `b[2]`),
           g_a = plogis(`b_zipp[1]`),
           g_b = plogis(`b_zipp[1]` + `b_zipp[2]`),
           e_a = plogis(`b_coi[1]`),
           e_b = plogis(`b_coi[1]` + `b_coi[2]`)) %>% 
    mutate(
      prop_bet_a = (g_a * e_a) + (g_a * (1-e_a) * mu_a),
      prop_bet_b = (g_b * e_b) + (g_b * (1-e_b) * mu_b)
    ) %>% 
    mutate(
      diff = prop_bet_b - prop_bet_a
    )
  res_prop_bet <- draws2 %>% 
    mean_qi(diff) %>% 
    rename(.mean = diff) %>% 
    mutate(parameter = "prop_bet")
  bind_rows(
    res_prop_bet,
    res_mu, 
    res_e, 
    res_g
  )
}

fit_one_dat <- function(data, compiled_model, 
                        formula, stanvars,
                        chains = 6, 
                        parallel_chains = 1,
                        refresh = 0, 
                        iter_warmup = 1000, iter_sampling = 1000) {
  standat <- make_standata(formula = formula, data = data, stanvars = stanvars)
  fit <- compiled_model$sample(
    data = standat, 
    chains = chains, 
    parallel_chains = parallel_chains,
    refresh = refresh, 
    iter_warmup = iter_warmup, iter_sampling = iter_sampling
  )
  draws <- fit$draws(format = "df")
  draws_to_summary(draws)
}

make_sim_df <- function(nsim, diff, 
                        N_cond = 500, 
                        g = 0.85, e = 0.2, mu = 0.4, phi = exp(1.11)) {
  tibble(
    index = seq_len(nsim), 
    diff = diff, 
    n = N_cond, 
    g = g, e = e, mu = mu, phi = phi
  )
}

run_one_sim <- function(index, sim_df, compiled_model, 
                        formula, stanvars,
                        chains = 6, 
                        parallel_chains = 1,
                        refresh = 0, 
                        iter_warmup = 1000, iter_sampling = 1000) {
  
  tmpdat <- gen_one_data(N_cond = sim_df$n[index], 
                         g = sim_df$g[index], 
                         e = sim_df$e[index], 
                         mu = sim_df$mu[index], 
                         diff = sim_df$diff[index], phi = sim_df$phi[index])
  obs_means <- tapply(tmpdat$response, tmpdat$condition, mean)
  cur_df <- sim_df[index,]
  cur_df$obs_diff <- obs_means[2] - obs_means[1]
  cur_df <- cur_df %>% 
    select(obs_diff, everything())
  fit_res <- fit_one_dat(data = tmpdat, compiled_model = compiled_model, 
              formula = zoib_model, stanvars = stanvars, 
              chains = chains, 
              parallel_chains = parallel_chains,
              refresh = refresh, 
              iter_warmup = iter_warmup, iter_sampling = iter_sampling)
  fit_res$index <- index
  fit_res %>% 
    select(index, .mean:.upper, parameter) %>% 
    mutate(exclude_0 = sign(.lower) == sign(.upper)) %>% 
    left_join(cur_df, by = "index")
}
