
# DEPENDENT PY PROCESS MIXTURE MODEL

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(tidyverse)
})

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if(length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  return(getwd())
}

PROJECT_DIR <- get_script_dir()


# Funzioni ausiliarie

## Ho lavarato sempre con log value per questioni di stabilità numerica

## log dnsità dei theta
log_beta_density <- function(x, a, b){
  if(any(a <= 0 | b <= 0)) return(rep(-Inf, length(x)))
  epsilon <- 1e-10
  x_safe <- pmax(epsilon, pmin(x, 1 - epsilon))
  val <- (a-1)*log(x_safe) + (b-1)*log(1-x_safe) + lgamma(a+b) - lgamma(a) - lgamma(b)
  return(val)
}

posterior_indices <- function(n_total, burnin_pct) {
  if(n_total < 2) stop("At least two MCMC iterations are required.")
  burnin_pct <- pmax(0, pmin(burnin_pct, 0.99))
  burnin <- floor(n_total * burnin_pct)
  return(seq.int(burnin + 1, n_total))
}

## creo quest funzione principalmente per sicurezza numerica
rinvgamma <- function(n, shape, rate) {
  shape <- pmax(shape, 0.001) ## Sempre per stabilità numerico ho dovuto inserire alcune protezioni
  rate  <- pmax(rate, 0.001) 
  vals <- tryCatch({
    rgamma(n, shape = shape, rate = rate)
  }, warning = function(w) {
    rep(1.0, n) 
  })
  vals <- pmax(vals, 1e-10)
  return(1 / vals)
}

## Calcolo Pesi
compute_dependent_weights <- function(V0, V1, V2){ ##Per semplicità implemento solo H1 ma 
  K <- length(V0)                                  ## è facilmente estendibile ad H2
  pi1 <- numeric(K); pi2 <- numeric(K)
  rem1 <- 1; rem2 <- 1
  
  for(k in 1:K){
    S1k <- V0[k] * V1[k]
    S2k <- V0[k] * V2[k]
    
    pi1[k] <- S1k * rem1
    pi2[k] <- S2k * rem2
    
    rem1 <- rem1 * (1 - S1k)
    rem2 <- rem2 * (1 - S2k)
  }
  return(list(pi1=pi1, pi2=pi2))
}

## Inizializzazione Sticks
initialize_sticks <- function(K, theta1, theta2, c){
  V0 <- rbeta(K, 1-c, theta1)
  V1 <- rbeta(K, 1-c + theta1, theta2 + c*(0:(K-1)))
  V2 <- rbeta(K, 1-c + theta1, theta2 + c*(0:(K-1)))
  return(list(V0=V0, V1=V1, V2=V2))
}

# Aggiornamenti MCMC

## Aggiornamento cluster 
update_clusters_anova <- function(Y1, Y2, z1, z2, 
                                  mu_0, sigma_0, # Parametri Comuni
                                  mu_1, sigma_1, # Parametri Specifici G1
                                  mu_2, sigma_2, # Parametri Specifici G2
                                  s0_sq=100, si_sq=100, # Varianze Prior
                                  lambda=1, epsilon=1   # Parametri IG Prior
){
  
  K <- length(mu_0)
  
  for(k in 1:K){
    idx1 <- which(z1 == k); n1 <- length(idx1)
    idx2 <- which(z2 == k); n2 <- length(idx2)
    
    y1k <- if(n1 > 0) Y1[idx1] else numeric(0)
    y2k <- if(n2 > 0) Y2[idx2] else numeric(0)
    
    ## Aggiornamento Specifico: la formula segue dal paper e dalla coniugatezza del modello N-IG
    if(n1 > 0){
      prec_post_mu1 <- (1/si_sq) + (n1 / (sigma_1[k] * sigma_0[k]))
      var_post_mu1 <- 1 / prec_post_mu1
      resid_1 <- sum(y1k - mu_0[k]) 
      mean_post_mu1 <- var_post_mu1 * (resid_1 / (sigma_1[k] * sigma_0[k]))
      mu_1[k] <- rnorm(1, mean_post_mu1, sqrt(var_post_mu1))
      
      shape_post_sig1 <- lambda/2 + n1/2
      sse_1 <- sum((y1k - (mu_0[k] + mu_1[k]))^2)
      rate_post_sig1 <- lambda/2 + sse_1 / (2 * sigma_0[k])
      sigma_1[k] <- rinvgamma(1, shape_post_sig1, rate_post_sig1)
    } else {
      mu_1[k] <- rnorm(1, 0, sqrt(si_sq))
      sigma_1[k] <- rinvgamma(1, lambda/2, lambda/2)
    }
    
    ## Aggiornamento Specifico
    if(n2 > 0){
      prec_post_mu2 <- (1/si_sq) + (n2 / (sigma_2[k] * sigma_0[k]))
      var_post_mu2 <- 1 / prec_post_mu2
      resid_2 <- sum(y2k - mu_0[k])
      mean_post_mu2 <- var_post_mu2 * (resid_2 / (sigma_2[k] * sigma_0[k]))
      mu_2[k] <- rnorm(1, mean_post_mu2, sqrt(var_post_mu2))
      
      shape_post_sig2 <- lambda/2 + n2/2
      sse_2 <- sum((y2k - (mu_0[k] + mu_2[k]))^2)
      rate_post_sig2 <- lambda/2 + sse_2 / (2 * sigma_0[k])
      sigma_2[k] <- rinvgamma(1, shape_post_sig2, rate_post_sig2)
    } else {
      mu_2[k] <- rnorm(1, 0, sqrt(si_sq))
      sigma_2[k] <- rinvgamma(1, lambda/2, lambda/2)
    }
    
    ## Aggiornamento Comune media 
    term1 <- if(n1>0) n1 / (sigma_1[k] * sigma_0[k]) else 0
    term2 <- if(n2>0) n2 / (sigma_2[k] * sigma_0[k]) else 0
    prec_post_mu0 <- (1/s0_sq) + term1 + term2
    var_post_mu0 <- 1 / prec_post_mu0
    
    sum_res_1 <- if(n1>0) sum(y1k - mu_1[k]) / (sigma_1[k] * sigma_0[k]) else 0
    sum_res_2 <- if(n2>0) sum(y2k - mu_2[k]) / (sigma_2[k] * sigma_0[k]) else 0
    
    mean_post_mu0 <- var_post_mu0 * (sum_res_1 + sum_res_2)
    mu_0[k] <- rnorm(1, mean_post_mu0, sqrt(var_post_mu0))
    
    ## Aggiornamento Comune var
    shape_post_sig0 <- epsilon/2 + (n1 + n2)/2
    sse_comb_1 <- if(n1>0) sum((y1k - (mu_0[k] + mu_1[k]))^2) / sigma_1[k] else 0
    sse_comb_2 <- if(n2>0) sum((y2k - (mu_0[k] + mu_2[k]))^2) / sigma_2[k] else 0
    rate_post_sig0 <- epsilon/2 + (sse_comb_1 + sse_comb_2) / 2
    sigma_0[k] <- rinvgamma(1, shape_post_sig0, rate_post_sig0)
  }
  
  return(list(mu_0=mu_0, sigma_0=sigma_0, mu_1=mu_1, sigma_1=sigma_1, mu_2=mu_2, sigma_2=sigma_2))
}

## Full conditional per MH Stick eq 38
log_Qk <- function(v_vec, A_k, B_k, psi, k, tol = 1e-9) {
  v0 <- v_vec[1]; v1 <- v_vec[2]; v2 <- v_vec[3]
  theta1 <- psi[1]; theta2 <- psi[2]; c <- psi[3]
  
  if(any(v_vec <= tol | v_vec >= (1 - tol))) return(-Inf) ## sicurezza numerica
  
  term_v0 <- (-c + A_k[1] + A_k[2]) * log(v0) + (theta1 - 1) * log(1 - v0) ##lavoro in log
  term_v1 <- (A_k[1] + theta1 - c) * log(v1) + (theta2 + c*(k-1) - 1) * log(1 - v1)
  term_v2 <- (A_k[2] + theta1 - c) * log(v2) + (theta2 + c*(k-1) - 1) * log(1 - v2)
  
  arg_c1 <- 1 - v0 * v1
  arg_c2 <- 1 - v0 * v2
  if(arg_c1 <= 0 || arg_c2 <= 0) return(-Inf) ## sicurezza numerica
  
  val <- term_v0 + term_v1 + term_v2 + B_k[1] * log(arg_c1) + B_k[2] * log(arg_c2) ## produttoria
  return(if(is.na(val)) -Inf else val)
}

## Aggiornamento Stick MH
update_sticks_MH <- function(V_list, D_list, psi, scale_prop = 0.3) {
  V0 <- V_list$V0; V1 <- V_list$V1; V2 <- V_list$V2
  K_curr <- length(V0)
  new_V0 <- V0
  new_V1 <- V1
  new_V2 <- V2 ## qui evito else dopo il test
  accepted <- 0
  
  for(k in 1:(K_curr-1)) { 
    A_k <- c(sum(D_list[[1]] == k), sum(D_list[[2]] == k))  ##numero di componenti in K
    B_k <- c(sum(D_list[[1]] > k), sum(D_list[[2]] > k)) ## numero dei componenti in un cluster maggiore K
    v_curr <- c(V0[k], V1[k], V2[k])
    
    logit_curr <- qlogis(pmax(pmin(v_curr, 1-1e-6), 1e-6))  ## lavriamo in Log per stabilità ma con sicurezze numeriche
    logit_prop <- rnorm(3, mean = logit_curr, sd = scale_prop)  ## proponiamo nuovi stick
    v_prop <- plogis(logit_prop) ##calcolo cdf logis della proposta (sicurezza numerica)
    
    log_q_curr <- log_Qk(v_curr, A_k, B_k, psi, k) ##posterior di V date le allocazioni, i parametri e i K
    log_q_prop <- log_Qk(v_prop, A_k, B_k, psi, k)
    
    log_jac_curr <- sum(log(pmax(v_curr, 1e-10)) + log(pmax(1 - v_curr, 1e-10))) ## abbiamo trasformato una r.v. serve jacobiano
    log_jac_prop <- sum(log(pmax(v_prop, 1e-10)) + log(pmax(1 - v_prop, 1e-10))) ## per riportare alla scala originale
    
    mh_ratio <- (log_q_prop + log_jac_prop) - (log_q_curr + log_jac_curr) ## questo è il classico MH ratio però
    
    if(is.finite(mh_ratio) && log(runif(1)) < mh_ratio) {   ## MH test in spazio log (log è trasformato monotona per valori positivi)
      new_V0[k] <- v_prop[1]; new_V1[k] <- v_prop[2]; new_V2[k] <- v_prop[3] 
      accepted <- accepted + 1
    }
  }
  new_V0[K_curr] <- 1; new_V1[K_curr] <- 1; new_V2[K_curr] <- 1
  
  result <- list(V0 = new_V0, V1 = new_V1, V2 = new_V2)
  attr(result, "acc_rate") <- accepted / (K_curr - 1)
  return(result)
}

## MH per Iperparametri
update_theta1 <- function(theta1, theta2, c, V0, V1, V2, a=0.01, b=0.01){
  theta1_star <- theta1 * exp(rnorm(1,0,0.1)) ##Log-normal RW
  K_eff <- length(V0) - 1; if(K_eff < 1) return(theta1) ## siamo in log sarebbe il prodotto dato indipendenza
                                                        ## diventa quinsi somma
  idx <- 1:K_eff
  
  logpost_star <- sum(log_beta_density(V0[idx], 1-c, theta1_star)) +
    sum(log_beta_density(V1[idx], 1-c+theta1_star, theta2 + c*(idx-1))) +
    sum(log_beta_density(V2[idx], 1-c+theta1_star, theta2 + c*(idx-1))) +
    dgamma(theta1_star, a, b, log=TRUE) ## sicurezza numerica
  
  logpost_curr <- sum(log_beta_density(V0[idx], 1-c, theta1)) +
    sum(log_beta_density(V1[idx], 1-c+theta1, theta2 + c*(idx-1))) +
    sum(log_beta_density(V2[idx], 1-c+theta1, theta2 + c*(idx-1))) +
    dgamma(theta1, a, b, log=TRUE)
  
  if(!is.finite(logpost_star) || !is.finite(logpost_curr)) return(theta1) ## MH
  log_accept <- (logpost_star + log(theta1_star)) - (logpost_curr + log(theta1))
  if(log(runif(1)) < log_accept) return(theta1_star) else return(theta1)
}

## Uguale a prima cambiano solo addendi perche theeta2 non appare in V0
update_theta2 <- function(theta1, theta2, c, V1, V2, a=0.01, b=0.01){
  theta2_star <- theta2 * exp(rnorm(1,0,0.1))
  K_eff <- length(V1) - 1; idx <- 1:K_eff
  
  logpost_star <- sum(log_beta_density(V1[idx], 1-c+theta1, theta2_star + c*(idx-1))) +
    sum(log_beta_density(V2[idx], 1-c+theta1, theta2_star + c*(idx-1))) +
    dgamma(theta2_star, a, b, log=TRUE)  ## Qui cambia perhè theta2 non appare in V0
  
  logpost_curr <- sum(log_beta_density(V1[idx], 1-c+theta1, theta2 + c*(idx-1))) +
    sum(log_beta_density(V2[idx], 1-c+theta1, theta2 + c*(idx-1))) +
    dgamma(theta2, a, b, log=TRUE)
  
  if(!is.finite(logpost_star) || !is.finite(logpost_curr)) return(theta2)
  log_accept <- (logpost_star + log(theta2_star)) - (logpost_curr + log(theta2))
  if(log(runif(1)) < log_accept) return(theta2_star) else return(theta2)
} 

update_c <- function(theta1, theta2, c, V0, V1, V2){
  logit_c <- log(c/(1-c)) ## inversa per essere in spazio R
  logit_c_star <- logit_c + rnorm(1,0,0.05) ## rw 
  c_star <- exp(logit_c_star)/(1+exp(logit_c_star)) ## trasf log
  if(c_star <= 1e-4 || c_star >= (1 - 1e-4)) return(c) ## sicurezza numerica
  
  K_eff <- length(V0) - 1; idx <- 1:K_eff
  
  logpost_star <- sum(log_beta_density(V0[idx], 1-c_star, theta1)) +
    sum(log_beta_density(V1[idx], 1-c_star+theta1, theta2 + c_star*(idx-1))) +
    sum(log_beta_density(V2[idx], 1-c_star+theta1, theta2 + c_star*(idx-1)))
  ## sarebbe prodotto ma data log diventa somma
  
  logpost_curr <- sum(log_beta_density(V0[idx], 1-c, theta1)) +
    sum(log_beta_density(V1[idx], 1-c+theta1, theta2 + c*(idx-1))) +
    sum(log_beta_density(V2[idx], 1-c+theta1, theta2 + c*(idx-1)))
  
  if(!is.finite(logpost_star) || !is.finite(logpost_curr)) return(c)
  log_jac_star <- log(c_star) + log1p(-c_star)
  log_jac_curr <- log(c) + log1p(-c)
  log_accept <- (logpost_star + log_jac_star) - (logpost_curr + log_jac_curr)
  if(log(runif(1)) < log_accept) return(c_star) else return(c)
}

# Ciclo principale


run_DPYP_ANOVA <- function(
    Y1, Y2,
    K_max = 15,
    n_iter = 6000,
    burnin_pct = 0.50,
    s0_sq = 100, si_sq = 0.5,  ## var a priori per la generazione di media e var dei cluster
    lambda = 1, epsilon = 1, 
    theta1_init = 1.0, 
    theta2_init = 1.0, 
    c_init = 0.3,
    PYP = TRUE, ## inseiro questo comando per fare DDP
    seed = NULL
) {
  
  if(!is.null(seed)) set.seed(seed)
  N_1 <- length(Y1); N_2 <- length(Y2)
  
  ## Init Parametri
  mu_0 <- rnorm(K_max, 0, sqrt(s0_sq)); sigma_0 <- rinvgamma(K_max, epsilon/2, epsilon/2)
  mu_1 <- rnorm(K_max, 0, sqrt(si_sq)); sigma_1 <- rinvgamma(K_max, lambda/2, lambda/2)
  mu_2 <- rnorm(K_max, 0, sqrt(si_sq)); sigma_2 <- rinvgamma(K_max, lambda/2, lambda/2)
  
  z_1 <- sample(1:K_max, N_1, replace=TRUE)
  z_2 <- sample(1:K_max, N_2, replace=TRUE)
  theta1 <- theta1_init; theta2 <- theta2_init; c <- c_init
  current_scale_prop <- 0.5 
  
  sticks <- initialize_sticks(K_max, theta1, theta2, c)
  V0 <- sticks$V0; V1 <- sticks$V1; V2 <- sticks$V2
  
  ## Storage
  all_pi1 <- array(NA, dim=c(n_iter, K_max)); all_pi2 <- array(NA, dim=c(n_iter, K_max))
  all_mu_1_tot <- array(NA, dim=c(n_iter, K_max)); all_sigma_1_tot <- array(NA, dim=c(n_iter, K_max))
  all_mu_2_tot <- array(NA, dim=c(n_iter, K_max)); all_sigma_2_tot <- array(NA, dim=c(n_iter, K_max))
  all_z_combined <- matrix(NA, nrow=n_iter, ncol=N_1+N_2)
  
  all_theta1 <- numeric(n_iter); all_theta2 <- numeric(n_iter); all_c <- numeric(n_iter)
  all_n_clusters <- matrix(NA, n_iter, 2); all_acc_rate <- numeric(n_iter)
  
  burnin_iter <- floor(n_iter * pmax(0, pmin(burnin_pct, 0.99)))
  
  ## Gibbs Loop
  for(it in 1:n_iter){
    
    ## Stick-breaking 
    weights <- compute_dependent_weights(V0, V1, V2)
    pi1 <- weights$pi1; pi2 <- weights$pi2
    
    ## Likelihood e Allocazioni
    tot_mu1 <- mu_0 + mu_1; tot_sd1 <- sqrt(sigma_0 * sigma_1)
    tot_mu2 <- mu_0 + mu_2; tot_sd2 <- sqrt(sigma_0 * sigma_2)
    
    log_probs_1 <- matrix(NA, nrow=N_1, ncol=K_max)
    for(k in 1:K_max) log_probs_1[, k] <- log(pi1[k] + 1e-100) + dnorm(Y1, tot_mu1[k], tot_sd1[k], log=TRUE)
    for(i in 1:N_1) { lp <- log_probs_1[i, ]; z_1[i] <- sample(1:K_max, 1, prob=exp(lp - max(lp))) }
    
    log_probs_2 <- matrix(NA, nrow=N_2, ncol=K_max)
    for(k in 1:K_max) log_probs_2[, k] <- log(pi2[k] + 1e-100) + dnorm(Y2, tot_mu2[k], tot_sd2[k], log=TRUE)
    for(i in 1:N_2) { lp <- log_probs_2[i, ]; z_2[i] <- sample(1:K_max, 1, prob=exp(lp - max(lp))) }
    
    ## Aggiornamento Atomi (ANOVA)
    atoms <- update_clusters_anova(Y1, Y2, z_1, z_2, mu_0, sigma_0, mu_1, sigma_1, mu_2, sigma_2, s0_sq, si_sq, lambda, epsilon)
    mu_0 <- atoms$mu_0; sigma_0 <- atoms$sigma_0
    mu_1 <- atoms$mu_1; sigma_1 <- atoms$sigma_1
    mu_2 <- atoms$mu_2; sigma_2 <- atoms$sigma_2
    
    ## Aggiornamento Sticks
    V_list <- list(V0 = V0, V1 = V1, V2 = V2); D_list <- list(z_1, z_2)
    V_update <- update_sticks_MH(V_list, D_list, c(theta1, theta2, c), scale_prop = current_scale_prop)
    V0 <- V_update$V0; V1 <- V_update$V1; V2 <- V_update$V2
    all_acc_rate[it] <- attr(V_update, "acc_rate")
    
    if(it <= burnin_iter && it %% 50 == 0) {
      recent_acc <- mean(all_acc_rate[max(1, it-49):it])
      if(recent_acc > 0.44) current_scale_prop <- current_scale_prop * 1.1
      if(recent_acc < 0.234) current_scale_prop <- current_scale_prop * 0.9
    }
    
    ## Aggiornamento Iperparametri
    theta1 <- update_theta1(theta1, theta2, c, V0, V1, V2)
    theta2 <- update_theta2(theta1, theta2, c, V1, V2)
    
    if(PYP) {
      c <- update_c(theta1, theta2, c, V0, V1, V2)
    }
    
    ## Salvataggio
    all_pi1[it, ] <- pi1; all_pi2[it, ] <- pi2
    all_mu_1_tot[it, ] <- mu_0 + mu_1; all_sigma_1_tot[it, ] <- sqrt(sigma_0 * sigma_1)
    all_mu_2_tot[it, ] <- mu_0 + mu_2; all_sigma_2_tot[it, ] <- sqrt(sigma_0 * sigma_2)
    all_z_combined[it, ] <- c(z_1, z_2)
    all_theta1[it] <- theta1; all_theta2[it] <- theta2; all_c[it] <- c
    all_n_clusters[it, ] <- c(length(unique(z_1)), length(unique(z_2)))
    
    if(it %% 1000 == 0) cat(sprintf("Iter %d | Acc: %.2f | c: %.4f\n", it, all_acc_rate[it], c))
  }
  
  return(list(
    data = list(Y1=Y1, Y2=Y2),
    config = list(n_iter=n_iter),
    trace = list(
      pi1 = all_pi1, pi2 = all_pi2,
      mu_1 = all_mu_1_tot, sigma_1 = all_sigma_1_tot,
      mu_2 = all_mu_2_tot, sigma_2 = all_sigma_2_tot,
      Z = all_z_combined,
      theta1 = all_theta1, theta2 = all_theta2, c = all_c,
      n_clusters = all_n_clusters, acc_rate = all_acc_rate
    )
  ))
}
# Grafici


plot_posterior_mean <- function(result, burnin_pct = 0.8){
  
  # --- CATTURA IL NOME DELL'OGGETTO ---
  obj_name <- deparse(substitute(result))
  
  Y1 <- result$data$Y1
  Y2 <- result$data$Y2
  n_total <- result$config$n_iter
  idx_keep <- posterior_indices(n_total, burnin_pct)
  
  y_grid_1 <- seq(min(Y1)-3, max(Y1)+3, length.out=500)
  y_grid_2 <- seq(min(Y2)-3, max(Y2)+3, length.out=500)
  
  compute_mean_density <- function(pi_trace, mu_trace, sigma_trace, grid, idxs) {
    dens_accum <- numeric(length(grid))
    for(i in seq_along(idxs)) {
      real_idx <- idxs[i]
      pi_t <- pi_trace[real_idx, ]
      mu_t <- mu_trace[real_idx, ]
      sig_t <- sigma_trace[real_idx, ]
      active <- which(pi_t > 1e-5)
      dens_t <- numeric(length(grid))
      for(k in active) {
        dens_t <- dens_t + pi_t[k] * dnorm(grid, mean=mu_t[k], sd=sig_t[k])
      }
      dens_accum <- dens_accum + dens_t
    }
    return(dens_accum / length(idxs))
  }
  
  dens_mean_1 <- compute_mean_density(result$trace$pi1, result$trace$mu_1, result$trace$sigma_1, y_grid_1, idx_keep)
  dens_mean_2 <- compute_mean_density(result$trace$pi2, result$trace$mu_2, result$trace$sigma_2, y_grid_2, idx_keep)
  
  df_hist <- rbind(data.frame(Val=Y1, Group="Group 1"), data.frame(Val=Y2, Group="Group 2"))
  df_dens <- rbind(data.frame(Val=y_grid_1, Dens=dens_mean_1, Group="Group 1"),
                   data.frame(Val=y_grid_2, Dens=dens_mean_2, Group="Group 2"))
  
  cols_fill <- c("Group 1" = "#6BAED6", "Group 2" = "#FD8D3C") 
  cols_line <- c("Group 1" = "#08306B", "Group 2" = "#D95F02")   
  p <- ggplot() +
	    geom_histogram(data=df_hist, aes(x=Val, y=after_stat(density), fill=Group),
	                   bins=100, color="white", alpha=0.65) +
    geom_line(data=df_dens, aes(x=Val, y=Dens, color=Group), size=1.2) +
    facet_wrap(~Group, scales="free") +
    scale_fill_manual(values=cols_fill) +
    scale_color_manual(values=cols_line) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      strip.text = element_text(face="bold", size=14),
      plot.title = element_text(size=16, face="bold")
    ) +
    labs(
      title = paste0("Posterior Mean Density - ", obj_name), 
      x = "Value",
      y = "Density"
    )
  
  print(p)
  return(p)
}

plot_diagnostics <- function(trace, title = NULL) {
  
  obj_name <- deparse(substitute(trace))
  
  if(is.null(title)) {
    final_title <- paste0("MCMC Diagnostics - ", obj_name)
  } else {
    final_title <- paste0(title, " (", obj_name, ")")
  }
  
  n_iter <- length(trace$theta1)
  seq_iter <- 1:n_iter
  
  df_params <- rbind(
    data.frame(Iter = seq_iter, Value = trace$theta1[seq_iter], Param = "theta1"),
    data.frame(Iter = seq_iter, Value = trace$theta2[seq_iter], Param = "theta2"),
    data.frame(Iter = seq_iter, Value = trace$c[seq_iter],      Param = "c")
  )
  
  df_params$Param <- factor(df_params$Param, levels = c("theta1", "theta2", "c"))
  
  p1 <- ggplot(df_params, aes(x = Iter, y = Value, color = Param)) +
    geom_line(alpha = 0.8, size = 0.6) +
    facet_wrap(~ Param, scales = "free_y", ncol = 1, strip.position = "top") +
    scale_color_manual(values = c("salmon", "#4daf4a", "#377eb8")) + 
    labs(title = "DPY Concentration Parameters & Dependency", y = "Value", x = NULL) +
    theme_minimal() +
    theme(legend.position = "none", strip.text = element_text(face = "bold", size = 10))
  
  df_k <- rbind(
    data.frame(Iter = seq_iter, K = trace$n_clusters[seq_iter, 1], Group = "Group 1"),
    data.frame(Iter = seq_iter, K = trace$n_clusters[seq_iter, 2], Group = "Group 2")
  )
  
  p2 <- ggplot(df_k, aes(x = Iter, y = K, color = Group)) +
    geom_line(alpha = 0.8, size = 0.6) +
    scale_color_manual(values = c("steelblue", "darkorange")) +
    labs(title = "Number of Active Components", y = "K", x = NULL) +
    theme_minimal() +
    theme(legend.position = c(0.9, 0.85), 
          legend.background = element_rect(fill = "white", color = NA))
  
  df_acc <- data.frame(Iter = seq_iter, Rate = trace$acc_rate[seq_iter])
  
  p3 <- ggplot(df_acc, aes(x = Iter)) +
    geom_line(aes(y = Rate), color = "darkgreen", alpha = 0.3, size = 0.3) + 
    geom_hline(yintercept = c(0.234, 0.44), linetype = "dashed", color = "red") +
    ylim(0, 1) +
    labs(title = "MH Acceptance Rate (Sticks)", y = "Acceptance Rate", x = "Iteration") +
    theme_minimal()
  
  grid.arrange(p1, p2, p3, ncol = 1, heights = c(1.3, 1, 1), top = final_title)
}

scatter <- function(result, burnin_pct = 0.8, n_top = 6) {
  
  obj_name <- deparse(substitute(result))
  
  n_total <- result$config$n_iter
  idx_keep <- posterior_indices(n_total, burnin_pct)
  
  pi1 <- result$trace$pi1[idx_keep, ]
  pi2 <- result$trace$pi2[idx_keep, ]
  mu1_tr <- result$trace$mu_1[idx_keep, ]
  mu2_tr <- result$trace$mu_2[idx_keep, ]
  sig1_tr <- result$trace$sigma_1[idx_keep, ]
  sig2_tr <- result$trace$sigma_2[idx_keep, ]
  
  mean_weights <- colMeans(pi1 + pi2) / 2
  n_top <- min(n_top, length(mean_weights))
  top_indices <- order(mean_weights, decreasing = TRUE)[1:n_top]
  
  df_labels_list <- list()
  
  cat("\n=== STATISTICHE CLUSTER (Media Posteriori) ===\n")
  cat(sprintf("%-10s | %-15s | %-15s\n", "Cluster", "Gruppo 1 (Mu/Var)", "Gruppo 2 (Mu/Var)"))
  cat(rep("-", 50), "\n")
  
  for(k in top_indices) {
    m1 <- mean(mu1_tr[, k]); v1 <- mean(sig1_tr[, k]^2) 
    m2 <- mean(mu2_tr[, k]); v2 <- mean(sig2_tr[, k]^2)
    
    cat(sprintf("Cluster %-2d | M: %5.2f V: %4.2f | M: %5.2f V: %4.2f\n", 
                k, m1, v1, m2, v2))
    
    # --- MODIFICA 1: Costruzione stringa per 'parse = TRUE' ---
    # atop(a, b): mette 'a' sopra 'b' (per andare a capo)
    # ~: aggiunge uno spazio
    # ==: disegna il simbolo uguale
    # mu / sigma: disegna le lettere greche
    # ^2: fa l'elevamento a potenza grafico
    
    label_text <- sprintf(
      "atop(G[1]:~mu==%s~~sigma^2==%s, G[2]:~mu==%s~~sigma^2==%s)",
      round(m1, 1), round(v1, 2),
      round(m2, 1), round(v2, 2)
    )
    
    df_labels_list[[length(df_labels_list)+1]] <- data.frame(
      Cluster = paste0("Cluster ", k),
      Label = label_text
    )
  }
  df_labels <- do.call(rbind, df_labels_list)
  df_labels$Cluster <- factor(df_labels$Cluster, levels = paste0("Cluster ", top_indices))
  
  idx_sub <- sample(seq_along(idx_keep), min(2000, length(idx_keep)))
  df_list <- list()
  
  for(i in seq_along(idx_sub)) {
    curr_idx <- idx_sub[i]
    for(k in top_indices) {
      w1 <- pi1[curr_idx, k]
      w2 <- pi2[curr_idx, k]
      if(w1 > 0.05 || w2 > 0.05) {
        df_list[[length(df_list)+1]] <- data.frame(
          Pi1 = w1,
          Pi2 = w2,
          Cluster = paste0("Cluster ", k)
        )
      }
    }
  }
  
  df_plot <- do.call(rbind, df_list)
  df_plot$Cluster <- factor(df_plot$Cluster, levels = paste0("Cluster ", top_indices))
  
  p <- ggplot(df_plot, aes(x = Pi1, y = Pi2)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60", size = 0.8) +
    geom_point(aes(color = Cluster), alpha = 0.5, size = 2) +
    
    # --- MODIFICA 2: parse = TRUE ---
    # Questo dice a ggplot di interpretare il testo come formula matematica
    geom_label(data = df_labels, aes(label = Label), 
               x = 0, y = 1, hjust = 0, vjust = 1, 
               size = 3, alpha = 0.8, 
               parse = TRUE,  # <--- FONDAMENTALE
               family = "sans") + # Usa font standard di sistema
    
    facet_wrap(~Cluster, ncol = 3) + 
    scale_color_brewer(palette = "Set1") +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      strip.text = element_text(face = "bold", size = 12, color = "black", margin = margin(b = 10)),
      strip.background = element_rect(fill = "gray95", color = NA), 
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray30")
    ) +
    labs(
      title = paste0("Shared Structure: ", obj_name, " (Top ", n_top, " Comp.)"),
      subtitle = "Scatterplot of weights with Cluster Parameters (Posterior Mean)",
      x = "Weight in Group 1",
      y = "Weight in Group 2"
    )
  
  print(p)
  return(p)
}

cluster_comparison <- function(res1, res2, label1 = "Model 1", label2 = "Model 2") {
  
  process_single_res <- function(res, label) {
    
    n_total <- dim(res$trace$Z)[1]
    idx_keep <- posterior_indices(n_total, 0.7)
    Z <- res$trace$Z[idx_keep, ]
    
    counts_list <- apply(Z, 1, function(z_row) {
      as.numeric(sort(table(z_row), decreasing = TRUE))
    })
    
    max_len <- max(sapply(counts_list, length))
    counts_matrix <- t(sapply(counts_list, function(x) c(x, rep(0, max_len - length(x)))))
    
    # 4. Media
    avg_size <- colMeans(counts_matrix)
    avg_size <- avg_size[avg_size >= 0.5]
    
    return(data.frame(
      Rank = 1:length(avg_size), 
      Size = avg_size,
      Model = label
    ))
  }
  
  df1 <- process_single_res(res1, label1)
  df2 <- process_single_res(res2, label2)

  df_final <- rbind(df1, df2)

  p <- ggplot(df_final, aes(x = log(Rank), y = log(Size), color = Model)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line(size = 1.2) +
    scale_color_manual(values = c("steelblue", "firebrick")) + 
    labs(
      title = "Confronto Rank-Size: DP vs PY",
      subtitle = "Distribuzione Log-Log delle dimensioni dei cluster",
      y = "Log(Dimensione Media)", 
      x = "Log(Rango)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = c(0.8, 0.8),
      legend.background = element_rect(fill = "white", color = NA)
    )
  
  print(p)
  
  return(invisible(df_final))
}

generate_mixture_data <- function(n_obs = 100, mix = "Mix1", seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  
  Y1 <- numeric(n_obs)
  Y2 <- numeric(n_obs)
  true_means_1 <- c()
  true_means_2 <- c()
  
  sample_comp <- function(probs) sample(1:length(probs), 1, prob = probs)
  
  for(t in 1:n_obs) {
    if(mix == "Mix1") {
      
      p <- c(1/3, 1/3, 1/3)
      mus <- c(-10, 0, 10)
      sds <- c(1, 1, 1)
      k1 <- sample_comp(p); k2 <- sample_comp(p)
      Y1[t] <- rnorm(1, mus[k1], sds[k1])
      Y2[t] <- rnorm(1, mus[k2], sds[k2])
      true_means_1 <- mus; true_means_2 <- mus
      
    } else if(mix == "Mix2") {
      # Unbalanced weights across groups
      p1 <- c(1/3, 1/3, 1/3)
      p2 <- c(1/6, 4/6, 1/6)
      mus <- c(-10, 0, 10)
      sds <- c(1, 1, 1)
      k1 <- sample_comp(p1); k2 <- sample_comp(p2)
      Y1[t] <- rnorm(1, mus[k1], sds[k1])
      Y2[t] <- rnorm(1, mus[k2], sds[k2])
      true_means_1 <- mus; true_means_2 <- mus
      
    } else if(mix == "Mix3") {
      # Different means and tighter variance
      p <- rep(0.25, 4)
      m1 <- c(0, 3, 2, 5)
      s1 <- sqrt(c(0.5, 0.25, 0.25, 0.5))
      m2 <- c(0, 3, -3, 7)
      s2 <- sqrt(c(0.5, 0.25, 0.25, 0.5))
      k1 <- sample_comp(p); k2 <- sample_comp(p)
      Y1[t] <- rnorm(1, m1[k1], s1[k1])
      Y2[t] <- rnorm(1, m2[k2], s2[k2])
      true_means_1 <- m1; true_means_2 <- m2
      
    } else if(mix == "Mix4") {
      # Many components
      p1 <- rep(1/7, 7)
      m1 <- c(-10, -5, -3, 0, 3, 5, 10)
      s1 <- sqrt(rep(0.4, 7))
      p2 <- rep(1/3, 3)
      m2 <- c(-10, 0, 10)
      s2 <- c(1, 1, 1)
      k1 <- sample_comp(p1); k2 <- sample_comp(p2)
      Y1[t] <- rnorm(1, m1[k1], s1[k1])
      Y2[t] <- rnorm(1, m2[k2], s2[k2])
      true_means_1 <- m1; true_means_2 <- m2
    }else if(mix == "Mix_LongTail") {
      # --- NUOVA AGGIUNTA: STRUTTURA "HEAVY TAIL" ---
      # 2 Cluster Giganti (80% massa) + 20 Cluster Piccoli (20% massa)
      
      # Parametri
      means_big <- c(-10, 10)
      means_small <- seq(-5, 5, length.out = 20)
      mus <- c(means_big, means_small)
      
      sds <- c(rep(1, 2), rep(0.5, 20)) # Piccoli cluster un po' più stretti
      
      # Pesi: 0.4 per i grandi, 0.01 per i piccoli
      prob_big <- 0.4
      prob_small <- 0.2 / 20 # = 0.01
      p <- c(rep(prob_big, 2), rep(prob_small, 20))
      
      # Generazione
      k1 <- sample_comp(p); k2 <- sample_comp(p)
      Y1[t] <- rnorm(1, mus[k1], sds[k1])
      Y2[t] <- rnorm(1, mus[k2], sds[k2])
      true_means_1 <- mus; true_means_2 <- mus
    }
  }
  
  return(list(
    Y1 = Y1, 
    Y2 = Y2, 
    true_means = list(m1 = unique(true_means_1), m2 = unique(true_means_2)),
    model = mix
  ))
}

# Dati simulati

RUN_FULL_ANALYSIS <- identical(tolower(Sys.getenv("RUN_FULL_ANALYSIS", "false")), "true")

if(RUN_FULL_ANALYSIS) {
data1 <- generate_mixture_data(n_obs=200, mix="Mix1", seed=123)
data2 <- generate_mixture_data(n_obs=200, mix="Mix2", seed=123)
data3 <- generate_mixture_data(n_obs=200, mix="Mix3", seed=123)
data4 <- generate_mixture_data(n_obs=200, mix="Mix4", seed=123)
data_longtail <- generate_mixture_data(n_obs = 500, mix = "Mix_LongTail", seed = 42)

res1 <- run_DPYP_ANOVA(
  Y1 = data1$Y1, 
  Y2 = data1$Y2,
  K_max = 15,
  n_iter = 10000,
  burnin_pct = 0.7,
  PYP = TRUE,
  seed = 42
)
p1_d <- plot_posterior_mean(res1)
plot_diagnostics(res1$trace)
p1_scat <- scatter(res1, 0.7, 3)

res2 <- run_DPYP_ANOVA(
  Y1 = data2$Y1, 
  Y2 = data2$Y2,
  K_max = 15,
  n_iter = 10000,
  burnin_pct = 0.7,
  PYP = TRUE,
  seed = 42
)
p2_d <- plot_posterior_mean(res2)
plot_diagnostics(res2$trace)
p2_scat <- scatter(res2, 0.7, 3)

res3 <- run_DPYP_ANOVA(
  Y1 = data3$Y1, 
  Y2 = data3$Y2,
  K_max = 15,
  n_iter = 10000,
  PYP = T,
  c_init = 0.3,
  burnin_pct = 0.7,
  seed = 42
)
p3_d <- plot_posterior_mean(res3)
plot_diagnostics(res3$trace)
p3_scat <- scatter(res3, 0.7, 4)

res4 <- run_DPYP_ANOVA(
  Y1 = data4$Y1, 
  Y2 = data4$Y2,
  K_max = 15,
  n_iter = 10000,
  burnin_pct = 0.7,
  PYP = T,
  c_init = 0.3,
  seed = 42
)
p4_d <- plot_posterior_mean(res4)
plot_diagnostics(res4$trace)
p4_scat <- scatter(res4, 0.7, 7)

Density <- grid.arrange( p1_d, p2_d, p3_d, p4_d,
  ncol = 2,
  nrow = 2,
  top = "Analisi Simulazioni: Densità Posteriori"
)

Corr1 <- grid.arrange( p1_scat, p2_scat,
                      ncol = 1,
                      nrow = 2,
                      top = "Analisi Simulazioni: Corr pesi"
)
Corr2 <- grid.arrange( p3_scat, p4_scat,
                       ncol = 2,
                       nrow = 1,
                       top = "Analisi Simulazioni: Corr pesi"
)

# Applicazione reale. Dati reali: https://www.kaggle.com/datasets/ealtman2019/credit-card-transactions/data

credit_file <- file.path(PROJECT_DIR, "IBM Credit Data", "User0_credit_card_transactions.csv")
if(file.exists(credit_file)) {
data <- data.frame(read_csv(credit_file, show_col_types = FALSE))
data$Amount <- as.numeric(gsub("[\\$,]", "", data$Amount))
data08 <- data %>%
  filter(Year == 2008)%>%
  pull(Amount)


data09 <- data %>%
  filter(Year == 2009)%>%
  pull(Amount)

if(length(data08) == 0 || length(data09) == 0) {
  stop("No observations found for one of the selected years: 2008 or 2009.")
}

n <- min(length(data08), length(data09))
data08 <- data08[1:n]
data09 <- data09[1:n]

credit_data <- list(Y1 = data08, Y2 = data09)
res_credit <- run_DPYP_ANOVA(
  Y1 = credit_data$Y1, 
  Y2 = credit_data$Y2,
  theta1_init = 3.0,
  theta2_init = 3.0,
  c_init = 0.9,
  n_iter = 10000,
  burnin_pct = 0.7,
  PYP = TRUE,
  seed = 42
)
plot_posterior_mean(res_credit)
plot_diagnostics(res_credit$trace)
scatter(res_credit, 0.7, 5)
} else {
  warning("Credit-card CSV not found at: ", credit_file, ". Skipping real-data analysis.")
}

# Cluster size

res_DPYP <- run_DPYP_ANOVA(
  Y1 = data_longtail$Y1, 
  Y2 = data_longtail$Y2,
  K_max = 30,
  n_iter = 10000,
  burnin_pct = 0.7,
  c_init = 0.6,
  PYP = TRUE,
  seed = 42
)

res_DDP <- run_DPYP_ANOVA(
  Y1 = data_longtail$Y1, 
  Y2 = data_longtail$Y2,
  K_max = 30,
  n_iter = 10000,
  burnin_pct = 0.7,
  PYP = FALSE,
  c_init = 0,
  seed = 42
)
cluster_comparison(res_DDP, res_DPYP, "DDP", "DPYP")
}
