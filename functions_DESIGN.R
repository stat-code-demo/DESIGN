get_event_rate <- function(v,
                           n_mc = 200000,
                           cov_dim,
                           beta_star,
                           shape_k,
                           scale_k,
                           target_prob_strata1 = 0.6,
                           Sigma) {
  set.seed(2026-03-16)
  
  Z <- rmvnorm(n = n_mc,
               mean = rep(0, cov_dim),
               sigma = Sigma)
  
  strata_index <- rbinom(n = n_mc, size = 1, prob = target_prob_strata1)
  
  shape_i <- ifelse(strata_index == 1, shape_k, 1)
  
  T0 <- numeric(n_mc)
  id1 <- which(strata_index == 1)
  id0 <- which(strata_index == 0)
  
  T0[id1] <- rweibull(length(id1), shape = shape_k, scale = scale_k)
  T0[id0] <- rweibull(length(id0), shape = 1, scale = scale_k)
  
  eta <- as.vector(Z %*% beta_star)
  
  T <- T0 * exp(- eta / shape_i)
  
  C <- runif(n_mc, min = 0, max = v)
  
  mean(T <= C)
}

calibrate_v_one_site <- function(target_rate,
                                 cov_dim,
                                 beta_star,
                                 shape_k,
                                 scale_k,
                                 Sigma,
                                 prob_strata1 = 0.6,
                                 n_mc = 200000,
                                 lower = 1e-4,
                                 upper = 50) {
  
  f <- function(v) {
    get_event_rate(v = v,
                   n_mc = n_mc,
                   cov_dim = cov_dim,
                   beta_star = beta_star,
                   shape_k = shape_k,
                   scale_k = scale_k,
                   target_prob_strata1 = prob_strata1,
                   Sigma = Sigma) - target_rate
  }
  
  # expand upper bound if needed
  while (f(upper) < 0) {
    upper <- upper * 2
    if (upper > 1e5) stop("Cannot find a sufficiently large upper bound.")
  }
  
  uniroot(f, lower = lower, upper = upper)$root
}


# this function is used for data generation
# the generated data set is a dataframe, it has the following columns
# column 1 (ID): patient ID
# column 2 (site): site ID
# column 3 (subcohort): a binary indicator showing whether the patient is included in the subcohort
# column 4 (t2e): the time-to-event outcome, i.e., the min of survival time and censoring time
# column 5 (censor_ind): the censoring indicator
# from column 6 to the second last column: covariates
# last_column (site_size): the full cohort size
data.generation <- function(center_num, center_size, cov_dim, beta_star, rho, Sigma,
                            prob_strata1 = 0.6,
                            shape_seq,
                            scale_seq,
                            v_seq,
                            my_seed){
  # K is number of sites 
  # center_size is site size for the K sites
  # cov_dim is covariate dimension
  # beta_star is regression parameter 
  # event_rate controls the event rate
  # rho is the proportion for the subcohort
  set.seed(my_seed)
  effective_center_size <- rep(NA, center_num)
  data_list <- as.list(rep(NA, center_num))

  for(k in 1:center_num){
    covariate_mat <- rmvnorm(n = center_size[k],
                             mean = rep(0, cov_dim),
                             sigma = Sigma)
    strata_index <- rbinom(n = center_size[k], size = 1, prob = prob_strata1)
    
    hazard_scaling <- c(covariate_mat %*% beta_star)
    
    baseline_t2e <- rep(NA, center_size[k])
    Cox_t2e <- rep(NA, center_size[k])
    for(yw in 1:center_size[k]){
      if(strata_index[yw] == 1){
        temp_shape <- shape_seq[k]
      }else{
        temp_shape <- 1
      }
      baseline_t2e[yw] <- rweibull(1, shape = temp_shape, scale = scale_seq[k])
      Cox_t2e[yw] <- baseline_t2e[yw] * exp(-hazard_scaling[yw] / temp_shape)
    }
    
    
    censor <- runif(center_size[k], min = 0, max = v_seq[k])
    trucation <- rep(0, center_size[k]) 
    
    un_truncated <- which(Cox_t2e > trucation)
    
    Cox_t2e <- Cox_t2e[un_truncated]
    censor <- censor[un_truncated]
    trucation <- trucation[un_truncated]
    covariate_mat <- covariate_mat[un_truncated, ]
    strata_index <- strata_index[un_truncated]
    effective_center_size[k] <- length(un_truncated)
    
    if(k == 1){
      local_data <- cbind(c(1:effective_center_size[k]),
                          rep(k, effective_center_size[k]), 
                          rep(NA, effective_center_size[k]),
                          trucation,
                          pmin(Cox_t2e, censor), 
                          as.numeric(Cox_t2e <= censor), 
                          covariate_mat,
                          strata_index)
    }else{
      local_data <- cbind(c(1:effective_center_size[k]) + sum(effective_center_size[1:(k-1)]),
                          rep(k, effective_center_size[k]), 
                          rep(NA, effective_center_size[k]),
                          trucation,
                          pmin(Cox_t2e, censor), 
                          as.numeric(Cox_t2e <= censor), 
                          covariate_mat,
                          strata_index)
    }
    
    
    local_data <- as.data.frame(local_data)
    names(local_data) <- c("ID", "site", "subcohort", "time_in", "time_out", "censor_ind", paste("V", 1:cov_dim, sep = ""), "yw_strata")
    local_data$subcohort <- rep(0, effective_center_size[k])
    local_data$subcohort[sample(1:effective_center_size[k], size = floor(rho * effective_center_size[k]))] <- 1 
    local_data <- local_data[which((local_data$subcohort == 1)|(local_data$censor_ind == 1)), ]
    data_list[[k]] <- local_data
  }
  return(do.call(rbind, data_list))
}

# the purpose of this function is to "pre-calculate" the weight before calculating the log-likelihood
# this would accelerate the subsequent calculation of log-likelihood
# currently, we only provide Prentice weight; more options will be provided later
weight_CC <- function(data_list, covariate_names, strata_names = c()){
  # for each site, pre-calculate the failure time points, the risk sets, and the respective weights
  # also, remove those sites with zero events
  K <- length(data_list)
  failure_num <- rep(NA, K)
  failure_times <- as.list(rep(NA, K))
  risk_sets <- as.list(rep(NA, K))
  covariate_list <- as.list(rep(NA, K))
  failure_position <- as.list(rep(NA, K))
  for(k in 1:K){
    local_data <- data_list[[k]]
    # prepare a list for covariates in matrix format so as to speed up computation of log partial likelihood, gradient, and hessian
    covariate_list[[k]] <- as.matrix(local_data[, covariate_names, drop = FALSE]) 
    # find over which position lies the failure times
    failure_position[[k]] <- which(local_data$censor_ind == 1)
    # find failure times
    failure_times[[k]] <- local_data$time_out[failure_position[[k]]]
    # the number of failures
    failure_num[k] <- length(failure_times[[k]])
    
    local_data[local_data$subcohort == 0, "time_in"] <- local_data[local_data$subcohort == 0, "time_out"] - 1e-6
    
    if (length(strata_names) == 0L) {
      strata_id <- rep.int(1L, nrow(local_data))
    } else {
      strata_vars <- local_data[, strata_names, drop = FALSE]
      if (ncol(strata_vars) == 1L) {
        strata_id <- as.integer(factor(strata_vars[[1]]))
      } else {
        strata_id <- as.integer(interaction(strata_vars, drop = TRUE, lex.order = TRUE))
      }
    }
    
    temp_risk <- as.list(rep(NA, failure_num[k]))
    for(j in 1:failure_num[k]){
      fpos_j <- failure_position[[k]][j]
      t_j    <- failure_times[[k]][j]
      
      idx_strata <- strata_id == strata_id[fpos_j]
      idx_time   <- (local_data$time_in <= t_j) & (local_data$time_out >= t_j)
      
      temp_risk[[j]] <- which(idx_strata & idx_time)
    }
    risk_sets[[k]] <- temp_risk
  }

  return(list(data_list = data_list,
              covariate_list = covariate_list,
              failure_position = failure_position,
              failure_num = failure_num,
              risk_sets = risk_sets,
              K = K))
}

### *******************************************************************************************************************
### *******************************************************************************************************************
### *******************************************************************************************************************
# meta analysis
meta_analysis <- function(data_list, covariate_names, strata_names){
  K <- length(data_list)
  cov_dim <- length(covariate_names)
  term1 <- rep(0, cov_dim)
  term2 <- matrix(0, cov_dim, cov_dim)
  
  for(k in 1:K){
    local_data <- data_list[[k]]
    local_data[local_data$subcohort == 0, "time_in"] <- local_data[local_data$subcohort == 0, "time_out"] - 1e-6
    
    if(is.null(strata_names)){
      formula_cox <- as.formula(paste0("Surv(time_in, time_out, censor_ind) ~",
                                       paste0(covariate_names, collapse = " + "), 
                                       "+ cluster(ID)"))
    }else{
      formula_cox <- as.formula(paste0("Surv(time_in, time_out, censor_ind) ~",
                                       paste0(covariate_names, collapse = " + "), 
                                       " + strata(",
                                       paste0(strata_names, collapse = " , "),
                                       ") ",
                                       "+ cluster(ID)"))
    }
    
    
    local_fit <- survival::coxph(formula_cox, data = local_data)
    inv_var <- solve(local_fit$var) 
    term1<- term1 + inv_var %*% local_fit$coefficients
    term2 <- term2 + inv_var
  }
  
  IVW_est <- c(solve(term2) %*% term1)
  
  return(list(IVW_est = IVW_est, IVW_var = solve(term2)))
}
### *******************************************************************************************************************
### *******************************************************************************************************************
### *******************************************************************************************************************



# =======================
# Optimized R Implementation (No Rcpp)
# =======================

log_plk <- function(beta, covariate, failure_position, failure_num, risk_sets) {
  eta <- covariate %*% beta
  exp_eta <- exp(eta)
  res <- sum(eta[failure_position])
  
  for (j in 1:failure_num) {
    idx <- risk_sets[[j]]
    res <- res - log(sum(exp_eta[idx]))
  }
  return(res)
}


grad_plk <- function(beta, X, failure_position, failure_num, risk_sets) {
  eta <- X %*% beta
  exp_eta <- exp(eta)
  
  grad <- colSums(X[failure_position, , drop = FALSE])
  
  for (j in 1:failure_num) {
    idx <- risk_sets[[j]]
    temp_w <- exp_eta[idx] 
    denom <- sum(temp_w)
    weighted_X <- sweep(X[idx, , drop = FALSE], 1, temp_w, '*')
    grad <- grad - colSums(weighted_X) / denom
  }
  return(grad)
}

hess_plk <- function(beta, X, failure_num, risk_sets) {
  eta <- X %*% beta
  exp_eta <- exp(eta)
  d <- ncol(X)
  H <- matrix(0, d, d)
  
  for (j in 1:failure_num) {
    idx <- risk_sets[[j]]
    temp_w <- exp_eta[idx]
    denom <- sum(temp_w)
    
    X_sub <- X[idx, , drop = FALSE]
    weighted_X <- sweep(X_sub, 1, temp_w, '*')
    mean_vec <- colSums(weighted_X)
    
    sqrt_wX <- sweep(X_sub, 1, sqrt(temp_w), '*')
    
    H <- H + (tcrossprod(mean_vec) / (denom^2)) - (crossprod(sqrt_wX) / denom)
  }
  return(H)
}

DESIGN_fun <- function(beta, covariate_lead, failure_position, failure_num, risk_sets, grad_list, hess_list, initial_beta) {
  surrogate <- sum(Reduce(`+`, grad_list) * beta) + 0.5 * crossprod(beta - initial_beta, Reduce(`+`, hess_list) %*% (beta - initial_beta))
  
  surrogate + log_plk(beta, covariate_lead, failure_position, failure_num, risk_sets)
}
