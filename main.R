rm(list = ls())
library(mvtnorm)
library(survival)
library(reshape2)
library(ggplot2)
library(Rcpp)
source("functions_DESIGN.R")

my_seed <- 2026-03-16
cov_dim <- 3
center_num <- 10
beta_star <- c(-0.5, 0.5, 1)

my_corr <- matrix(0.3, cov_dim, cov_dim)
diag(my_corr) <- 1

shape_seq <- seq(3, 2, length.out = center_num)
scale_seq <- seq(5, 10, length.out = center_num)

target_rate <- 0.02   # or 0.05

v_seq <- sapply(1:center_num, function(k) {
  calibrate_v_one_site(target_rate = target_rate,
                       cov_dim = cov_dim,
                       beta_star = beta_star,
                       shape_k = shape_seq[k],
                       scale_k = scale_seq[k],
                       Sigma = my_corr,
                       prob_strata1 = 0.6,
                       n_mc = 200000)
})

v_seq <- round(v_seq, digits = 2)

# first, generate a synthetic dataset, where we introduce left truncation (so it is a three-column format)
synthetic_data <- data.generation(center_num = center_num, # number of sites 
                                  center_size = sample(1000:3000, center_num), # site size for the K sites
                                  cov_dim = cov_dim, # covariate dimension
                                  beta_star = beta_star, # regression parameter 
                                  rho = 0.3,
                                  Sigma = my_corr,
                                  prob_strata1 = 0.6,
                                  shape_seq = shape_seq,
                                  scale_seq = scale_seq,
                                  v_seq = v_seq,
                                  my_seed = my_seed)

# we have ten covariates: V1~V10, one site indicator called "site", and one center indicator called "center"
# each site includes two centers; e.g., centers 1 and 2 belong to site 1, centers 3 and 4 belong to site 2
# we do stratification on V10 and center
covariate_names = paste("V", 1:cov_dim, sep = "")
strata_names = c("site", "yw_strata")


start_time <- proc.time()
precision <- 1e-6

synthetic_data_prentice <- synthetic_data

synthetic_data_prentice[synthetic_data_prentice$subcohort == 0, "time_in"] <- synthetic_data_prentice[synthetic_data_prentice$subcohort == 0, "time_out"] - precision


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


pooled_cox <- survival::coxph(formula_cox, data = synthetic_data_prentice)
pool_time <- (proc.time() - start_time)["elapsed"]


strata_names <- c("yw_strata")
# for ease of subsequent analysis, we transform synthetic dataset into a list
site_names <- unique(synthetic_data$site)
# data_list is a list of length K, each element contains the data frame for a site
data_list <- as.list(rep(NA, length(site_names)))
for(i in 1:length(site_names)){
  data_list[[i]] <- synthetic_data[synthetic_data$site == site_names[i], ]
  print(c(i, mean(data_list[[i]]$censor_ind)))  
}

# we implement meta analysis to find one initial estimator
start_time <- proc.time()

meta_res <- meta_analysis(data_list, covariate_names = covariate_names, strata_names = strata_names)

initial_beta <- meta_res$IVW_est

meta_time <- (proc.time() - start_time)["elapsed"]

# implementation of DESIGN
# first, we will use the function weight_CC to calculate the weights for the case-cohort design
# the purpose here is to "pre-calculate" the weight before calculating the log-likelihood
# this would accelerate the subsequent calculation of log-likelihood
# the argument "covariate_names" is required to indicate the names of the covariates in the data frame 
# currently, we only provide Prentice weight; more options will be provided later
start_time <- proc.time()
pre_processing <- weight_CC(data_list, covariate_names = covariate_names, strata_names = strata_names)

grad_list <- as.list(rep(NA, pre_processing$K))
hess_list <- as.list(rep(NA, pre_processing$K))
meat_list <- as.list(rep(NA, pre_processing$K))
for(k in 1:pre_processing$K){
  grad_list[[k]] <- grad_plk(beta = initial_beta, 
                             X = pre_processing$covariate_list[[k]],
                             failure_position = pre_processing$failure_position[[k]],
                             failure_num = pre_processing$failure_num[[k]],
                             risk_sets = pre_processing$risk_sets[[k]])
  
  hess_list[[k]] <- hess_plk(beta = initial_beta, 
                             X = pre_processing$covariate_list[[k]],
                             failure_num = pre_processing$failure_num[[k]],
                             risk_sets = pre_processing$risk_sets[[k]])
  
  ### then, calculate meat in each site
  local_data <- data_list[[k]]
  

  local_data[local_data$subcohort == 0, "time_in"] <- local_data[local_data$subcohort == 0, "time_out"] - precision
  
  
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
  
  local_cox <- survival::coxph(formula_cox,
                               data = local_data,
                               init = initial_beta,
                               iter = 0)
  
  score_resid <- resid(local_cox, type = "score")  # n x p matrix
  
  # Compute meat
  meat_list[[k]] <- crossprod(score_resid)
}

lead_site <- 1
grad_list_without_lead <- grad_list[-lead_site]
hess_list_without_lead <- hess_list[-lead_site]


result <- optim(par = initial_beta, fn = DESIGN_fun, 
                control = list(fnscale = -1), method = "BFGS", 
                covariate_lead = pre_processing$covariate_list[[lead_site]],
                failure_position = pre_processing$failure_position[[lead_site]],
                failure_num = pre_processing$failure_num[[lead_site]],
                risk_sets = pre_processing$risk_sets[[lead_site]],
                grad_list = grad_list_without_lead,
                hess_list = hess_list_without_lead,
                initial_beta = initial_beta)
DESIGN_est <- result$par

### variance estimate for DESIGN
bread <- Reduce("+", hess_list)
DESIGN_var <- solve(bread) %*% Reduce("+", meat_list) %*% solve(bread)
DESIGN_time <- (proc.time() - start_time)["elapsed"]

pool_time
meta_time
DESIGN_time


# pooled estimates
pooled_cox$coefficients
# DESIGN estimates
DESIGN_est
# meta estimates
meta_res$IVW_est

# variance of pooled estimates
diag(pooled_cox$var) 
# variance of DESIGN estimates
diag(DESIGN_var)
# variance of meta estimates
diag(meta_res$IVW_var)
