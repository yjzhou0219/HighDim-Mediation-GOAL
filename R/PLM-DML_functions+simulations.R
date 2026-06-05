# ==============================================================================
# Simulation Code for Comparison Method: PLM-DML (Yang et al., 2025)
# Framework: Partially Linear Mediation Models (PLM) + Double Machine Learning (DML)
#
# Model Specification：
#   M = alpha*A + g(X) + e1    --> 'alpha' is estimated via DoubleMLPLR
#   Y = theta*A + beta*M + f(X) + e2 --> 'theta' and 'beta' are estimated in stages
#
# Target Causal Parameters:
#   theta : Natural Direct Effect (NDE)
#   alpha : Effect of Exposure (A) on Mediator (M)
#   beta  : Effect of Mediator (M) on Outcome (Y)
#   NIE   : Natural Indirect Effect, calculated as (alpha * beta). 
#           Standard Error (SE) is derived via the Delta method.
#   TE    : Total Effect, calculated as (theta + alpha * beta).
#
# Estimand Note:
#   Unlike the proposed method which evaluates finite-difference counterfactuals dependent on specific 'a1' values, the PLM-DML framework estimates global linear coefficients. 
#   Thus, its results align with our framework under the scenario where (a1 - a0) = 1 (i.e., a unit change in exposure).
# ==============================================================================



# install.packages("DoubleML")
# install.packages("mlr3")
# install.packages("mlr3learners")
# install.packages("data.table")

library(DoubleML)
library(mlr3)
library(mlr3learners)  
library(data.table)
library(doParallel)
library(foreach)
library(parallel)




# ---- 1. Simulation Parameter Setup ----


n <- 5000
p <- 200

## S1

yita <- c(0.3,1,0.3) # Coefficients for M~A and Y~A+M
coe_yx <- c(0.3,0.3,0.3,0.3,0,0) # X1, X6 (confounders), X11, X16 (predictors) relate to Y
coe_ax <- c(0.3,0.3,0,0,1,1) # X1, X6 (confounders), X21, X26 (IVs) relate to A
coe_mx <- c(0.3,0.3,0,0,0,0) # X1, X6 (confounders) relate to M


## S2

yita <- c(0.3,1,0.3) 
coe_yx <- c(0.3,0.3,0.3,0.3,0,0) 
coe_ax <- c(0.15,0.15,0,0,1,1) 
coe_mx <- c(0.3,0.3,0,0,0,0) 


## S3

yita <- c(0.3,1,0.3) 
coe_yx <- c(0.15,0.15,0.3,0.3,0,0) 
coe_ax <- c(0.3,0.3,0,0,1,1) 
coe_mx <- c(0.3,0.3,0,0,0,0) 


n.simulations <- 100



# True parameter values for PLM-DML (Global linear coefficients)

true_vals_PLM <- c(
  theta = yita[2],                    # NDE = 1
  alpha = yita[1],                    # A->M = 0.3
  beta  = yita[3],                    # M->Y = 0.3
  NIE   = yita[1] * yita[3],          # NIE = 0.09
  TE    = yita[2] + yita[1]*yita[3]   # TE = 1.09
)





# ---- 2. Data-Generating Process (DGP) ----


sim_data_list <- vector("list", n.simulations)
for (i in 1:n.simulations) {
  set.seed(i)
  
  # Generate independent uniform covariates, taking rho=0 as an example
  X_mat   <- matrix(runif(n * p, -1, 1), nrow = n, ncol = p)
  
  colnames(X_mat) <- paste0("X", 1:p)
  TrueVar <- X_mat[, c(1, 6, 11, 16, 21, 26)]
  
  a  <- as.numeric(TrueVar %*% coe_ax + runif(n, -2, 2))
  m  <- as.numeric(TrueVar %*% coe_mx + yita[1] * a + runif(n, -2, 2))
  y  <- as.numeric(TrueVar %*% coe_yx + yita[2] * a + yita[3] * m + runif(n, -2, 2))
  
  sim_data_list[[i]] <- list(Y = y, A = a, M = m, X = X_mat)
}

save(sim_data_list,file = "sim_data_list.Rdata")






# ---- 3. Core Estimation Function: PLM-DML ----


plm_dml_mediation <- function(Y, A, M, X_mat, n_folds = 5) {
  
  
  # --- Step 1: Estimate alpha (Effect of A on M) ---
  # Model: M = alpha*A + g(X) --> PLR Model (outcome=M, treatment=A, controls=X)
  
  data_M <- data.table(M = M, A = A, X_mat)
  dml_data_M <- DoubleMLData$new(
    data = data_M, y_col = "M", d_cols = "A", x_cols = paste0("X", 1:p)
  )
  
  # Nuisance parameter estimators: Cross-validated LASSO
  learner_l <- lrn("regr.cv_glmnet", nfolds = 5, alpha = 1)  # LASSO for E[M|X]
  learner_m <- lrn("regr.cv_glmnet", nfolds = 5, alpha = 1)  # LASSO for E[A|X]
  
  # Cross-fitting and Partialling out: OLS on residuals
  plr_M <- DoubleMLPLR$new(
    data  = dml_data_M, ml_l = learner_l$clone(), ml_m = learner_m$clone(), # ml_l: 预测结局变量; ml_m: 预测处理变量
    n_folds = n_folds, score = "partialling out" # partialling out: 用拟合的M对A做OLS
  )
  plr_M$fit(store_predictions = FALSE) 
  
  alpha_hat <- plr_M$coef["A"]
  alpha_se  <- plr_M$se["A"]
  
  
  # --- Step 2: Estimate theta (NDE) and beta (Effect of M on Y) ---
  # Model: Y = theta*A + beta*M + f(X) --> PLR Model (outcome=Y, treatment=c(A,M), controls=X)
  
  data_Y <- data.table(Y = Y, A = A, M = M, X_mat)
  dml_data_Y <- DoubleMLData$new(
    data = data_Y, y_col = "Y", d_cols = c("A", "M"), x_cols = paste0("X", 1:p) 
    # d_cols中同时估计theta和beta
  )
  
  # Note: DoubleML automatically fits independent nuisance models for multiple treatment variables
  plr_Y <- DoubleMLPLR$new(
    data = dml_data_Y,
    ml_l = lrn("regr.cv_glmnet", nfolds = 5, alpha = 1)$clone(), # E[Y|X]
    ml_m = lrn("regr.cv_glmnet", nfolds = 5, alpha = 1)$clone(), # E[M|X] & E[A|X]
    n_folds = n_folds,
    score = "partialling out" 
  )
  plr_Y$fit(store_predictions = FALSE)
  
  theta_hat <- plr_Y$coef["A"]   # NDE
  beta_hat  <- plr_Y$coef["M"]   # Effect of M on Y
  theta_se  <- plr_Y$se["A"]
  beta_se   <- plr_Y$se["M"]
  
  
  # --- Step 3: Calculate NIE and TE, deriving SEs via Delta Method ---
  # NIE = alpha * beta (Product method for linear mediation)
  # Var(NIE) ≈ beta^2 * Var(alpha) + alpha^2 * Var(beta) (Taylor expansion approximation)
  
  NIE_hat <- alpha_hat * beta_hat
  NIE_se  <- sqrt(beta_hat^2 * alpha_se^2 + alpha_hat^2 * beta_se^2)
  
  # TE = NDE + NIE = theta + alpha*belta
  # Var(TE) ≈ Var(theta) + Var(NIE) (Assuming independent estimation steps)
  TE_hat <- theta_hat + NIE_hat
  TE_se  <- sqrt(theta_se^2 + NIE_se^2)
  

  list(
    coef = c(TE = TE_hat, NDE1 = theta_hat, NDE0 = theta_hat, NIE1 = NIE_hat, NIE0 = NIE_hat),
    se = c(TE = TE_se, NDE1 = theta_se, NDE0 = theta_se, NIE1 = NIE_se, NIE0 = NIE_se)
  )
  
}






# ---- 4. Parallel Monte Carlo Simulation Loop ----


load("sim_data_list.Rdata")

cat(sprintf("=== PLM-DML Simulation Started [%s] ===\n", Sys.time()))

n_cores <- detectCores() - 3
cl <- makeCluster(n_cores)
registerDoParallel(cl)

coef_results_PLM <- foreach(i = 1:n.simulations, .combine = rbind, .packages = c("DoubleML", "mlr3", "mlr3learners", "data.table"), 
                            .errorhandling = "pass") %dopar% {

  set.seed(i)
  
  dat_i  <- sim_data_list[[i]]
  Y      <- dat_i$Y
  A      <- dat_i$A
  M      <- dat_i$M
  X_mat  <- dat_i$X
  
  # Execute PLM-DML estimation
  res <- tryCatch(
    plm_dml_mediation(Y, A, M, X_mat, n_folds = 5),
    error = function(e) NULL
  )
  
  if(is.null(res)) {
    res_list <- list(
      estimates  = c(TE=NA, NDE1=NA, NDE0=NA, NIE1=NA, NIE0=NA),
      se = c(TE=NA, NDE1=NA, NDE0=NA, NIE1=NA, NIE0=NA),
      failed = TRUE
    )
  } else {
    res_list <- list(estimates = res$coef, se = res$se, failed = FALSE)
  }
  
  return(res_list)
  
}

stopCluster(cl)

cat(sprintf("=== PLM-DML Simulation Completed [%s] ===\n", Sys.time()))

save(coef_results_PLM,file = "coef_results_PLM.Rdata")





