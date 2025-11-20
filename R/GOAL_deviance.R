remotes::install_github("cran/lqa")

install.packages("remotes")
install.packages("devtools")
install.packages("stats")
install.packages("MASS")
remotes::install_github("strengejacke/sjstats")
install.packages("car")
install.packages("glmnet")

library(remotes)
library(devtools)
library(lqa) 
library(stats)
library(MASS)
library(sjstats) 
library(car) 
library(glmnet) 




## Function: Adaptive LASSO cross-validation deviance evaluation

#' @param n Number of observation
#' @param folds A vector indicating fold assignment for each observation in cross-validation
#' @param nfold Number of folds to use in cross-validation
#' @param covar Names of known confounders
#' @param il Exponent index selected from a predefined lambda vector
#' @param ig Exponent index selected from a predefined gamma vector
#' @param Data Data (data.frame) that include all covariates, outcome, treatment and mediator
#' @param beta_ax  Unpenalized coefficients of the covariates in the ‘full’ outcome regression model, used in the treatment GPS model given specified covariates
#' @param beta_amx Unpenalized coefficients of the covariates in the ‘full’ outcome regression model, used in the treatment GPS model given specified covariates and the mediator
#' @param w.full.form_ax Formula of the treatment GPS model given specified covariates
#' @param w.full.form_amx Formula of the treatment GPS model given specified covariates and the mediator


adalasso_cv_deviance <- function(n, folds, covar, nfold, il, ig, Data, beta_ax, beta_amx, w.full.form_ax, w.full.form_amx ) {
  
  fold_results <- lapply(1:nfold, function(fold) { 
    train_idx <- which(folds != fold)
    test_idx <- which(folds == fold)
    train_data <- Data[train_idx, ]
    test_data <- Data[test_idx, ]
    
    oal_pen_ax <- adaptive.lasso(lambda = n^(il), al.weights = c(rep(0, length(covar)), abs(beta_ax)^(-1*ig)))
    alasso_ax <- lqa.formula(formula = w.full.form_ax, data = train_data, penalty = oal_pen_ax, family = gaussian())
    pred_ax <- as.matrix(test_data[, names(coef(alasso_ax)[-1]), drop = FALSE]) %*% coef(alasso_ax)[-1] + coef(alasso_ax)[1]
    deviance_ax <- sum((test_data$Trt - pred_ax)^2)
    
    oal_pen_amx <- adaptive.lasso(lambda = n^(il), al.weights = c(0, rep(0, length(covar)), abs(beta_amx)^(-1*ig)))
    alasso_amx <- lqa.formula(formula = w.full.form_amx, data = train_data, penalty = oal_pen_amx, family = gaussian())
    pred_amx <- as.matrix(test_data[, names(coef(alasso_amx)[-1]), drop = FALSE]) %*% coef(alasso_amx)[-1] + coef(alasso_amx)[1]
    deviance_amx <- sum((test_data$Trt - pred_amx)^2)
    
    list(deviance_ax = deviance_ax, deviance_amx = deviance_amx)
    
  }) 
  
  mean_deviance_ax <- mean(sapply(fold_results, function(res) res$deviance_ax))
  mean_deviance_amx <- mean(sapply(fold_results, function(res) res$deviance_amx))
  
  list(mean_deviance_ax = mean_deviance_ax, mean_deviance_amx = mean_deviance_amx)
  
}











## Function: Generalized outcome-adaptive LASSO (GOAL) estimation based on cross-validation

#' @param n Number of observation
#' @param Data Data (data.frame) that include all covariates, outcome, treatment and mediator.
#' @param var.list Names of potential unknown confounders
#' @param covar Names of known confounders
#' @param lambda_vec  A vector of possible lambda values and default is c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
#' @param gamma_convergence The value of gamma convergence and default is 2



GOAL_deviance <- function(n, Data, var.list, covar = NULL, lambda_vec = c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49), gamma_convergence = 2) {
  
  if (is.null(covar)) {
    allvar <- var.list
  } else {
    allvar <- c(covar, var.list)
  }
  
  Data <- Data[,c(allvar,"Trt","Med","Y")]
  temp.mean = colMeans(Data[,c(allvar)])  #每个X计算均值（列均值）
  Temp.mean = matrix(temp.mean,ncol=length(allvar),nrow=nrow(Data),byrow=TRUE)  #构建一个以X均值组成的N*p的矩阵
  Data[,allvar] = Data[,allvar] - Temp.mean #构建X观察值和均值的差值
  temp.sd = apply(Data[allvar], MARGIN=2, FUN=sd) #margin=2:列 → 按列计算sd
  Temp.sd = matrix(temp.sd,ncol=length(allvar),nrow=nrow(Data),byrow=TRUE) #构建一个以sd(X)组成的N*p的矩阵
  Data[allvar] = Data[,allvar] / Temp.sd #Data数据框中的X值为mean/sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd")) #从系统environment中移除这些data
  
  y.form <- formula(paste("Y","~",paste(c("Trt","Med",allvar),collapse="+")))
  table_y <- as.data.frame(table(Data$Y))
  if(nrow(table_y)==2){
    lm.Y <- glm(y.form,data=Data,family = binomial(link = "logit"))  
  }else{
    lm.Y <- lm(y.form,data=Data)
  }
  betaXY <- coef(lm.Y)[var.list] #提取未知混杂的系数
  
  names(lambda_vec) <- as.character(lambda_vec)
  gamma_vals <- 2*( gamma_convergence - lambda_vec + 1 ) #此处lambda_vec=log(lambda,base=n)！！已经是幂了！
  names(gamma_vals) <- names(lambda_vec)
  
  w.full.form_ax <- formula(paste("Trt","~",paste(allvar,collapse="+")))  #A~X的模型
  w.full.form_amx <- formula(paste("Trt","~",paste(c("Med",allvar),collapse="+")))   #拟合模型：A~M,X
  
  set.seed(123)  # 保证结果可重复
  nfold=10
  folds <- sample(rep(1:nfold, length.out = nrow(Data)))
  
  results <- mapply(function(il, ig) adalasso_cv_deviance(n=n, folds=folds, covar=covar, nfold=nfold, il, ig, Data=Data, beta_ax=betaXY, beta_amx=betaXY, w.full.form_ax=w.full.form_ax, w.full.form_amx=w.full.form_amx), lambda_vec, gamma_vals, SIMPLIFY = FALSE)
  
  cv_results_ax <- data.frame(lambda = lambda_vec, gamma = gamma_vals, deviance = sapply(results, function(res) res$mean_deviance_ax))
  cv_results_amx <- data.frame(lambda = lambda_vec, gamma = gamma_vals, deviance = sapply(results, function(res) res$mean_deviance_amx))
  
  sorted_idx_ax <- order(cv_results_ax$deviance)
  for (i in sorted_idx_ax) {
    lambda_ax <- cv_results_ax$lambda[i]
    gamma_ax <- cv_results_ax$gamma[i]
    oal_pen_ax <- adaptive.lasso(lambda = n^(lambda_ax), al.weights = c(rep(0, length(covar)), abs(betaXY)^(-1*gamma_ax)))
    alasso_ax <- lqa.formula(formula = w.full.form_ax, data = Data, penalty = oal_pen_ax, family = gaussian())
    GOAL_ax_coef <- coef(alasso_ax)[var.list]
    Svar_ax <- c(covar, names(GOAL_ax_coef)[which(round(GOAL_ax_coef, 5) != 0)])
    if (length(Svar_ax) > 0) {break}
  }
  
  sorted_idx_amx <- order(cv_results_amx$deviance)
  for (j in sorted_idx_amx) {
    lambda_amx <- cv_results_amx$lambda[j]
    gamma_amx <- cv_results_amx$gamma[j]
    oal_pen_amx <- adaptive.lasso(lambda = n^(lambda_amx), al.weights = c(0, rep(0, length(covar)), abs(betaXY)^(-1*gamma_amx)))
    alasso_amx <- lqa.formula(formula = w.full.form_amx, data = Data, penalty = oal_pen_amx, family = gaussian())
    GOAL_amx_coef <- coef(alasso_amx)[var.list]
    Svar_amx <- c(covar, names(GOAL_amx_coef)[which(round(GOAL_amx_coef, 5) != 0)])
    if (length(Svar_amx)>0) {break}
  }
  
  GOAL_results<-list(   
    GOAL_AX=list(
      lambda_ax=lambda_ax,
      Svar_ax=Svar_ax
    ),
    GOAL_AMX=list(
      lambda_amx=lambda_amx,
      Svar_amx=Svar_amx
    )
  )
  
  return (GOAL_results)
  
}










