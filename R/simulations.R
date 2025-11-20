library(dplyr)
library(stringr)
library(remotes)
library(devtools)
library(MASS)
library(sjstats) 
library(np) 
library(stats)
library(car) 
library(dplyr)
library(lqa) 
library(glmnet) 
library(future)
library(foreach)
library(parallel)
library(doParallel)




## simulation initiation ####


### S1, SoSt
yita <- c(0.3,1,0.3) 
coe_yx <- c(0.3,0.3,0.3,0.3,0,0) 
coe_ax <- c(0.3,0.3,0,0,1,1) 
coe_mx <- c(0.3,0.3,0,0,0,0) 


### S2, SoWt
yita <- c(0.3,1,0.3) 
coe_yx <- c(0.3,0.3,0.3,0.3,0,0) 
coe_ax <- c(0.15,0.15,0,0,1,1) 
coe_mx <- c(0.3,0.3,0,0,0,0) 


### S2, WoSt
yita <- c(0.3,1,0.3) 
coe_yx <- c(0.15,0.15,0.3,0.3,0,0)
coe_ax <- c(0.3,0.3,0,0,1,1)
coe_mx <- c(0.3,0.3,0,0,0,0) 



### treatment value range
a1_all_20 <- c(seq(-1,-0.1,0.1),seq(0.1,1,0.1))
a1_all_20_char <- as.character(a1_all_20)
a1_all_20_char[1] <- "-1.0"
a1_all_20_char[20] <- "1.0"


### sample size (n) and covariate dimension (p), taking (n,p)=(2000,100) as example

n<-2000
p<-100
n.simulations <- 100







## procedure 1: GOAL-based variable selection for treatment GPS ####


cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

GOAL_S1_res <- foreach(i = 1:n.simulations, .combine = rbind, .packages = c('MASS','lqa','sjstats','glmnet','parallel')) %dopar% {
  ## rho0
  set.seed(i)
  X <- matrix(data=NA,nrow = n,ncol = p)
  for (nc in c(1:p)){
    X[,nc] <- runif(n,min=-1,max=1)
  }
  var.list<-paste("X",1:p,sep="")
  colnames(X) <- var.list
  
  TrueVar<-X[,c(1,6,11,16,21,26)]
  MyData <- as.data.frame(X)
  MyData$Trt<-as.numeric(TrueVar%*%coe_ax + runif(n=n,min=-2,max=2))
  MyData$Med<-as.numeric(TrueVar%*%coe_mx + yita[1]*MyData$Trt + runif(n=n,min=-2,max=2))
  MyData$Y<-as.numeric(TrueVar%*%coe_yx + yita[2]*MyData$Trt + yita[3]*MyData$Med + runif(n=n,min=-2,max=2))
  
  GOAL_deviance_res <- GOAL_deviance(n=n,Data=MyData,var.list=var.list,covar=NULL,gamma_convergence = 2) 
  Svar_ax_final<- GOAL_deviance_res$GOAL_AX$Svar_ax
  Svar_amx_final<- GOAL_deviance_res$GOAL_AMX$Svar_amx
  lambda_ax<-as.numeric(GOAL_deviance_res$GOAL_AX$lambda_ax)
  lambda_amx<-as.numeric(GOAL_deviance_res$GOAL_AMX$lambda_amx)
  
  y <- as.numeric(MyData$Y)
  m <- as.numeric(MyData$Med)
  a <- as.numeric(MyData$Trt)
  X_ax <- as.matrix(X[,Svar_ax_final])
  X_amx <- as.matrix(X[,Svar_amx_final])
  
  list(
    y=y,
    m=m,
    a=a,
    X_ax=X_ax,
    X_amx=X_amx,
    Svar_ax_final=Svar_ax_final,
    Svar_amx_final=Svar_amx_final,
    lambda_ax=lambda_ax,
    lambda_amx=lambda_amx
  )
  
}
stopCluster(cl)

save(GOAL_S1_res,file = "S1_sost_GOAL_res.Rdata")






## procedure 2: GPS-weighted mediation effect estimation ####


y_df <- as.data.frame(do.call(cbind,GOAL_S1_res[,1]))
m_df <- as.data.frame(do.call(cbind,GOAL_S1_res[,2]))
a_df <- as.data.frame(do.call(cbind,GOAL_S1_res[,3]))
X_ax_list <- GOAL_S1_res[,4]
X_amx_list <- GOAL_S1_res[,5]


sim_res_GOAL_S1 <- list()
ntrim_df_GOAL_S1 <- data.frame(matrix(NA,n.simulations,1)) 


for (j in 1:length(a1_all_20)){
  cat(sprintf("########################## cycle starts: a1 = %s : %s ###\n", a1_all_20[j], Sys.time()))
  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)
  
  boot_results_a1 <- foreach(i=1:n.simulations,.combine = rbind,.packages = c('np','stats')) %dopar% {
    y <- as.numeric(y_df[,i])
    m <- as.numeric(m_df[,i])
    a <- as.numeric(a_df[,i])
    X_ax <- X_ax_list[[i]]
    X_amx <- X_amx_list[[i]]
    
    res_list <- medweightcont(y = y, a = a, m = m, x_ax = X_ax, x_amx = X_amx, a0 = 0, a1 = a1_all_20[j], ATET = FALSE, 
                              trim = 0.05, lognorm = FALSE, bw = 2.34*(n^(-0.25)), boot = 500, cluster = NULL) 
    boot_res_ai <- as.data.frame(res_list$boot_results)
    rownames(boot_res_ai) <- paste(rownames(boot_res_ai), a1_all_20_char[j], sep = "_")
    ntrim <- res_list$ntrimmed
    list(ntrim=ntrim,boot_res_ai=boot_res_ai)
    
  }
  
  stopCluster(cl)
  sim_res_GOAL_S1[[j]] <- do.call(rbind,boot_results_a1[,2]) 
  ntrim <- as.data.frame(do.call(rbind,boot_results_a1[,1]))
  colnames(ntrim) <- paste("ntrim",a1_all_20_char[j],sep = "_")
  ntrim_df_GOAL_S1 <- cbind(ntrim_df_GOAL_S1,ntrim)
  
  cat(sprintf("#################### a1 = %s cycle ends. \n", a1_all_20[j]))
  
}

save(sim_res_GOAL_S1,file = "S1_sost_sim_res_GOAL_wp.Rdata")
ntrim_df_GOAL_S1 <- ntrim_df_GOAL_S1[,-1]
write.csv(ntrim_df_AdaLASSO_S1,file = "S1_sost_sim_ntrim_GOAL_wp.csv")









