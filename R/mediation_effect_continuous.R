#library(devtools)
#devtools::install_github("kaskarn/causamed")

#setwd("/Users/mymymy/Desktop/biostatistics_fdu/3 统计学课题/1 proj_1/code/")

install.packages("remotes")
install.packages("devtools")
install.packages("MASS")
remotes::install_github("strengejacke/sjstats")
remotes::install_github("JeffreyRacine/R-Package-np")
install.packages("stats")
install.packages("car")
install.packages("dplyr")

library(remotes)
library(devtools)
library(MASS)
library(sjstats) 
library(np) 
library(stats)
library(car) 
library(dplyr)








## Function 1. Mediation effects estimation for continuous treatment

#' @param y Dependent variable, must not contain missings.
#' @param a Continuous treatment, must not contain missings.
#' @param m Continuous mediator, may be a scalar or a vector, must not contain missings.
#' @param x_ax Data matrix of pre-treatment covariates selected in the treatment GPS model given specified covariates.
#' @param x_amx Data matrix of pre-treatment covariates selected in the treatment GPS model given specified covariates and the mediator.
#' @param a0 Value of \code{a} under non-treatment. Effects are based on pairwise comparisons, i.e. differences in potential outcomes evaluated at \code{a1} and \code{a0}.
#' @param a1 Value of \code{a} under treatment. Effects are based on pairwise comparisons, i.e. differences in potential outcomes evaluated at \code{a1} and \code{a0}.
#' @param ATET If FALSE, the average treatment effect (ATE) and the corresponding direct and indirect effects are estimated. If TRUE, the average treatment effect on the treated (ATET)  and the corresponding direct and indirect effects are estimated. Default is FALSE.
#' @param trim Trimming rule for discarding observations with extreme generalized propensity scores. If \code{lognorm=FALSE}, observations with f(D=\code{a1}|M,X)<\code{trim} or f(D=\code{a0}|M,X)<\code{trim} are dropped, with f denoting the generalized propensity score (or conditional density of treatment). If \code{lognorm=TRUE}, then \code{trim} corresponds to the share of lowest f(D=\code{a1}|M,X) or f(D=\code{a0}|M,X), respectively, that are dropped.
#' @param lognorm If FALSE, a linear model with normally distributed errors is assumed for generalized propensity score estimation. If TRUE, a lognormal model is assumed. Default is FALSE.
#' @param bw Bandwith for the second order Epanechnikov kernel functions of the treatment. If set to NULL, the rule of thumb for Epanechnikov kernels is used for bandwidth computation. Default is NULL.
#' @details Estimation of causal mechanisms (natural direct and indirect effects) of a continuous treatment under a selection on observables assumption assuming that all confounders of the treatment and the mediator, the treatment and the outcome, or the mediator and the outcome are observed. Units are weighted by the inverse of their conditional treatment densities (known as generalized propensity scores) given the mediator and/or observed confounders, which are estimated by linear or loglinear regression.



mediation.cont<-function(y,a,m,x_ax,x_amx,a0,a1,ATET=FALSE,trim=0.05,lognorm=FALSE,bw){ #本研究中此处输入的x已经是经GMAL方法选择后的协变量子集V
  
  if(lognorm==TRUE){
    aa=a;aa[a==0]=0.00001  
    ggg=glm(log(aa)~x_ax) 
    if (a0==0) a0=0.00001; if (a1==0) a1=0.00001
    pscore1a0=(dnorm( (log(a0)-cbind(1,x_ax)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2)))/a0)
    pscore1a1=(dnorm( (log(a1)-cbind(1,x_ax)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2)))/a1)
    ggg=glm(log(aa)~cbind(x_amx,m))
    pscore2a0=(dnorm( (log(a0)-cbind(1,x_amx,m)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2)))/a0)
    pscore2a1=(dnorm( (log(a1)-cbind(1,x_amx,m)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2)))/a1)
  }
  
  if(lognorm==FALSE){
    ggg=glm(a~x_ax)
    pscore1a0=(dnorm( (a0-cbind(1,x_ax)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2))))
    pscore1a1=(dnorm( (a1-cbind(1,x_ax)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2))))
    ggg=glm(a~cbind(x_amx,m))
    pscore2a0=(dnorm( (a0-cbind(1,x_amx,m)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2))))
    pscore2a1=(dnorm( (a1-cbind(1,x_amx,m)%*%ggg$coefficients)/sqrt(mean(ggg$residuals^2))))
  }
  
  kernwgta0=npksum(bws=bw, txdat = a, tydat = y, exdat = a0, return.kernel.weights=TRUE, ckertype="epanechnikov", ckerorder=2)$kw
  kernwgta1=npksum(bws=bw, txdat = a, tydat = y, exdat = a1, return.kernel.weights=TRUE, ckertype="epanechnikov", ckerorder=2)$kw
  
  if(lognorm==FALSE) ind=( (pscore2a0>=trim) & (pscore2a1>=trim) )  
  if(lognorm==TRUE)  ind=( (pscore2a0>=quantile(pscore2a0, trim)) & (pscore2a1>=quantile(pscore2a1, trim)) )
  y=y[ind];  pscore1a0=pscore1a0[ind];pscore1a1=pscore1a1[ind];pscore2a0=pscore2a0[ind]; pscore2a1=pscore2a1[ind];
  kernwgta0=kernwgta0[ind]; kernwgta1= kernwgta1[ind];
  
  if (ATET==FALSE){
    ya1m1=sum(y*kernwgta1/pscore1a1)/sum(kernwgta1/pscore1a1)
    ya1m0=sum(y*kernwgta1*pscore2a0/(pscore2a1*pscore1a0))/sum(kernwgta1*pscore2a0/(pscore2a1*pscore1a0))
    ya0m0=sum(y*kernwgta0/pscore1a0)/sum(kernwgta0/pscore1a0)
    ya0m1=sum(y*kernwgta0*pscore2a1/(pscore2a0*pscore1a1))/sum(kernwgta0*pscore2a1/(pscore2a0*pscore1a1))
    relative_wgt_a1m1 <- as.data.frame((kernwgta1/pscore1a1)/sum(kernwgta1/pscore1a1)) 
    relative_wgt_a1m0 <- as.data.frame(((kernwgta1*pscore2a0)/(pscore2a1*pscore1a0))/sum((kernwgta1*pscore2a0)/(pscore2a1*pscore1a0)))
    relative_wgt_a0m0 <- as.data.frame((kernwgta0/pscore1a0)/sum(kernwgta0/pscore1a0))
    relative_wgt_a0m1 <- as.data.frame(((kernwgta0*pscore2a1)/(pscore2a0*pscore1a1))/sum((kernwgta0*pscore2a1)/(pscore2a0*pscore1a1)))
  }
  
  if (ATET==TRUE){
    ya1m1=sum(y*kernwgta1)/sum(kernwgta1)
    ya1m0=sum(y*kernwgta1*pscore2a0*pscore1a1/(pscore2a1*pscore1a0))/sum(kernwgta1*pscore2a0*pscore1a1/(pscore2a1*pscore1a0))
    ya0m0=sum(y*kernwgta0*pscore1a1/pscore1a0)/sum(kernwgta0*pscore1a1/pscore1a0)
    ya0m1=sum(y*kernwgta0*pscore2a1/pscore2a0)/sum(kernwgta0*pscore2a1/pscore2a0)
    relative_wgt_a1m1 <- as.data.frame(kernwgta1/sum(kernwgta1))
    relative_wgt_a1m0 <- as.data.frame(((kernwgta1*pscore2a0*pscore1a1)/(pscore2a1*pscore1a0))/sum((kernwgta1*pscore2a0*pscore1a1)/(pscore2a1*pscore1a0)))
    relative_wgt_a0m1 <- as.data.frame(((kernwgta0*pscore1a1)/pscore1a0)/sum((kernwgta0*pscore1a1)/pscore1a0))
    relative_wgt_a0m0 <- as.data.frame(((kernwgta0*pscore2a1)/pscore2a0)/sum((kernwgta0*pscore2a1)/pscore2a0))
  }
  
  relative_wgt=data.frame(relative_wgt_a1m1,relative_wgt_a1m0,relative_wgt_a0m1,relative_wgt_a0m0)
  colnames(relative_wgt) <- c("relative_wgt_a1m1","relative_wgt_a1m0","relative_wgt_a0m1","relative_wgt_a0m0")
  
  results <- list(
    # return average total effect, natural direct effect under treatment and non-treatment, natural indirect effect under treatment and non-treatment, and the number of abnormal GPS values
    effect=c(ya1m1-ya0m0, ya1m1-ya0m1, ya1m0-ya0m0, ya1m1-ya1m0, ya0m1-ya0m0, sum(1-ind)),
    relative_wgt=relative_wgt
    ) 
  
  return(results)
  
}









## Function 2. Bootstrap samples generation

#' @param y Dependent variable, must not contain missings.
#' @param a Continuous treatment, must not contain missings.
#' @param m Continuous mediator, may be a scalar or a vector, must not contain missings.
#' @param x_ax Data matrix of pre-treatment covariates selected in the treatment GPS model given specified covariates.
#' @param x_amx Data matrix of pre-treatment covariates selected in the treatment GPS model given specified covariates and the mediator.
#' @param a0 Value of \code{a} under non-treatment. Effects are based on pairwise comparisons, i.e. differences in potential outcomes evaluated at \code{a1} and \code{a0}.
#' @param a1 Value of \code{a} under treatment. Effects are based on pairwise comparisons, i.e. differences in potential outcomes evaluated at \code{a1} and \code{a0}.
#' @param ATET If FALSE, the average treatment effect (ATE) and the corresponding direct and indirect effects are estimated. If TRUE, the average treatment effect on the treated (ATET)  and the corresponding direct and indirect effects are estimated. Default is FALSE.
#' @param trim Trimming rule for discarding observations with extreme generalized propensity scores. If \code{lognorm=FALSE}, observations with f(A=\code{a1}|M,X)<\code{trim} or f(A=\code{a0}|M,X)<\code{trim} are dropped, with f denoting the generalized propensity score (or conditional density of treatment). If \code{lognorm=TRUE}, then \code{trim} corresponds to the share of lowest f(A=\code{a1}|M,X) or f(A=\code{a0}|M,X), respectively, that are dropped.
#' @param lognorm If FALSE, a linear model with normally distributed errors is assumed for generalized propensity score estimation. If TRUE, a lognormal model is assumed. Default is FALSE.
#' @param bw Bandwith for the second order Epanechnikov kernel functions of the treatment. If set to NULL, the rule of thumb for Epanechnikov kernels is used for bandwidth computation. Default is NULL.
#' @param boot Number of bootstrap replications for estimating standard errors. Default is 1999.
#' @param cluster A cluster ID for block or cluster bootstrapping when units are clustered rather than iid. Must be numerical. Default is NULL (standard bootstrap without clustering).





bootstrap.mediation.cont<-function(y,a,m,x_ax,x_amx,a0,a1,ATET=FALSE,trim=0.05,lognorm=FALSE,bw,boot,cluster=NULL){ #boot=1999
  
  if (is.null(cluster)){
    obs<-length(y) 
    bsamples=matrix(NA,boot,6) 
    for(i in 1:boot){
      sboot<-sample(1:obs,obs,TRUE) 
      yb=y[sboot] 
      ab<-a[sboot] 
      if (is.null(ncol(m))) mb<-m[sboot]
      if (is.null(ncol(m))==0) mb<-m[sboot,]
      if (is.null(ncol(x_ax))) x_ax_b<-x_ax[sboot]
      if (is.null(ncol(x_ax))==0) x_ax_b<-x_ax[sboot,]
      if (is.null(ncol(x_amx))) x_amx_b<-x_amx[sboot]
      if (is.null(ncol(x_amx))==0) x_amx_b <- x_amx[sboot,]
      bsamples[i,]=c(mediation.cont(y=yb,a=ab,m=mb,x_ax=x_ax_b,x_amx=x_amx_b,a0=a0,a1=a1,ATET=ATET,trim=trim,lognorm=lognorm,bw=bw)$effect)
    }
  }
  
  if (is.null(cluster)==0){
    temp<-sort(cluster); clusters<-min(cluster)
    for (i in 1:length(temp)){
      if (temp[i]>max(clusters)) clusters=c(clusters,temp[i])
    }
    key=cluster; bsamples=c(); temp=c()
    obs<-length(clusters)
    while(length(temp)<boot){
      sboot<-sample(clusters,obs,TRUE)
      ab<-c(); yb<-c(); xb<-c() ; mb=c()
      for (k in 1:length(sboot)) {
        ab<-c(ab,d[key==sboot[k]]); yb<-c(yb,y[key==sboot[k]])
        if (is.null(ncol(m))) mb<-c(mb,m[key==sboot[k]])
        if (is.null(ncol(m))==0) mb=rbind(mb,m[key==sboot[k],])
        if (is.null(ncol(x_ax))) x_ax_b<-c(x_ax_b,x_ax[key==sboot[k]])
        if (is.null(ncol(x_ax))==0) x_ax_b=rbind(x_ax_b,x_ax[key==sboot[k],])
        if (is.null(ncol(x_amx))) x_amx_b<-c(x_amx_b,x_amx[key==sboot[k]])
        if (is.null(ncol(x_amx))==0) x_amx_b=rbind(x_amx_b,x_amx[key==sboot[k],])
      }
      est=c(mediation.cont(y=yb,a=ab,m=mb,x_ax=x_ax_b,x_amx=x_amx_b,a0=a0,a1=a1,ATET=ATET,trim=trim,lognorm=lognorm,bw=bw)$effect)
      bsamples<-rbind(bsamples, est)
      temp<-c(temp,1)
    }
  }
  
  bna=apply(bsamples, 1, sum) 
  bsamples=bsamples[is.na(bna)==0,] 
  if (sum(is.na(bna))>0) cat("Warning: ",sum(is.na(bna)>0)," bootstrap sample(s) dropped due to NA's")
  
  return(bsamples)
  
}











## Function 3. GPS-weighted mediation effect estimates 

#' @param y Dependent variable, must not contain missings.
#' @param a Continuous treatment, must not contain missings.
#' @param m Continuous mediator, may be a scalar or a vector, must not contain missings.
#' @param x_ax Data matrix of pre-treatment covariates selected in the treatment GPS model given specified covariates.
#' @param x_amx Data matrix of pre-treatment covariates selected in the treatment GPS model given specified covariates and the mediator.
#' @param a0 Value of \code{a} under non-treatment. Effects are based on pairwise comparisons, i.e. differences in potential outcomes evaluated at \code{a1} and \code{a0}.
#' @param a1 Value of \code{a} under treatment. Effects are based on pairwise comparisons, i.e. differences in potential outcomes evaluated at \code{a1} and \code{a0}.
#' @param ATET If FALSE, the average treatment effect (ATE) and the corresponding direct and indirect effects are estimated. If TRUE, the average treatment effect on the treated (ATET)  and the corresponding direct and indirect effects are estimated. Default is FALSE.
#' @param trim Trimming rule for discarding observations with extreme generalized propensity scores. If \code{lognorm=FALSE}, observations with f(A=\code{a1}|M,X)<\code{trim} or f(A=\code{a0}|M,X)<\code{trim} are dropped, with f denoting the generalized propensity score (or conditional density of treatment). If \code{lognorm=TRUE}, then \code{trim} corresponds to the share of lowest f(A=\code{a1}|M,X) or f(A=\code{a0}|M,X), respectively, that are dropped.
#' @param lognorm If FALSE, a linear model with normally distributed errors is assumed for generalized propensity score estimation. If TRUE, a lognormal model is assumed. Default is FALSE.
#' @param bw Bandwith for the second order Epanechnikov kernel functions of the treatment. If set to NULL, the rule of thumb for Epanechnikov kernels is used for bandwidth computation. Default is NULL.
#' @param boot Number of bootstrap replications for estimating standard errors. Default is 1999.
#' @param cluster A cluster ID for block or cluster bootstrapping when units are clustered rather than iid. Must be numerical. Default is NULL (standard bootstrap without clustering).
#' Standard errors are obtained by bootstrapping the effects.
#' @return A medweightcont object contains two components, \code{results} and \code{ntrimmed}:
#' @return \code{results}: a 3X5 matrix containing the effect estimates in the first row ("effects"), standard errors in the second row ("se"), and p-values in the third row ("p-value").
  #' The first column provides the total effect, namely the average treatment effect (ATE) if \code{ATET=FALSE} or the average treatment effect on the treated (ATET), i.e. those with D=\code{d1}, if \code{ATET=TRUE}.
  #' The second and third columns provide the direct effects under treatment and control, respectively ("dir.treat", "dir.control"). The fourth and fifth columns provide the indirect effects under treatment and control, respectively ("indir.treat", "indir.control").
#' @return \code{ntrimmed}: number of discarded (trimmed) observations due to extreme propensity score values.
#' @references Hsu, Y.-C., Huber, M., Lee, Y.-Y., Pipoz, L.: (2018): "Direct and indirect effects of continuous treatments based on generalized propensity score weighting",  SES working paper 495, University of Fribourg.
#' @examples # A little example with simulated data (10000 observations)

#' n=10000
#' x=runif(n=n,min=-1,max=1)
#' a=0.25*x+runif(n=n,min=-2,max=2)
#' a=a-min(a)
#' m=0.5*a+0.25*x+runif(n=n,min=-2,max=2)
#' y=0.5*a+m+0.25*x+runif(n=n,min=-2,max=2)
#' # The true direct and indirect effects are all equal to 0.5
#' output=medweightcont(y,a,m,x, a0=2, a1=3, ATET=FALSE, trim=0.05, lognorm=FALSE, bw=NULL, boot=999)
#' round(output$results,3)
#' output$ntrimmed
#' @importFrom stats binomial fitted.values glm lm pnorm sd rnorm dnorm quantile
#' @importFrom np npksum
#' @export



medweightcont<-function(y,a,m,x_ax,x_amx,a0,a1,ATET=FALSE,trim=0.05,lognorm=FALSE,bw=NULL,boot=1999,cluster=NULL){
  
  if(is.null(bw)) bw=sd(a)*2.34/(length(a)^0.2) 
  
  mediation_cont_res <- mediation.cont(y=y,a=a,m=m,x_ax=x_ax,x_amx=x_amx,a0=a0,a1=a1,ATET=ATET,trim=trim,lognorm=lognorm,bw=bw)
  
  temp=mediation_cont_res$effect
  ntrimmed=temp[length(temp)] 
  temp=temp[1:(length(temp)-1)] 
  temp2=bootstrap.mediation.cont(y=y,a=a,m=m,x_ax=x_ax,x_amx=x_amx,a0=a0,a1=a1,ATET=ATET,trim=trim,lognorm=lognorm,bw=bw,boot=boot,cluster=cluster)
  
  se=apply(temp2[,1:(ncol(temp2)-1)], 2, sd) 
  temp3=2*pnorm(-abs(temp/se)) 
  boot_results=rbind(temp, se, temp3)
  
  if (ATET==FALSE) colnames(boot_results)=c("ATE", "dir.treat", "dir.control", "indir.treat", "indir.control")
  if (ATET==TRUE) colnames(boot_results)=c("ATET", "dir.treat", "dir.control", "indir.treat", "indir.control")
  rownames(boot_results)=c("effect", "se", "p-value")
  
  res_list <- list(
    boot_results=boot_results,
    ntrimmed=ntrimmed,
    relative_wgt=as.data.frame(mediation_cont_res$relative_wgt)
    )
  
  return(res_list)
  
}








  



