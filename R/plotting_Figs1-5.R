## Figures 1-5 ####
library(ggplot2)
library(plyr)
library(data.table)
library(psych)




### 1 result_organization ####

## Using GOAL as an example here. The same applies to AdaLASSO and LASSO results.

load("S1_sost_GOAL_deviance_res.Rdata")

Svar_ax_final <- GOAL_S1_res[,6]
Svar_amx_final <- GOAL_S1_res[,7]

variable_count <- as.character(1:30)

count_result_ax <- sapply(variable_count, function(x) {
  sum(sapply(Svar_ax_final, function(y) sum(y == x)))
})

count_result_amx <- sapply(variable_count, function(x) {
  sum(sapply(Svar_amx_final, function(y) sum(y == x)))
})

count_result_ax <- data.frame(count_result_ax)
count_result_ax$proportion <- count_result_ax$count_result_ax/100

count_result_amx <- data.frame(count_result_amx)
count_result_amx$proportion <- count_result_amx$count_result_amx/100

write.csv(count_result_ax,file = "S1_sost_GOAL_deviance_ax.csv")
write.csv(count_result_amx,file = "S1_sost_GOAL_deviance_amx.csv")



### 2 plot ####

colors <- c("GOAL"="#DE6E53","AdaLASSO"="#10C7FC","LASSO"="#229019","Reference"="#0032A0")
group <- c("GOAL", "AdaLASSO", "LASSO")

S1_amx_list <- list.files(pattern = "S1.*_amx\\.csv$")
S1_amx_list <- S1_amx_list[c(2,1,3)] # Adjust the order to match the 'group' vector
dat_df <- data.frame()

for (i in 1:3){
  dat <- fread(S1_amx_list[i])
  dat$Varname <- factor(dat$V1,levels = paste0("X",c(1,6,11,16,21,26,2:5,7:10,12:15,17:20,22:25,27:30)))
  dat$group <- group[i]
  dat_df <- rbind(dat_df,dat)
}

ref <- as.data.frame(matrix(c(paste0("X",1:30),c(rep(c(1,rep(0,4)),4),rep(0,10))),ncol = 2,byrow = F))
names(ref)[names(ref)=="V2"] <- "proportion"
ref$Varname <- factor(ref$V1,levels = paste0("X",c(1,6,11,16,21,26,2:5,7:10,12:15,17:20,22:25,27:30)))
ref$group <- "Reference"


dat_df_full <- rbind(dat_df[,-2],ref)
dim(dat_df_full)
dat_df_full$proportion <- as.numeric(dat_df_full$proportion)
head(dat_df_full)


## plot

pdf("selection_plot_S1_sost_amx_n2000rho0.pdf",width = 6,height = 4)
 
ggplot(dat_df_full,aes(x=Varname,y=proportion,group=group,color=as.factor(group)))+
  geom_line(linewidth=0.7)+
  scale_x_discrete(name="Covariate")+
  scale_y_continuous(name = "Proportion of covariates selected",breaks = seq(0,1,by=0.1))+
  scale_color_manual(values = colors)+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA,color = "grey30",linewidth = 0.6),
    axis.text.x = element_text(colour = "black",size = 7.2),
    axis.text.y = element_text(colour = "black",size = 8),
    axis.title.x = element_text(size=10,face = "bold"),
    axis.title.y = element_text(size=8,face = "bold"),
    legend.position = "inside",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.background = element_blank(),
    plot.margin = unit(c(0.3,0.1,0.2,0.2),"cm") ##top, right, bottom, left
  )

dev.off()







## Density histogram: Fig 3 ####


## exposure
describe(Data_realexamp$Trt)

breaks_exp=seq(0, 26, by = 2)

hist(
  Data_realexamp$Trt,
  breaks=breaks_exp, 
  xlim = c(0,25),
  xlab = "FINDRISC",
  ylab = "Frequency",
  right = T,
  col = "grey75",
  border = "white",
  main = title(mian=NULL)
)









## real-data analysis: effect distribution plot: Figs 4-5 ####


### 1 result_organization ####

setwd("./1 GOAL_deviance/")

est <- fread("final_res_est.csv")
est <- est[,-1]

se <- fread("final_res_se_mean.csv")
se <- se[-12,-1]

colnames(est)=colnames(se)=c("NDE1","NDE0","NIE1","NDE1")


### 2 plot ####

## Using NDE1 as an example here. The same applies to NDE0, NIE1, and NIE0 results.

NDE1_GOAL_deviance_res <- fread("NDE1_GOAL_deviance_res.csv")

pdf("RE_plot_NDE1_res_GOAL_deviance.pdf",width = 3,height = 4)

ggplot(NDE1_GOAL_deviance_res,aes(x=exp,y=est))+
  geom_hline(yintercept = 0,color="red",linetype = "dashed",linewidth=0.4)+
  geom_line(linewidth=0.4)+ 
  geom_point(size=1,colour="black",shape=21,fill="white")+
  geom_line(aes(y = est-1.96*se), colour = 'grey30', linetype = "dotted",linewidth=0.6) +  
  geom_line(aes(y = est+1.96*se), colour = 'grey30', linetype = "dotted",linewidth=0.6) + 
  scale_x_continuous(name="FINDRISC",limits = c(3,13),breaks = seq(4,12,2))+ 
  
  ## NDE
  scale_y_continuous(limits=c(-0.13,0.22), breaks = round(seq(-0.1,0.2,by=0.1),2))+  
  
  ## NIE
  #scale_y_continuous(limits=c(-0.022,0.048),breaks = seq(-0.02,0.04,by=0.02))+ 
  
  theme(
    panel.border=element_rect(color = "black", fill = NA),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = "black",size = 10), 
    axis.text.y = element_text(colour = "black",size = 10), 
    axis.title.x = element_text(size=10, face = "bold"),
    axis.title.y = element_blank(),
    plot.margin = unit(c(0.2,0.2,0.1,0.1),"cm") ##top, right, bottom, left
  )

dev.off()


