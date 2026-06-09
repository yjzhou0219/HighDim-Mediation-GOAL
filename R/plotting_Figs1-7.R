## Figures 1-7

#install.packages(ggplot2)
#install.packages(plyr)
#install.packages(data.table)
#install.packages(psych)
#install.packages(tidyr)
#install.packages(ggh4x)



library(ggplot2)
library(plyr)
library(data.table)
library(psych)
library(tidyr)
library(ggh4x)




## 1 simulation analysis: variable selection -- Figs 1-2 ####


###  result_organization  ####

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





###  plot  ####

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







## 2 simulation analysis: bias estimation comparison -- Figs 3-4 ####


### Successively read in the previously calculated abias data, which consists of 20 rows (corresponding to 20 exposure assessment points) and four columns (corresponding to different causal mediation effects, namely NDE1, NDE0, NIE1, and NIE0) for each specific scenario and method.

## Repeat the import process for all scenarios and methods. 

scenario <- c("S1","S2","S3")
sample <- c("n2000_rho0","n5000_rho0")
methods <- c("GOAL","DoubleML","AdaLASSO","LASSO")


res_abias_comb <- data.frame()
for (sf_idx in scenario) {
  for (s in sample) {
    for (met in methods) {
      ## Successively read in the previously calculated abias data, which consists of 20 rows (corresponding to 20 exposure assessment points) and four columns (corresponding to different causal mediation effects, namely NDE1, NDE0, NIE1, and NIE0) for each specific scenario and method.
      res_abias <- read.csv(paste0(sf_idx,"_sim_res_abias_",met,"_",s,".csv"),row.names = F)
      
      estimand_name <- c("NDE1","NDE0","NIE1","NIE0")
      colnames(res_abias) <- estimand_name
      res_abias$Exposure <- a1_all_20
      
      res_abias_long <- res_abias%>%
        pivot_longer(
          cols = contains(estimand_name), # 选中所有指标列
          names_to = "Estimand", 
          values_to = "ABias"
        ) %>%
        select(Estimand, Exposure, ABias)
      
      res_abias_long <- data.frame(Scenario=sf_idx,Sample=s,Method=met,res_abias_long)
      res_abias_comb <- rbind(res_abias_comb,res_abias_long)
    }
  }
}



###  plot  #### 


my_colors <- c("GOAL" = "#CB250F", "DoubleML" = "#3E6EC7", "AdaLASSO"="#E7A255", "LASSO" = "#238D3E") #, "Benchmark (Full)"="#757575"


## Using the setting of S1 and n=2000 as an example
df_subset <- res_abias_comb %>% filter(Scenario == "S1" & Sample == "n2000_rho0")
head(df_subset)
nrow(df_subset)
df_subset <- df_subset[,-c(1,2)]

method_levels <- c("GOAL", "DoubleML", "AdaLASSO", "LASSO")
estimand_levels <- c("NDE1", "NDE0", "NIE1", "NIE0")

df_subset$Method   <- factor(df_subset$Method, levels = method_levels)
df_subset$Estimand <- factor(df_subset$Estimand, levels = estimand_levels)



## plot

pdf("plot_S1_n2000_bias.pdf",width = 6,height = 2)

ggplot(df_subset, aes(x = Method, y = abs(ABias), fill = Method)) +
  
  # Add a horizontal reference line at 0 (zero line)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey10", linewidth = 0.5) +
  
  geom_boxplot(width = 0.75, alpha = 0.7, linewidth=0.3,
               aes(fill = Method),color = "grey20",
               outlier.shape = NA) +
  
  
  # Layout: 1 row and 4 columns
  facet_wrap(~ Estimand, nrow = 1,scales = "free_y") +
  
  # Set the y-axis range for the NDE and NIE panels separately
  facetted_pos_scales(
    y = list(
      Estimand %in% c("NDE1","NDE0") ~
        scale_y_continuous(
          limits = c(-0.01,0.1), 
          breaks = seq(0,0.1,0.05)
        ),
      
      Estimand %in% c("NIE1","NIE0") ~
        scale_y_continuous(
          limits = c(-0.002,0.015), 
          breaks = seq(0,0.015,0.005)
        )
    )
  )+
  
  scale_fill_manual(values = my_colors) +
  
  labs(y = NULL, x = NULL) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none", 
    axis.text.x = element_blank(), 
    axis.text.y = element_text(colour = "gray10",size = 5), 
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "white", color = "gray30"),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80",size=0.3),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.25,"lines")
  )


dev.off()









## 3 real-data analysis: effect distribution plot -- Figs 5-6 ####


###  result_organization  ####

setwd("./1 GOAL_deviance/")

est <- fread("final_res_est.csv")
est <- est[,-1]

se <- fread("final_res_se_mean.csv")
se <- se[-12,-1]

colnames(est)=colnames(se)=c("NDE1","NDE0","NIE1","NDE1")





###  plot  ####

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









## Density histogram: Fig 7 ####


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



