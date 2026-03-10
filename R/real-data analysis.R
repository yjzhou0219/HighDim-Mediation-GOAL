source("GOAL_deviance.R")
source("mediation_effect_continuous.R")

library(data.table)
library(dplyr)
library(tidyverse)
library(psych)
library(openxlsx)
library(lubridate)
library(lattice)
library(nnet)
library(MASS)
library(VIM)
library(mediation)
library(foreach)
library(parallel)
library(doParallel)
library(devtools)
library(remotes)
library(sjstats)







##___________________________________________________________________________________________
# 1 exclusion ####


### 1.1 withdraw ####

withdraw_id <- fread("withdraw98410_214_20240627.txt") 
ukb_all$exc <- 0
ukb_all$exc[match(withdraw_id$V1,ukb_all$eid)] <- 1

ukb_data_allpar <- subset(ukb_all,ukb_all$exc==0) 



### 1.2 pregnancy ####

ukb_data_allpar$preg_exc <- 0
ukb_data_allpar$preg_exc[ukb_data_allpar$n.3140.0.0=="Yes"|ukb_data_allpar$n.3140.0.0=="Unsure"] <- 1
sum(ukb_data_allpar$preg_exc==1)
ukb_data_allpar <- subset(ukb_data_allpar,preg_exc==0) 



# 2 remove cases with missing data on FINDRISC (exposure) ####


### 1 age: n.21022.0.0
ukb_data_exp <- filter(ukb_data_allpar,n.21022.0.0>=40 & n.21022.0.0<70) 

### 2 bmi
ukb_data_exp <- ukb_data_allpar[complete.cases(ukb_data_allpar$n.21001.0.0),] 

### 3 waist circumference
ukb_data_exp <- ukb_data_exp[complete.cases(ukb_data_exp$n.48.0.0),] 

### 4 activity
ukb_data_exp <- ukb_data_exp[complete.cases(ukb_data_exp$n.22034.0.0),] 

### 5 vegetable and fruit intake, n.1289.0.0, n.1299.0.0 ,n.1309.0.0
recoding_rule_vage_fruit <- function(x){
  case_when(
    x<0&x>-4~NA,
    x==-10~0.5,
    TRUE~x
  )
}

ukb_data_exp <- ukb_data_exp %>% mutate(
  across(all_of(c(28:30)),recoding_rule_vage_fruit) 
)

ukb_data_exp <- ukb_data_exp[complete.cases(ukb_data_exp[,c(28:30)]),]

## 6 history of hypertension medication, n.6153.0.0, n.6177.0.2
recoding_rule_ht <- function(x){  
  case_when(
    x=="Cholesterol lowering medication" ~ 1,
    x=="Blood pressure medication" ~ 2,
    x=="Insulin" ~ 3,
    x=="Hormone replacement therapy" ~ 4,
    x=="Oral contraceptive pill or minipill" ~ 5,
    x=="None of the above" ~ -7,
    x=="Do not know" ~ -1,
    x=="Prefer not to answer" ~ -3
  )
}

ukb_data_exp <- ukb_data_exp %>% 
  mutate(
    across(all_of(c(361:367)), recoding_rule_ht)
  )

ukb_data_exp <- subset(ukb_data_exp,ukb_data_exp$n.31.0.0==1|(ukb_data_exp$n.31.0.0==0&!is.na(ukb_data_exp$n.6153.0.0)))
ukb_data_exp <- subset(ukb_data_exp,ukb_data_exp$n.31.0.0==0|(ukb_data_exp$n.31.0.0==1&!is.na(ukb_data_exp$n.6177.0.0)))


### 7 history of hyperglycemia: n.2443.0.0, n.120007.0.0

ukb_data_exp$n.2443.0.0[ukb_data_exp$n.2443.0.0<0] <- NA

recoding_rule_120007 <- function(x){  
  case_when(
    x=="Yes"~1,
    x=="No"~0,
    x=="Do not know" ~ -1,
    x=="Prefer not to answer" ~ -3,
    x=="" ~ NA
  )
}
which(names(ukb_data_exp)=="n.120007.0.0")
ukb_data_exp <- ukb_data_exp %>% 
  mutate(
    across(all_of(c(1035)), recoding_rule_120007)
  )
class(ukb_data_exp$n.120007.0.0)
table(ukb_data_exp$n.120007.0.0)

ukb_data_exp$n.120007.0.0[ukb_data_exp$n.120007.0.0<0] <- NA


## integrate
ukb_data_exp$hyperglycemia_history <- apply(ukb_data_exp[,c(64,1035)], 1, function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else if (any(x == 1, na.rm = TRUE)) {
    return(1)
  } else if (any(x == 0, na.rm = TRUE)) {
    return(0)
  } 
})

ukb_data_exp <- ukb_data_exp[complete.cases(ukb_data_exp$hyperglycemia_history),] 


### 8 family history of diabetes, n.20107.0.0, n.20110.0.0, n.20111.0.0

recoding_rule_family_diabetes <- function(x){ 
  case_when(
    x=="Hip fracture" ~ 14,
    x=="Prostate cancer" ~ 13,
    x=="Severe depression" ~ 12,
    x=="Parkinsons disease" ~ 11,
    x=="Alzheimers disease/dementia" ~ 10,
    x=="Diabetes" ~ 9,
    x=="High blood pressure" ~ 8,
    x=="Chronic bronchitis/emphysema" ~ 6,
    x=="Breast cancer" ~ 5,
    x=="Bowel cancer" ~ 4,
    x=="Lung cancer" ~ 3,
    x=="Stroke" ~ 2,
    x=="Heart disease" ~ 1,
    x=="Do not know (group 1)" ~ -11,
    x=="Prefer not to answer (group 1)" ~ -13,
    x=="None of the above (group 1)" ~-17,
    x=="Do not know (group 2)" ~ -21,
    x=="Prefer not to answer (group 2)" ~ -23,
    x=="None of the above (group 2)" ~-27
  )
}

ukb_data_exp <- ukb_data_exp %>% 
  mutate(
    across(all_of(c(380,390,401)), recoding_rule_family_diabetes)
  )

ukb_data_exp <- ukb_data_exp[complete.cases(ukb_data_exp[,c(380,390,401)]),] 

ukb_data_exp <- subset(ukb_data_exp,n.20107.0.0!=-11&n.20107.0.0!=-13&n.20110.0.0!=-11&n.20110.0.0!=-13&n.20111.0.0!=-11&n.20111.0.0!=-13)

ukb_data_exp <-  ukb_data_exp %>% 
  mutate(family_diabetes=ifelse(n.20107.0.0==9|n.20110.0.0==9|n.20111.0.0==9,1,0))








##___________________________________________________________________________________________
# 3 remove cases with missing data on the outcome (cancers) and cases with baseline diseases ####


ukb_data_exp$recruit_date <- ymd(ukb_data_exp$s.53.0.0, format = "%d-%b-%y") 

data_cancer <- ukb_data_exp[,c(1,65,134:139,368:379,413:1032,226,1043)]





### diag_icd10, s.41270.0.0 : s.41270.0.258 ####

## indicator

cancer_code_icd10 <- c(as.character(paste("C0",0:9,sep = "")),as.character(paste("C",c(10:43,45:97),sep = "")))

cancer_icd10 <- apply(data_cancer[,29:287], 1, function(x) {
  any(sapply(x, function(val) {
    ifelse (!is.na(val) && nchar(val) >= 3L, substr(val, 1, 3) %in% cancer_code_icd10, FALSE)
  }))
})

cancer_icd10[apply(data_cancer[,29:287]=="", 1, all)] <- NA  
data_cancer$cancer_icd10 <- as.integer(cancer_icd10)


## date, s.41280.0.0 : s.41280.0.258

data_cancer[,c(335:593)] <- lapply(data_cancer[,c(335:593)], function(x) {
  dmy(x)
})

pattern_icd10 <- paste0("^(", paste(cancer_code_icd10, collapse = "|"), ")") 
indices_icd10 <- apply(data_cancer[,29:287], 1, function(x) {
  which(grepl(pattern_icd10, x)) 
})

date_min_41280 <- rep(NA, nrow(data_cancer))
date_min_41280 <- as.Date(date_min_41280)

for(i in 1:length(indices_icd10)){
  if(length(indices_icd10[[i]])!=0){
    pos_date <- indices_icd10[[i]]+334 
    date_number_df <- unlist(lapply(data_cancer[i,pos_date],as.Date))
    col_index <- which(date_number_df==min(date_number_df,na.rm = T)) 
    col_index <- as.numeric(col_index[1]) 
    if(length(indices_icd10[[i]])>1){
      date_min_41280[i] <- data_cancer[i,pos_date][,col_index]
    }else{
      date_min_41280[i] <- data_cancer[i,pos_date]
    }
  }else{
    date_min_41280[i] <- NA
  }
}

data_cancer$date_cancer_i10_first <- date_min_41280





### diag_icd9, s.41271.0.0 : s.41271.0.46 ####


## indicator
pattern_icd9 <- paste0("^(", paste(c(140:172,174:208), collapse = "|"), ")")  

cancer_icd9 <- apply(data_cancer[,c(288:325,329:334)], 1, function(x) {
  non_empty_values <- x[!is.na(x) & x != ""]
  if (length(non_empty_values) == 0) {
    return(NA)
  }
  matches <- grepl(pattern_icd9, non_empty_values, ignore.case = FALSE)
  return(any(matches))
})

data_cancer$cancer_icd9_indicate <- as.integer(cancer_icd9)


## date_icd9

data_cancer[,c(594:640)] <- lapply(data_cancer[,c(594:640)], function(x) {
  dmy(x)
})

cancer_icd9 <- data_cancer[!is.na(data_cancer$cancer_icd9_indicate),]

indices_icd9 <- apply(cancer_icd9[,c(1,594:640)], 1, function(x) {
  which(grepl(pattern_icd9, x)) 
})

date_min_41281 <- rep(NA, nrow(cancer_icd9))
date_min_41281 <- as.Date(date_min_41281)

for(i in 1:length(indices_icd9)){
  if(length(indices_icd9[[i]])!=0){
    pos_date <- indices_icd9[[i]]+593 
    date_number_df <- unlist(lapply(cancer_icd9[i,pos_date],as.Date))
    col_index <- which(date_number_df==min(date_number_df,na.rm = T)) 
    col_index <- as.numeric(col_index[1]) 
    if(length(indices_icd9[[i]])>1){
      date_min_41281[i] <- cancer_icd9[i,pos_date][,col_index]
    }else{
      date_min_41281[i] <- cancer_icd9[i,pos_date]
    }
  }else{
    date_min_41281[i] <- NA
  }
}

cancer_icd9$date_cancer_i9_first <- date_min_41281
data_cancer <- left_join(data_cancer,cancer_icd9[,c(1,646)],by="eid")




### self_reported ####


#### n.20001.0.0, n.20007.0.0

## indicate
cancer_code_self <- c(1001:1059,1063:1072,1074:1088)

cancer_self <- apply(data_cancer[,3:8], 1, function(x) {
  any(sapply(x, function(val) {
    ifelse (!is.na(val) && nchar(val) >= 4L, substr(val, 1, 4) %in% cancer_code_self, FALSE)
  }))
})

cancer_self[apply(is.na(data_cancer[,3:8]), 1, all)] <- NA  
data_cancer$cancer_self <- as.integer(cancer_self)


## age

data_cancer[,c(15:20)] <- lapply(data_cancer[,c(15:20)], function(x) {
  as.numeric(x)
})

cancer_self <- data_cancer[!is.na(data_cancer$cancer_self),]

pattern_self <- paste0("^(", paste0(cancer_code_self,collapse = "|"),")")

indices_self <- apply(cancer_self[,c(1,3:8)], 1, function(x) {
  which(grepl(pattern_self, x)) 
})

age_min_self <- rep(NA, nrow(cancer_self))

for(i in 1:length(indices_self)){
  if(length(indices_self[[i]])!=0){
    pos_age <- indices_self[[i]]+14 
    age_number_df <- cancer_self[i,pos_age]
    col_index <- which(age_number_df==min(age_number_df,na.rm = T)) 
    col_index <- as.numeric(col_index[1]) 
    if(!is.na(col_index)){
      if(length(indices_self[[i]])>1){
        age_min_self[i] <- cancer_self[i,pos_age][,col_index]
      }else{
        age_min_self[i] <- cancer_self[i,pos_age]
      }
    }else{
      age_min_self[i] <- NA
    }
  }
}

cancer_self$age_cancer_self_first <- age_min_self
data_cancer <- left_join(data_cancer,cancer_self[,c(1,648)],by="eid")





#### 2453 Cancer diagnosed by doctor

data_cancer$n.2453.0.0[!is.na(data_cancer$n.2453.0.0)&data_cancer$n.2453.0.0<0] <- NA



### other: Date of cancer diagnosis, s.40005.0.0 ####

data_cancer$s.40005.0.0 <- dmy(data_cancer$s.40005.0.0)

data_cancer$diagnosis_other <- ifelse(is.na(data_cancer$s.40005.0.0),0,1)





### integration of cancer data ####


which(names(data_cancer)=="cancer_icd10") # 643
which(names(data_cancer)=="cancer_icd9_indicate") # 645
which(names(data_cancer)=="cancer_self") # 647
which(names(data_cancer)=="diagnosis_other") # 649

data_cancer$overall_cancer <- apply(data_cancer[,c(643,645,647,649)], 1, function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else if (any(x == 1, na.rm = TRUE)) {
    return(1)
  } else if (any(x == 0, na.rm = TRUE)) {
    return(0)
  } 
})


## date

data_cancer <- data_cancer %>%
  mutate(date_first_cancer = case_when(
    is.na(date_cancer_i10_first) & is.na(date_cancer_i9_first) & is.na(s.40005.0.0) ~ NA,        
    TRUE ~ pmin(date_cancer_i10_first, date_cancer_i9_first, s.40005.0.0, na.rm = TRUE)  
  ))

data_cancer$overall_cancer[!is.na(data_cancer$date_first_cancer)] <- 1


## age
data_cancer <- data_cancer %>%
  mutate(diff = case_when(
    is.na(date_first_cancer) & is.na(age_cancer_self_first) ~ NA,   
    !is.na(date_first_cancer) & is.na(age_cancer_self_first) ~ as.numeric(difftime(date_first_cancer, recruit_date, units = "days"))/365,  
    is.na(date_first_cancer) & !is.na(age_cancer_self_first) ~ age_cancer_self_first - n.21022.0.0,  
    !is.na(date_first_cancer) & !is.na(age_cancer_self_first) ~ pmin(as.numeric(difftime(date_first_cancer, recruit_date, units = "days"))/365, age_cancer_self_first - n.21022.0.0)  
  ))



data_cancer$overall_cancer[!is.na(data_cancer$diff)] <- 1



###  final integration ####


ukb_data_out <- ukb_data_exp[,-c(65,134:139,368:379,413:1032)]
ukb_data_out <- left_join(ukb_data_out,data_cancer[,c(1,650,652)],by="eid")

ukb_data_out <- subset(ukb_data_out,is.na(diff)|(!is.na(diff)&diff>=(-10)))



## 
ukb_data_simp <- ukb_data_out[,-c(3:11,13,16:20,47,64:70,72:73,75:76,78:87,89,91,94:95,99:103,105:114,
                                  116:123,125:127,128:180,183,188:199,201:205,218,220:221,223,236,250,259:269,285:286,289:291,293:295,
                                  338:344,346,351:353,355:357,359:398,405)]









##___________________________________________________________________________________________
# 4 complete cases for covariates ####


# race, n.21000.0.0
ukb_data_cov <- ukb_data_simp[complete.cases(ukb_data_simp$n.21000.0.0),]
ukb_data_cov <- subset(ukb_data_simp,n.21000.0.0>=0) 


# income, n.738.0.0
ukb_data_cov <- ukb_data_cov[complete.cases(ukb_data_cov$n.738.0.0),] 
ukb_data_cov <- subset(ukb_data_cov,n.738.0.0>=0) 


# education, n.6138.0.0
ukb_data_cov <- subset(ukb_data_cov,n.6138.0.0!=-3) 


# TDI, n.22189.0.0
ukb_data_cov <- ukb_data_cov[complete.cases(ukb_data_cov$n.22189.0.0),] 


# drinking, n.1558.0.0
ukb_data_cov <- subset(ukb_data_cov,n.1558.0.0>=0) 


# occupation
ukb_data_cov$occup <- ukb_data_cov$n.6142.0.0
ukb_data_cov$occup[!is.na(ukb_data_cov$n.20119.0.0)] <- ukb_data_cov$n.20119.0.0[!is.na(ukb_data_cov$n.20119.0.0)]
ukb_data_cov <- subset(ukb_data_cov,occup!=-3) 



ukb_data_cov_simp <- ukb_data_cov[,-c(12,53:56,58,63:64,138:139,163,165)]




### remove missingness of covariates ####

cols_for_rule1 <- c(4,6,10:11,16:23,29:44,46,54,56,57:71,107:126,172:174)
cols_for_rule2 <- c(7:9,12:15,24:28)
cols_for_rule3 <- 45
cols_for_rule4 <- c(52:53,168:169,176)

recoding_rule_1 <- function(x){  
  case_when(
    x<0 ~ NA,
    TRUE ~ x
  )
}

recoding_rule_2 <- function(x){  
  case_when(
    x<0&x>-4 ~ NA,
    x==-10 ~ 0.5,
    TRUE ~ x
  )
}

recoding_rule_3 <- function(x){  
  case_when(
    x<0 ~ NA,
    x==99 ~ 2,
    TRUE ~ x
  )
}

recoding_rule_4 <- function(x){  
  case_when(
    x<0&x>-4 ~ NA,
    x==-7 ~ 0,
    TRUE ~ x
  )
}

ukb_data_cov_recode <- ukb_data_cov_simp %>% 
  mutate(
    across(all_of(cols_for_rule1), recoding_rule_1),
    across(all_of(cols_for_rule2), recoding_rule_2),
    across(all_of(cols_for_rule3), recoding_rule_3),
    across(all_of(cols_for_rule4), recoding_rule_4)
  )


## record missingness

missingness <- rep(NA,ncol(ukb_data_cov_recode)-1)
for(i in 2:ncol(ukb_data_cov_recode)){
  missingness[i-1] <- sum(is.na(ukb_data_cov_recode[,i])|ukb_data_cov_recode[,i]=="")
}
missingness_df <- as.data.frame(missingness)

## manually examine missingness
missing_cols_cov <- c(3:4,6:11,15:28,30:32,35:38,41:51,53,55:70,74:146,148:155,157:160,162,164:167,177:178) 

ukb_data_comp <- subset(ukb_data_comp,ukb_data_comp$n.31.0.0==1|(ukb_data_comp$n.31.0.0==0&!is.na(ukb_data_comp$n.6153.0.0)))
ukb_data_comp <- subset(ukb_data_comp,ukb_data_comp$n.31.0.0==0|(ukb_data_comp$n.31.0.0==1&!is.na(ukb_data_comp$n.6177.0.0)))

ukb_data_comp <- ukb_data_cov_recode[complete.cases(ukb_data_cov_recode[,missing_cols_cov]),] 




### recode covariates and construct composite variables ####

### NLR
ukb_data_comp$NLR <- ukb_data_comp$n.30140.0.0/ukb_data_comp$n.30120.0.0

### physical activity
ukb_data_comp$physical_activity <- apply(ukb_data_comp[,c(172:174)], 1, function(x) {
  if (all(is.na(x))) {
    NA
  } else {
    sum(x, na.rm = TRUE)
  }
})

ukb_data_comp <- subset(ukb_data_comp,!is.na(ukb_data_comp$physical_activity))

ukb_data_comp <- ukb_data_comp[,-c(57,147,156,167,172:174,178)]



### age
ukb_data_comp <-  ukb_data_comp %>% 
  mutate(age_score=ifelse(n.21022.0.0<45,0,ifelse(n.21022.0.0<=54,2,ifelse(n.21022.0.0<=64,3,4))))


### family_hist_diabetes
ukb_data_comp$famliy_diabetes_hist_score <- ifelse(ukb_data_comp$family_diabetes==1,5,0)


### fruit/vegetable
ukb_data_comp$fru_veg_intake <- apply(ukb_data_comp[,c(12:14)],1,function(x){
  ifelse(any(x>0),1,0)
})
ukb_data_comp$fru_veg_intake_score <- ifelse(ukb_data_comp$fru_veg_intake==1,0,1)


## physical_act
ukb_data_comp <- ukb_data_comp %>%
  mutate(PA_indicator=ifelse((physical_activity*7)/60<4,0,1))
ukb_data_comp$PA_score <- ifelse(ukb_data_comp$physical_activity==1,0,2)


## blood pressure medication
ukb_data_comp$bp_medication <- ifelse((ukb_data_comp$n.31.0.0==0&ukb_data_comp$n.6153.0.0==2)|(ukb_data_comp$n.31.0.0==1&ukb_data_comp$n.6177.0.0==2),1,0)
ukb_data_comp$bp_medication_score <- ifelse(ukb_data_comp$bp_medication==1,2,0)


## history of high blood glucose 
ukb_data_comp$hyperglycemia_history_score <- ifelse(ukb_data_comp$hyperglycemia_history==1,5,0)


## BMI
ukb_data_comp <- ukb_data_comp%>%
  mutate(bmi_score=ifelse(n.21001.0.0<=25,0,ifelse(n.21001.0.0<=30,1,3)))


## waist
ukb_data_comp <- ukb_data_comp %>%
  mutate(waist_score = case_when(
    n.31.0.0 == 0 & n.48.0.0 < 80 ~ 0,
    n.31.0.0 == 0 & n.48.0.0 >= 80 & n.48.0.0 < 88 ~ 3,
    n.31.0.0 == 0 & n.48.0.0 >= 88 ~ 4,
    n.31.0.0 == 1 & n.48.0.0 < 94 ~ 0,
    n.31.0.0 == 1 & n.48.0.0 >= 94 & n.48.0.0 < 102 ~ 3,
    n.31.0.0 == 1 & n.48.0.0 >= 102 ~ 4
  ))


## FINDRISC 
ukb_data_comp$FINDRISC_score <- apply(ukb_data_comp[,c(173:174,176,178,180:183)],1,sum)
ukb_data_comp_exp <- ukb_data_comp


## education
ukb_data_comp_exp$edu <- ifelse(ukb_data_comp_exp$n.6138.0.0==1,3,ifelse(ukb_data_comp_exp$n.6138.0.0==0,1,2))


## income
ukb_data_comp_exp$income <- ifelse(ukb_data_comp_exp$n.738.0.0<=2,1,ifelse(ukb_data_comp_exp$n.738.0.0==3,2,3))


ukb_data_final <- ukb_data_comp_exp[,-c(4:6,12:25,33:45,47:49,53:54,57:71,78,119:123,132,151:153,156:158,160,162:167,169:170,172:183)]






##___________________________________________________________________________________________
# 5 final analysis ####


## 0 data processing ####


## exp: FINDRISC_score
## med: Apolipoprotein B
## out: cancer


Data_realexamp <- ukb_data_final[,c(2:78,80:97,99,101:103,98,100,79)] 

n=as.numeric(nrow(Data_realexamp)) 
p=as.numeric(ncol(Data_realexamp)-3) 

colnames(Data_realexamp) <- c(as.character(paste0("X",1:p)),"Y","Trt","Med")

var.list <- names(Data_realexamp)[-c((ncol(Data_realexamp)-2):ncol(Data_realexamp))] 



## 1 exploratory analysis ####

## interaction term

formula1 <- formula(paste("Y~Trt+Med+Trt*Med+",paste(var.list,collapse = "+")))
model1 <- lm(formula1, data = Data_realexamp)
summary(model1)

formula.out <- formula(paste("Y~Med+Trt+Trt*Med+",paste(var.list,collapse = "+")))
out.fit <- glm(formula.out,data=Data_realexamp,family = binomial("probit"))
formula.med <- formula(paste("Med~Trt+",paste(var.list,collapse = "+")))
med.fit <- lm(formula.med,data=Data_realexamp)
med.out <- mediate(med.fit,out.fit,treat = "Trt",mediator="Med",boot=T, treat.value = 6,control.value = 2) 
test.TMint(med.out)


## Nonlinear polynomial term
formula2 <- formula(paste("Y~Trt+Med+Trt*Med+Trt2+Med2+",paste(var.list,collapse = "+")))
model2 <- glm(formula2, data = Data_realexamp_temp)
summary(model2)

formula.out2 <- formula(paste("Y~Med+Trt+Trt*Med+Trt2+Med2+",paste(var.list,collapse = "+")))
out.fit2 <- glm(formula.out2,data=Data_realexamp,family = binomial("probit"))
formula.med2 <- formula(paste("Med~Trt+",paste(var.list,collapse = "+")))
med.fit2 <- lm(formula.med2,data=Data_realexamp)
med.out2 <- mediate(med.fit2,out.fit2,treat = "Trt",mediator="Med",boot=T, treat.value = 6,control.value = 2) 
summary(med.out2)




## 2 variable selection ####


n <- nrow(Data_realexamp)
p <- ncol(Data_realexamp)-3
var.list <- names(Data_realexamp)[1:p]


## GOAL_deviance

cat(sprintf("########################## 开始循环: %s ## \n", Sys.time()))
GOAL_deviance_res <- GOAL_deviance(n=n,Data=Data_realexamp,var.list=var.list,lambda_vec=c(-10,-5,-2,-1,-0.75,-0.5,-0.25,0.25,0.49),covar=NULL,gamma_convergence=2) 
cat(sprintf("## 结束循环: %s ##\n", Sys.time()))
Svar_ax_final<- GOAL_deviance_res$GOAL_AX$Svar_ax
Svar_amx_final<- GOAL_deviance_res$GOAL_AMX$Svar_amx

save(GOAL_deviance_res,file = "GOAL_deviance_res_RE.Rdata")




## 3 estimation ####

y=as.numeric(Data_realexamp[,"Y"])
a=as.numeric(Data_realexamp[,"Trt"])
m=as.numeric(Data_realexamp[,"Med"])


X_ax_var=GOAL_deviance_res$GOAL_AX$Svar_ax
X_amx_var=GOAL_deviance_res$GOAL_AMX$Svar_amx
X_ax <- as.matrix(Data_realexamp[,X_ax_var])
X_amx <- as.matrix(Data_realexamp[,X_amx_var])

exp_score <- seq(3,13,1) 

cat(sprintf("########################## 开始循环: %s ## \n", Sys.time()))
cl <- makeCluster(5)
registerDoParallel(cl)

boot_results_a1 <- foreach(j = 1:length(exp_score), .combine = rbind, .packages = c('np')) %dopar% {
  res_list <- medweightcont(y = y, a = a, m = m, x_ax = X_ax, x_amx = X_amx, a0 = 2, a1 = exp_score[j], ATET = FALSE, 
                            trim = 0.05, lognorm = F, bw = 2.34*(n^(-0.25)), boot = 500, cluster = NULL) 
  boot_res_ai <- as.data.frame(res_list$boot_results)
  rownames(boot_res_ai) <- paste(rownames(boot_res_ai), exp_score[j], sep = "_")
  ntrim <- res_list$ntrimmed
  relative_wgt <- as.data.frame(res_list$relative_wgt)
  list(
    ntrim=ntrim,
    boot_res_ai=boot_res_ai,
    relative_wgt_=relative_wgt
  )
}

stopCluster(cl)
cat(sprintf("## 结束循环: %s ##\n", Sys.time()))

GOAL_deviance_est_wp <- do.call(rbind,boot_results_a1[,2])
relative_wgt_list <- boot_results_a1[,3]
ntrim <- as.data.frame(do.call(rbind,boot_results_a1[,1]))
rownames(ntrim) <- paste("ntrim",exp_score,sep = "_")

save(GOAL_deviance_est_wp,file = "GOAL_deviance_est_wp.Rdata")
save(relative_wgt_list,file = "relative_wgt_list.Rdata")
write.csv(ntrim,file = "ntrim_GOAL_deviance_wp.csv",quote = F)



final_res_est <- GOAL_deviance_est_wp[grepl("effect", row.names(GOAL_deviance_est_wp)),] 
final_res_se <- GOAL_deviance_est_wp[grepl("se", row.names(GOAL_deviance_est_wp)),]
final_res_p <- GOAL_deviance_est_wp[grepl("p-value", row.names(GOAL_deviance_est_wp)),]
final_res_p

final_res_se_mean <- rbind(final_res_se,colMeans(final_res_se,na.rm = T)) #忽略缺失值
rownames(final_res_se_mean) <- c(rownames(final_res_se_mean)[1:11],"mean")

write.csv(final_res_est,file = "final_res_est.csv",quote = F)
write.csv(final_res_se_mean,file = "final_res_se_mean.csv",quote = F)
write.csv(final_res_p,file = "final_res_p.csv",quote = F)



