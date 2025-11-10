library(Matrix)
library(lme4)
library(performance) 
library(broom.mixed) 
library(dplyr)

###Table 1 FAS
library(tableone)

#data1<-read_csv("/data/20250511_RCT_data.csv")
data1<-X20250511_RCT_data
data1$group<-ifelse(data1$group==2,0,data1$group)

data1$group<-as.factor(data1$group)
data1$inferity_type<-as.factor(data1$inferity_type)
data1$fre_trans<-as.factor(data1$fre_trans)

data1 <- transform(data1, E2 = E2 / 3)

data1$inferity_type <- factor(data1$inferity_type,
                              levels = c(1, 2),
                              labels = c("old", "DOR"))
table(data1$inferity_type)

attach(data1)
convars1<-c("age","BMI","AMH","FSH","LH","E2","TT","PRL","TSH","GLU","IN","AFC","inferity")
cate<-c("type","inferity_type")
var<-c(convars1,cate)
a<-CreateTableOne(vars =var,strata="group",data1,factorVars = cate)
b<-print(a, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)
write.csv(b,"C:/Users/table1.csv")

###Table 2
####PPS
library(tableone)
data1$endogroup <- ifelse(is.na(data1$Endo_th), "NA",
                          as.character(cut(data1$Endo_th,
                                           breaks = c(-Inf, 5.999, 7.001, Inf),
                                           labels = c("<6", "6-7", ">7"),
                                           right = FALSE,
                                           include.lowest = TRUE)))
data1$endogroup <- factor(data1$endogroup)

dataex<-subset(data1,data1$exclu==0)

attach(dataex)
convars1<-c("Gn_day1","Gn_dose1","oocytes","MII","FOR","Endo_th",
            "emb_num","high_embr")
cate<-c("endogroup")
var<-c(convars1,cate)
a<-CreateTableOne(vars =var,strata="group",dataex,factorVars = cate)
b<-print(a, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)
write.csv(b,"C:/Users/table2.csv")

library(stats)
##Proportion of rest embryos-%(n/total no.)
aa<-matrix(c(70,198,121,273),nrow=2)
chisq.test(aa)

###Endometrium 6 to 7 mm on trigger day -%
bb<-matrix(c(26,127,10,140),nrow=2)
chisq.test(bb)

###Endometrium <6 mm on trigger day -%
cc<-matrix(c(7,146,6,144),nrow=2)
fisher.test(cc)

##Fresh embryo transfer rate-%
datatra<-subset(data1,data1$emb_num!=0)
attach(datatra)
cate<-c("fre_trans")
a1<-CreateTableOne(vars =cate,strata="group",datatra,factorVars = cate)
b1<-print(a1, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)
write.csv(b1,"C:/Users/table2-1.csv")


###No. of embryos transferred
library(dplyr)
library(tidyr)

dataemb <- bind_rows(
  dataex %>%
    filter(fre_trans == 1) %>%
    select(group, No_embr_fre) %>%
    rename(No_embr_tran = No_embr_fre) %>%
    mutate(transfer_type = "fresh"),
  
  dataex %>%
    filter(frozen_f == 1) %>%
    select(group, No_embr_fro_f) %>%
    rename(No_embr_tran = No_embr_fro_f) %>%
    mutate(transfer_type = "frozen"),
  
  dataex %>%
    filter(frozen_s == 1) %>%
    select(group, No_embr_fro_s) %>%
    rename(No_embr_tran = No_embr_fro_s) %>%
    mutate(transfer_type = "frozen"),
  
  dataex %>%
    filter( (is.na(fre_trans) | fre_trans != 1) & 
              (is.na(frozen_f) | frozen_f != 1) ) %>%
    select(group, No_embr_fre) %>%
    rename(No_embr_tran = No_embr_fre) %>%
    mutate(transfer_type = "none")
)

attach(dataemb)
convars1<-c("No_embr_tran")
a2<-CreateTableOne(vars =convars1,strata="group",dataemb)
b2<-print(a2, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)
write.csv(b2,"C:/Users/table2-2.csv")

##########table3  Table S7##################################
#Cumulative live birth rate -%
#Cumulative clinical pregnancy rate -%

data1$LBR<-as.factor(data1$LBR)
data1$CCP<-as.factor(data1$CCP)

attach(data1)
cate<-c("LBR","CCP")
a3<-CreateTableOne(vars =cate,strata="group",data1,factorVars = cate)
b3<-print(a3, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate

  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])

  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 1),  
    Group1_Rate = round(group1_rate * 100, 1),  
    Risk_Ratio = round(rr, 2),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

lbr_result <- calculate_risk_measures(data1, "group", "LBR")
ccp_result <- calculate_risk_measures(data1, "group", "CCP")
final_results <- rbind(lbr_result, ccp_result)
print(final_results)

library(Matrix)
library(lme4)
library(dplyr)
library(performance)

data1$center<-as.factor(data1$center)
attach(data1)
####LBR
model <- glmer(LBR ~ group + age + (1 | center), 
               family = binomial, data = data1)
summary(model)
fixed_effects <- as.data.frame(coef(summary(model)))
fixed_effects$OR <- exp(fixed_effects$Estimate)
fixed_effects$OR_CI <- sprintf(
  "%.2f (%.2f–%.2f)", 
  fixed_effects$OR,
  exp(fixed_effects$Estimate - 1.96 * fixed_effects$`Std. Error`),
  exp(fixed_effects$Estimate + 1.96 * fixed_effects$`Std. Error`)
)
print(fixed_effects[, c("Estimate", "OR", "OR_CI", "Pr(>|z|)")])
icc(model)


####CPR
model <- glmer(CPR ~ group + age + (1 | center), 
               family = binomial, data = data1)
summary(model)
fixed_effects <- as.data.frame(coef(summary(model)))
fixed_effects$OR <- exp(fixed_effects$Estimate)
fixed_effects$OR_CI <- sprintf(
  "%.2f (%.2f–%.2f)", 
  fixed_effects$OR,
  exp(fixed_effects$Estimate - 1.96 * fixed_effects$`Std. Error`),
  exp(fixed_effects$Estimate + 1.96 * fixed_effects$`Std. Error`)
)

print(fixed_effects[, c("Estimate", "OR", "OR_CI", "Pr(>|z|)")])
icc(model)

# Benjamini-Hochberg校正
p_values <- c(0.991,0.283)
adjusted_p <- p.adjust(p_values, method = "BH")  
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)


#########PPS
#Cumulative live birth rate -%
#Cumulative clinical pregnancy rate -%
dataex<-subset(data1,data1$exclu==0)

attach(dataex)
cate<-c("LBR","CCP")
a4<-CreateTableOne(vars =cate,strata="group",dataex,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 1),  
    Group1_Rate = round(group1_rate * 100, 1),  
    Risk_Ratio = round(rr, 2),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

lbr_result <- calculate_risk_measures(dataex, "group", "LBR")
ccp_result <- calculate_risk_measures(dataex, "group", "CCP")
final_results <- rbind(lbr_result, ccp_result)
print(final_results)

attach(dataex)
dataex$center<-as.factor(dataex$center)
####LBR
model <- glmer(LBR ~ group + age + (1 | center), 
               family = binomial, data = dataex)
summary(model)
fixed_effects <- as.data.frame(coef(summary(model)))
fixed_effects$OR <- exp(fixed_effects$Estimate)
fixed_effects$OR_CI <- sprintf(
  "%.2f (%.2f–%.2f)", 
  fixed_effects$OR,
  exp(fixed_effects$Estimate - 1.96 * fixed_effects$`Std. Error`),
  exp(fixed_effects$Estimate + 1.96 * fixed_effects$`Std. Error`)
)
print(fixed_effects[, c("Estimate", "OR", "OR_CI", "Pr(>|z|)")])
icc(model)

####CPR
model <- glmer(CPR ~ group + age + (1 | center), 
               family = binomial, data = dataex)
summary(model)
fixed_effects <- as.data.frame(coef(summary(model)))
fixed_effects$OR <- exp(fixed_effects$Estimate)
fixed_effects$OR_CI <- sprintf(
  "%.2f (%.2f–%.2f)", 
  fixed_effects$OR,
  exp(fixed_effects$Estimate - 1.96 * fixed_effects$`Std. Error`),
  exp(fixed_effects$Estimate + 1.96 * fixed_effects$`Std. Error`)
)

print(fixed_effects[, c("Estimate", "OR", "OR_CI", "Pr(>|z|)")])
icc(model)

#Cumulative live birth rate per embryo-available COS cycle -%
#Cumulative clinical pregnancy rate per embryo-available COS cycle -%
datacos<-subset(dataex,!is.na(dataex$No_trans))
attach(datacos)
cate<-c("LBR","CCP")
a4<-CreateTableOne(vars =cate,strata="group",datacos,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 1),  
    Group1_Rate = round(group1_rate * 100, 1),  
    Risk_Ratio = round(rr, 2),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

lbr_result <- calculate_risk_measures(datacos, "group", "LBR")
ccp_result <- calculate_risk_measures(datacos, "group", "CCP")
final_results <- rbind(lbr_result, ccp_result)
print(final_results)


####LBR
dataET<-subset(data1,data1$emb_num!=0)
model <- glmer(LBR ~ group + age + (1 | center), 
               family = binomial, data = dataET)
summary(model)
fixed_effects <- as.data.frame(coef(summary(model)))
fixed_effects$OR <- exp(fixed_effects$Estimate)
fixed_effects$OR_CI <- sprintf(
  "%.2f (%.2f–%.2f)", 
  fixed_effects$OR,
  exp(fixed_effects$Estimate - 1.96 * fixed_effects$`Std. Error`),
  exp(fixed_effects$Estimate + 1.96 * fixed_effects$`Std. Error`)
)
print(fixed_effects[, c("Estimate", "OR", "OR_CI", "Pr(>|z|)")])
icc(model)


####CPR
model <- glmer(CPR ~ group + age + (1 | center), 
               family = binomial, data = dataET)
summary(model)
fixed_effects <- as.data.frame(coef(summary(model)))
fixed_effects$OR <- exp(fixed_effects$Estimate)
fixed_effects$OR_CI <- sprintf(
  "%.2f (%.2f–%.2f)", 
  fixed_effects$OR,
  exp(fixed_effects$Estimate - 1.96 * fixed_effects$`Std. Error`),
  exp(fixed_effects$Estimate + 1.96 * fixed_effects$`Std. Error`)
)

print(fixed_effects[, c("Estimate", "OR", "OR_CI", "Pr(>|z|)")])
icc(model)

# Benjamini-Hochberg
p_values <- c(0.991,0.275,0.508,0.794)
adjusted_p <- p.adjust(p_values, method = "BH") 
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)



#####fresh ET
####Live birth rate in fresh ET-% 
####Clinical pregnancy rate in fresh ET-% 
#data3<-read_xlsx("/data/RCT_fresh_data.xlsx")
data3<-RCT_fresh_data
attach(data3)
cate<-c("LBR","CPR")
a4<-CreateTableOne(vars =cate,strata="group",data3,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 1),  
    Group1_Rate = round(group1_rate * 100, 1),  
    Risk_Ratio = round(rr, 2),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

lbr_result <- calculate_risk_measures(data3, "group", "LBR")
ccp_result <- calculate_risk_measures(data3, "group", "CPR")
final_results <- rbind(lbr_result, ccp_result)
print(final_results)


####LBR
data3<-subset(data1,data1$fre_trans==1)
model <- glmer(LBR ~ group + age + (1 | center), 
               family = binomial, data = data3)
summary(model)
fixed_effects <- as.data.frame(coef(summary(model)))
fixed_effects$OR <- exp(fixed_effects$Estimate)
fixed_effects$OR_CI <- sprintf(
  "%.2f (%.2f–%.2f)", 
  fixed_effects$OR,
  exp(fixed_effects$Estimate - 1.96 * fixed_effects$`Std. Error`),
  exp(fixed_effects$Estimate + 1.96 * fixed_effects$`Std. Error`)
)
print(fixed_effects[, c("Estimate", "OR", "OR_CI", "Pr(>|z|)")])
icc(model)

####CPR
data3<-subset(data1,data1$fre_trans==1)
model <- glmer(CPR ~ group + age + (1 | center), 
               family = binomial, data = data3)
summary(model)
fixed_effects <- as.data.frame(coef(summary(model)))
fixed_effects$OR <- exp(fixed_effects$Estimate)
fixed_effects$OR_CI <- sprintf(
  "%.2f (%.2f–%.2f)", 
  fixed_effects$OR,
  exp(fixed_effects$Estimate - 1.96 * fixed_effects$`Std. Error`),
  exp(fixed_effects$Estimate + 1.96 * fixed_effects$`Std. Error`)
)

print(fixed_effects[, c("Estimate", "OR", "OR_CI", "Pr(>|z|)")])
icc(model)

###########Pregnancy loss rate in fresh ET-%-Among biochemical pregnancy

data3<-RCT_fresh_data
####miscarr
data4<-subset(data3,data3$BIO==1)
attach(data4)
cate<-c("miscarr1")
a4<-CreateTableOne(vars =cate,strata="group",data4,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 1),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}
mis_result <- calculate_risk_measures(data4, "group", "miscarr1")
final_results <- rbind(mis_result)
print(final_results)

####miscarr
data3<-RCT_fresh_data
data4<-subset(data3,data3$BIO==1)
model <- glmer(miscarr1 ~ group + age + (1 | center), 
               family = binomial, data = data4)
summary(model)
fixed_effects <- as.data.frame(coef(summary(model)))
fixed_effects$OR <- exp(fixed_effects$Estimate)
fixed_effects$OR_CI <- sprintf(
  "%.2f (%.2f–%.2f)", 
  fixed_effects$OR,
  exp(fixed_effects$Estimate - 1.96 * fixed_effects$`Std. Error`),
  exp(fixed_effects$Estimate + 1.96 * fixed_effects$`Std. Error`)
)
print(fixed_effects[, c("Estimate", "OR", "OR_CI", "Pr(>|z|)")])
icc(model)

###########Pregnancy loss rate in fresh ET-%-Among clinical pregnancy
####miscarr
data3<-RCT_fresh_data
data5<-subset(data3,data3$CPR==1)
attach(data5)
cate<-c("miscarr1")
a4<-CreateTableOne(vars =cate,strata="group",data5,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 1),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}
mis_result <- calculate_risk_measures(data5, "group", "miscarr1")
final_results <- rbind(mis_result)
print(final_results)


data3<-RCT_fresh_data
data5<-subset(data3,data3$CPR==1)
model <- glmer(miscarr ~ group + age + (1 | center), 
               family = binomial, data = data5)
summary(model)
fixed_effects <- as.data.frame(coef(summary(model)))
fixed_effects$OR <- exp(fixed_effects$Estimate)
fixed_effects$OR_CI <- sprintf(
  "%.2f (%.2f–%.2f)", 
  fixed_effects$OR,
  exp(fixed_effects$Estimate - 1.96 * fixed_effects$`Std. Error`),
  exp(fixed_effects$Estimate + 1.96 * fixed_effects$`Std. Error`)
)
print(fixed_effects[, c("Estimate", "OR", "OR_CI", "Pr(>|z|)")])
icc(model)

# Benjamini-Hochberg
p_values <- c(0.164,0.753,0.244,0.665)
adjusted_p <- p.adjust(p_values, method = "BH") 
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)




#####table4 
##DOR
data3<-RCT_fresh_data
datador<-subset(data3,data3$inferity_type=="Diminished ovarian reserve")
attach(datador)
cate<-c("LBR","CPR")
a4<-CreateTableOne(vars =cate,strata="group",datador,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

lbr_result <- calculate_risk_measures(datador, "group", "LBR")
ccp_result <- calculate_risk_measures(datador, "group", "CPR")
final_results <- rbind(lbr_result, ccp_result)
print(final_results)




#######Clinical pregnancy rate in fresh ET-% -1 cleavage embryo transferred

datacle<-subset(datador,datador$`No. of embryos transferred in fresh transfer`==1)
attach(datacle)
cate<-c("CPR")
a4<-CreateTableOne(vars =cate,strata="group",datacle,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

ccp_result <- calculate_risk_measures(datacle, "group", "CPR")
final_results <- rbind(ccp_result)
print(final_results)


#######Clinical pregnancy rate in fresh ET-% -2 cleavage embryos transferred

datacle1<-subset(datador,datador$`No. of embryos transferred in fresh transfer`==2&
                   datador$`Type of embryos of  fresh transfer`=="cleavage embryo transferred")
attach(datacle1)
cate<-c("CPR")
a4<-CreateTableOne(vars =cate,strata="group",datacle1,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

ccp_result <- calculate_risk_measures(datacle1, "group", "CPR")
final_results <- rbind(ccp_result)
print(final_results)

#####Pregnancy loss rate among clinical pregnancy in fresh ET-%
datadorcpr<-subset(datador,datador$CPR==1)
attach(datadorcpr)
cate<-c("miscarr")
a4<-CreateTableOne(vars =cate,strata="group",datadorcpr,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

ccp_result <- calculate_risk_measures(datadorcpr, "group", "miscarr")
final_results <- rbind(ccp_result)
print(final_results)


###################Age with 40-45 years
data3<-RCT_fresh_data
dataold<-subset(data3,data3$inferity_type=="Old")
attach(dataold)
cate<-c("LBR","CPR")
a4<-CreateTableOne(vars =cate,strata="group",dataold,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

###Live birth rate in fresh ET-%-----P value
library(stats)
a1<-matrix(c(9,34,6,38),nrow=2)
fisher.test(a1)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

lbr_result <- calculate_risk_measures(dataold, "group", "LBR")
ccp_result <- calculate_risk_measures(dataold, "group", "CPR")
final_results <- rbind(lbr_result, ccp_result)
print(final_results)


#######Clinical pregnancy rate in fresh ET-% -1 cleavage embryo transferred
dataoldcle<-subset(dataold,dataold$`No. of embryos transferred in fresh transfer`==1&
                     dataold$`Type of embryos of  fresh transfer`=="cleavage embryo transferred")
attach(dataoldcle)
cate<-c("CPR")
a4<-CreateTableOne(vars =cate,strata="group",dataoldcle,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

##Clinical pregnancy rate in fresh ET-% ---1 cleavage embryo transferred-----P value
library(stats)
a1<-matrix(c(1,9,3,11),nrow=2)
fisher.test(a1)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

ccp_result <- calculate_risk_measures(dataoldcle, "group", "CPR")
final_results <- rbind(ccp_result)
print(final_results)


#######Clinical pregnancy rate in fresh ET-% -2 cleavage embryos transferred

dataoldcle1<-subset(dataold,dataold$`No. of embryos transferred in fresh transfer`==2&
                      dataold$`Type of embryos of  fresh transfer`=="cleavage embryo transferred")
attach(dataoldcle1)
cate<-c("CPR")
a4<-CreateTableOne(vars =cate,strata="group",dataoldcle1,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

ccp_result <- calculate_risk_measures(dataoldcle1, "group", "CPR")
final_results <- rbind(ccp_result)
print(final_results)


#####Pregnancy loss rate among clinical pregnancy in fresh ET-%
dataoldcpr<-subset(dataold,dataold$CPR==1)
attach(dataoldcpr)
cate<-c("miscarr")
a4<-CreateTableOne(vars =cate,strata="group",dataoldcpr,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

###Pregnancy loss rate among clinical pregnancy in fresh ET-%-----P value
library(stats)
a1<-matrix(c(3,9,5,6),nrow=2)
fisher.test(a1)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

ccp_result <- calculate_risk_measures(dataoldcpr, "group", "miscarr")
final_results <- rbind(ccp_result)
print(final_results)


# Benjamini-Hochberg校正
p_values <- c(0.132,0.056,1,0.015,1,0.408,0.949,0.615,0.64,0.4)
adjusted_p <- p.adjust(p_values, method = "BH") 
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)


######Supplemental data table S4##########################
data1$fre_trans[is.na(data1$fre_trans)] <- "NA" 
data1$fre_trans<-as.factor(data1$fre_trans)
attach(data1)
convars1<-c("Gn_day","Gn_dose","oocytes","MII","FOR","Endo_th",
            "emb_num","high_embr")
cate<-c("fre_trans")
var<-c(convars1,cate)
a<-CreateTableOne(vars =var,strata="group",data1,factorVars = cate)
b<-print(a, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)
write.csv(b,"C:/Users/tableS4.csv")

library(stats)
aa<-matrix(c(94,65,99,60),nrow=2)
chisq.test(aa)

###No. of embryos transferred
library(dplyr)
library(tidyr)

dataemb <- bind_rows(
  dataex %>%
    filter(fre_trans == 1) %>%
    select(group, No_embr_fre) %>%
    rename(No_embr_tran = No_embr_fre) %>%
    mutate(transfer_type = "fresh"),
  
  dataex %>%
    filter(frozen_f == 1) %>%
    select(group, No_embr_fro_f) %>%
    rename(No_embr_tran = No_embr_fro_f) %>%
    mutate(transfer_type = "frozen"),
  
  dataex %>%
    filter(frozen_s == 1) %>%
    select(group, No_embr_fro_s) %>%
    rename(No_embr_tran = No_embr_fro_s) %>%
    mutate(transfer_type = "frozen"),
  
  dataex %>%
    filter( (is.na(fre_trans) | fre_trans != 1) & 
              (is.na(frozen_f) | frozen_f != 1) ) %>%
    select(group, No_embr_fre) %>%
    rename(No_embr_tran = No_embr_fre) %>%
    mutate(transfer_type = "none")
)

attach(dataemb)
convars1<-c("No_embr_tran")
a2<-CreateTableOne(vars =convars1,strata="group",dataemb)
b2<-print(a2, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)
write.csv(b2,"C:/Users/tableS4-2.csv")


##############Table S5#############
library(stats)
a1<-matrix(c(7,146,1,149),nrow=2)
fisher.test(a1)

a2<-matrix(c(52,94,50,99),nrow=2)
chisq.test(a2)

a3<-matrix(c(9,43,5,45),nrow=2)
fisher.test(a3)

a4<-matrix(c(16,36,6,44),nrow=2)
chisq.test(a4)

a5<-matrix(c(4,48,7,43),nrow=2)
fisher.test(a5)

a6<-matrix(c(10,42,12,38),nrow=2)
chisq.test(a6)

a7<-matrix(c(2,50,2,48),nrow=2)
fisher.test(a7)

a8<-matrix(c(4,48,6,44),nrow=2)
fisher.test(a8)

a9<-matrix(c(1,51,4,46),nrow=2)
fisher.test(a9)

a10<-matrix(c(6,46,8,42),nrow=2)
chisq.test(a10)

##############Table S6#############
attach(dataex)
cate<-c("IVF_type")
a<-CreateTableOne(vars =cate,strata="group",dataex,factorVars = cate)
b<-print(a, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

library(stats)
aa<-matrix(c(109,25,3,119,21,4),nrow=3)
fisher.test(aa)


#############Table S7 see code of Table 3#######


#############Table S8   center correction#################
#LBR CLBR per COS cycle-%
data <- data.frame(
  center = rep(1:6, each = 2),
  treatment = rep(c("mLP", "GnRH-ant"), 6),
  responder = c(6,6,10,14,2,0,1,0,18,16,2,0),
  total = c(16,18,60,48,11,12,12,9,54,70,6,2)
)
fit <- glmer(cbind(responder, total - responder) ~ treatment + (1 | center),
             family = binomial, data = data)

fe <- summary(fit)$coefficients
fe_df <- data.frame(
  term = rownames(fe),
  estimate = fe[, "Estimate"],
  std.error = fe[, "Std. Error"],
  statistic = fe[, "z value"],
  p.value = fe[, "Pr(>|z|)"],
  OR = exp(fe[, "Estimate"]),
  OR_CI_lower = exp(fe[, "Estimate"] - 1.96 * fe[, "Std. Error"]),
  OR_CI_upper = exp(fe[, "Estimate"] + 1.96 * fe[, "Std. Error"])
)

re <- as.data.frame(VarCorr(fit))

icc_result <- icc(fit)  
vc <- as.data.frame(VarCorr(fit))
var_random <- vc$vcov[1]  
se_random <- vc$sdcor[1]   
var_residual <- pi^2 / 3 
icc_value <- var_random / (var_random + var_residual)
se_icc <- se_random * var_residual / (var_random + var_residual)^2
icc_ci_lower <- max(0, icc_value - 1.96 * se_icc)
icc_ci_upper <- min(1, icc_value + 1.96 * se_icc)

cat("===================== ICC (Adjusted) =====================\n")
cat(sprintf(
  "ICC = %.3f (95%% CI: [%.3f, %.3f])\n",
  icc_value,
  icc_ci_lower,
  icc_ci_upper
))

cat("\n对比 performance::icc() 的结果：\n")
print(icc_result)
cat("===================== Fixed Effects =====================\n")
print(fe_df)
cat("\n=================== Random Effects ===================\n")
print(re)
cat("\n=================== ICC with 95% CI ===================\n")
print(icc_result)


#CPR  CCPR per COS cycle 
data <- data.frame(
  center = rep(1:6, each = 2),
  treatment = rep(c("mLP", "GnRH-ant"), 6),
  responder = c(6,11,13,17,2,1,4,1,22,24,4,0),
  total = c(16,18,60,48,11,12,12,9,54,70,6,2)
)

fit <- glmer(cbind(responder, total - responder) ~ treatment + (1 | center),
             family = binomial, data = data)

fe <- summary(fit)$coefficients
fe_df <- data.frame(
  term = rownames(fe),
  estimate = fe[, "Estimate"],
  std.error = fe[, "Std. Error"],
  statistic = fe[, "z value"],
  p.value = fe[, "Pr(>|z|)"],
  OR = exp(fe[, "Estimate"]),
  OR_CI_lower = exp(fe[, "Estimate"] - 1.96 * fe[, "Std. Error"]),
  OR_CI_upper = exp(fe[, "Estimate"] + 1.96 * fe[, "Std. Error"])
)

re <- as.data.frame(VarCorr(fit))

icc_result <- icc(fit)  
vc <- as.data.frame(VarCorr(fit))
var_random <- vc$vcov[1]  
se_random <- vc$sdcor[1]  
var_residual <- pi^2 / 3   
icc_value <- var_random / (var_random + var_residual)
se_icc <- se_random * var_residual / (var_random + var_residual)^2
icc_ci_lower <- max(0, icc_value - 1.96 * se_icc)
icc_ci_upper <- min(1, icc_value + 1.96 * se_icc)

cat("===================== ICC (Adjusted) =====================\n")
cat(sprintf(
  "ICC = %.3f (95%% CI: [%.3f, %.3f])\n",
  icc_value,
  icc_ci_lower,
  icc_ci_upper
))

cat("\n对比 performance::icc() 的结果：\n")
print(icc_result)
cat("===================== Fixed Effects =====================\n")
print(fe_df)
cat("\n=================== Random Effects ===================\n")
print(re)
cat("\n=================== ICC with 95% CI ===================\n")
print(icc_result)




#CPR
data <- data.frame(
  center = rep(1:6, each = 2),
  treatment = rep(c("mLP", "GnRH-ant"), 6),
  responder = c(6,11,13,17,2,1,4,1,22,24,4,0),
  total = c(16,18,54,40,11,11,12,9,54,70,6,2)
)

fit <- glmer(cbind(responder, total - responder) ~ treatment + (1 | center),
             family = binomial, data = data)

# OR 95%CI
fe <- summary(fit)$coefficients
fe_df <- data.frame(
  term = rownames(fe),
  estimate = fe[, "Estimate"],
  std.error = fe[, "Std. Error"],
  statistic = fe[, "z value"],
  p.value = fe[, "Pr(>|z|)"],
  OR = exp(fe[, "Estimate"]),
  OR_CI_lower = exp(fe[, "Estimate"] - 1.96 * fe[, "Std. Error"]),
  OR_CI_upper = exp(fe[, "Estimate"] + 1.96 * fe[, "Std. Error"])
)

re <- as.data.frame(VarCorr(fit))
# ICC  95%CI
icc_result <- icc(fit) 

vc <- as.data.frame(VarCorr(fit))
var_random <- vc$vcov[1]  
se_random <- vc$sdcor[1]   
var_residual <- pi^2 / 3  
icc_value <- var_random / (var_random + var_residual)
se_icc <- se_random * var_residual / (var_random + var_residual)^2

icc_ci_lower <- max(0, icc_value - 1.96 * se_icc)
icc_ci_upper <- min(1, icc_value + 1.96 * se_icc)

cat("===================== ICC (Adjusted) =====================\n")
cat(sprintf(
  "ICC = %.3f (95%% CI: [%.3f, %.3f])\n",
  icc_value,
  icc_ci_lower,
  icc_ci_upper
))

cat("\n对比 performance::icc() 的结果：\n")
print(icc_result)
cat("===================== Fixed Effects =====================\n")
print(fe_df)
cat("\n=================== Random Effects ===================\n")
print(re)
cat("\n=================== ICC with 95% CI ===================\n")
print(icc_result)



#IR
data <- data.frame(
  center = rep(1:6, each = 2),
  treatment = rep(c("mLP", "GnRH-ant"), 6),
  responder = c(6,11,16,18,2,1,4,1,25,25,4,0),
  total = c(16,18,54,40,11,11,12,9,54,70,6,2)
)

fit <- glmer(cbind(responder, total - responder) ~ treatment + (1 | center),
             family = binomial, data = data)

#OR 95%CI
fe <- summary(fit)$coefficients
fe_df <- data.frame(
  term = rownames(fe),
  estimate = fe[, "Estimate"],
  std.error = fe[, "Std. Error"],
  statistic = fe[, "z value"],
  p.value = fe[, "Pr(>|z|)"],
  OR = exp(fe[, "Estimate"]),
  OR_CI_lower = exp(fe[, "Estimate"] - 1.96 * fe[, "Std. Error"]),
  OR_CI_upper = exp(fe[, "Estimate"] + 1.96 * fe[, "Std. Error"])
)

re <- as.data.frame(VarCorr(fit))
icc_result <- icc(fit)  
vc <- as.data.frame(VarCorr(fit))
var_random <- vc$vcov[1]  
se_random <- vc$sdcor[1]   
var_residual <- pi^2 / 3  
icc_value <- var_random / (var_random + var_residual)
se_icc <- se_random * var_residual / (var_random + var_residual)^2
icc_ci_lower <- max(0, icc_value - 1.96 * se_icc)
icc_ci_upper <- min(1, icc_value + 1.96 * se_icc)

cat("===================== ICC (Adjusted) =====================\n")
cat(sprintf(
  "ICC = %.3f (95%% CI: [%.3f, %.3f])\n",
  icc_value,
  icc_ci_lower,
  icc_ci_upper
))

cat("\n对比 performance::icc() 的结果：\n")
print(icc_result)
cat("===================== Fixed Effects =====================\n")
print(fe_df)
cat("\n=================== Random Effects ===================\n")
print(re)
cat("\n=================== ICC with 95% CI ===================\n")
print(icc_result)



#LBR
data <- data.frame(
  center = rep(1:6, each = 2),
  treatment = rep(c("mLP", "GnRH-ant"), 6),
  responder = c(6,6,10,14,2,0,1,0,18,16,2,0),
  total = c(16,18,54,40,11,11,12,9,54,70,6,2)
)

fit <- glmer(cbind(responder, total - responder) ~ treatment + (1 | center),
             family = binomial, data = data)

# OR 95%CI
fe <- summary(fit)$coefficients
fe_df <- data.frame(
  term = rownames(fe),
  estimate = fe[, "Estimate"],
  std.error = fe[, "Std. Error"],
  statistic = fe[, "z value"],
  p.value = fe[, "Pr(>|z|)"],
  OR = exp(fe[, "Estimate"]),
  OR_CI_lower = exp(fe[, "Estimate"] - 1.96 * fe[, "Std. Error"]),
  OR_CI_upper = exp(fe[, "Estimate"] + 1.96 * fe[, "Std. Error"])
)

re <- as.data.frame(VarCorr(fit))

icc_result <- icc(fit) 
vc <- as.data.frame(VarCorr(fit))
var_random <- vc$vcov[1]  
se_random <- vc$sdcor[1]   
var_residual <- pi^2 / 3   
icc_value <- var_random / (var_random + var_residual)
se_icc <- se_random * var_residual / (var_random + var_residual)^2

icc_ci_lower <- max(0, icc_value - 1.96 * se_icc)
icc_ci_upper <- min(1, icc_value + 1.96 * se_icc)

cat("===================== ICC (Adjusted) =====================\n")
cat(sprintf(
  "ICC = %.3f (95%% CI: [%.3f, %.3f])\n",
  icc_value,
  icc_ci_lower,
  icc_ci_upper
))

cat("\n对比 performance::icc() 的结果：\n")
print(icc_result)
cat("===================== Fixed Effects =====================\n")
print(fe_df)

cat("\n=================== Random Effects ===================\n")
print(re)

cat("\n=================== ICC with 95% CI ===================\n")
print(icc_result)




#IR
data <- data.frame(
  center = rep(1:6, each = 2),
  treatment = rep(c("mLP", "GnRH-ant"), 6),
  responder = c(6,11,16,18,2,1,4,1,25,25,4,0),
  total = c(16,18,60,48,11,12,12,9,54,70,6,2)
)
fit <- glmer(cbind(responder, total - responder) ~ treatment + (1 | center),
             family = binomial, data = data)

fe <- summary(fit)$coefficients
fe_df <- data.frame(
  term = rownames(fe),
  estimate = fe[, "Estimate"],
  std.error = fe[, "Std. Error"],
  statistic = fe[, "z value"],
  p.value = fe[, "Pr(>|z|)"],
  OR = exp(fe[, "Estimate"]),
  OR_CI_lower = exp(fe[, "Estimate"] - 1.96 * fe[, "Std. Error"]),
  OR_CI_upper = exp(fe[, "Estimate"] + 1.96 * fe[, "Std. Error"])
)

re <- as.data.frame(VarCorr(fit))

icc_result <- icc(fit)  

vc <- as.data.frame(VarCorr(fit))
var_random <- vc$vcov[1]  
se_random <- vc$sdcor[1]   
var_residual <- pi^2 / 3   
icc_value <- var_random / (var_random + var_residual)
se_icc <- se_random * var_residual / (var_random + var_residual)^2
icc_ci_lower <- max(0, icc_value - 1.96 * se_icc)
icc_ci_upper <- min(1, icc_value + 1.96 * se_icc)

cat("===================== ICC (Adjusted) =====================\n")
cat(sprintf(
  "ICC = %.3f (95%% CI: [%.3f, %.3f])\n",
  icc_value,
  icc_ci_lower,
  icc_ci_upper
))

cat("\n对比 performance::icc() 的结果：\n")
print(icc_result)
cat("===================== Fixed Effects =====================\n")
print(fe_df)
cat("\n=================== Random Effects ===================\n")
print(re)
cat("\n=================== ICC with 95% CI ===================\n")
print(icc_result)







###########################Table S9#############
########IVF cycle =1
data3<-RCT_fresh_data
datacycfre<-subset(data3,data3$cycle==1)
attach(datacycfre)
cate<-c("LBR","CPR")
a3<-CreateTableOne(vars =cate,strata="group",datacycfre,factorVars = cate)
b3<-print(a3, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

lbr_result <- calculate_risk_measures(datacycfre, "group", "LBR")
ccp_result <- calculate_risk_measures(datacycfre, "group", "CPR")
final_results <- rbind(lbr_result, ccp_result)
print(final_results)

#######Clinical pregnancy rate in fresh ET-% -1 cleavage embryo transferred

datacyccle<-subset(datacycfre,datacycfre$`No. of embryos transferred in fresh transfer`==1&
                     datacycfre$`Type of embryos of  fresh transfer`=="cleavage embryo transferred")
attach(datacyccle)
cate<-c("CPR")
a4<-CreateTableOne(vars =cate,strata="group",datacyccle,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

ccp_result <- calculate_risk_measures(datacyccle, "group", "CPR")
final_results <- rbind(ccp_result)
print(final_results)


#######Clinical pregnancy rate in fresh ET-% -2 cleavage embryos transferred

datacyccle1<-subset(datacycfre,datacycfre$`No. of embryos transferred in fresh transfer`==2&
                      datacycfre$`Type of embryos of  fresh transfer`=="cleavage embryo transferred")
attach(datacyccle1)
cate<-c("CPR")
a4<-CreateTableOne(vars =cate,strata="group",datacyccle1,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

ccp_result <- calculate_risk_measures(datacyccle1, "group", "CPR")
final_results <- rbind(ccp_result)
print(final_results)


########IVF cycle > 1
data3<-RCT_fresh_data
datacycfre1<-subset(data3,!data3$cycle==1)
attach(datacycfre1)
cate<-c("LBR","CPR")
a3<-CreateTableOne(vars =cate,strata="group",datacycfre1,factorVars = cate)
b3<-print(a3, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

lbr_result <- calculate_risk_measures(datacycfre1, "group", "LBR")
ccp_result <- calculate_risk_measures(datacycfre1, "group", "CPR")
final_results <- rbind(lbr_result, ccp_result)
print(final_results)

#######Clinical pregnancy rate in fresh ET-% -1 cleavage embryo transferred

datacyc1cle1<-subset(datacycfre1,datacycfre1$`No. of embryos transferred in fresh transfer`==1&
                     datacycfre1$`Type of embryos of  fresh transfer`=="cleavage embryo transferred")
attach(datacyc1cle1)
cate<-c("CPR")
a4<-CreateTableOne(vars =cate,strata="group",datacyc1cle1,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

####P-value
aa<-matrix(c(2,4,1,9),nrow=2)
fisher.test(aa)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

ccp_result <- calculate_risk_measures(datacyc1cle1, "group", "CPR")
final_results <- rbind(ccp_result)
print(final_results)


#######Clinical pregnancy rate in fresh ET-% -2 cleavage embryos transferred

datacyc1cle2<-subset(datacycfre1,datacycfre1$`No. of embryos transferred in fresh transfer`==2&
                      datacycfre1$`Type of embryos of  fresh transfer`=="cleavage embryo transferred")
attach(datacyc1cle2)
cate<-c("CPR")
a4<-CreateTableOne(vars =cate,strata="group",datacyc1cle2,factorVars = cate)
b4<-print(a4, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)

calculate_risk_measures <- function(data, group_var, outcome_var) {
  tab <- table(data[[group_var]], data[[outcome_var]])
  prop_tab <- prop.table(tab, margin = 1)
  
  group0_rate <- prop_tab[1, "1"]
  group1_rate <- prop_tab[2, "1"]
  rr <- group1_rate / group0_rate
  rd <- group1_rate - group0_rate
  
  n0 <- sum(tab[1, ])
  n1 <- sum(tab[2, ])
  
  se_log_rr <- sqrt((1/tab[1, "1"] - 1/n0) + (1/tab[2, "1"] - 1/n1))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  se_rd <- sqrt(group0_rate*(1-group0_rate)/n0 + group1_rate*(1-group1_rate)/n1)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100  
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100  
  result <- data.frame(
    Outcome = outcome_var,
    Group0_Rate = round(group0_rate * 100, 2),  
    Group1_Rate = round(group1_rate * 100, 2),  
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0("(", round(rr_ci_lower, 2), ", ", round(rr_ci_upper, 2), ")"),
    Risk_Difference = round(rd * 100, 2), 
    RD_CI = paste0("(", round(rd_ci_lower, 2), ", ", round(rd_ci_upper, 2), ")")  
  )
  return(result)
}

ccp_result <- calculate_risk_measures(datacyc1cle2, "group", "CPR")
final_results <- rbind(ccp_result)
print(final_results)


#CLBR per embryo-available COS cycle-% IVF cycle =1
calculate_risk_from_table <- function(treatment_success, treatment_total, control_success, control_total, outcome_name = "Outcome") {
  tab <- matrix(c(treatment_success, treatment_total - treatment_success,
                  control_success, control_total - control_success), 
                nrow = 2, byrow = TRUE)

  treatment_rate <- treatment_success / treatment_total
  control_rate <- control_success / control_total

  rr <- treatment_rate / control_rate
  rd <- treatment_rate - control_rate

  n_treatment <- treatment_total
  n_control <- control_total

  se_log_rr <- sqrt((1/treatment_success - 1/n_treatment) + (1/control_success - 1/n_control))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)

  se_rd <- sqrt(control_rate*(1-control_rate)/n_control + treatment_rate*(1-treatment_rate)/n_treatment)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100
  
  chisq_test <- chisq.test(tab)
  fisher_test <- fisher.test(tab)
  
  result <- data.frame(
    Outcome = outcome_name,
    Treatment_Rate = paste0(round(treatment_rate * 100, 2), "% (", treatment_success, "/", treatment_total, ")"),
    Control_Rate = paste0(round(control_rate * 100, 2), "% (", control_success, "/", control_total, ")"),
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0(round(rr_ci_lower, 2), " - ", round(rr_ci_upper, 2)),
    Risk_Difference = paste0(round(rd * 100, 2), "%"),
    RD_CI = paste0(round(rd_ci_lower, 2), "% to ", round(rd_ci_upper, 2), "%"),
    Chisq_pvalue = format.pval(chisq_test$p.value, digits = 3),
    Fisher_pvalue = format.pval(fisher_test$p.value, digits = 3)
  )
  
  return(result)
}
result <- calculate_risk_from_table(
  treatment_success = 22, 
  treatment_total = 73, 
  control_success = 26, 
  control_total = 90,
  outcome_name = "CPR"
)

print(result)

#CLBR per embryo-available COS cycle-% IVF cycle>1
calculate_risk_from_table <- function(treatment_success, treatment_total, control_success, control_total, outcome_name = "Outcome") {
  tab <- matrix(c(treatment_success, treatment_total - treatment_success,
                  control_success, control_total - control_success), 
                nrow = 2, byrow = TRUE)
  
  treatment_rate <- treatment_success / treatment_total
  control_rate <- control_success / control_total
  
  rr <- treatment_rate / control_rate
  rd <- treatment_rate - control_rate
  
  n_treatment <- treatment_total
  n_control <- control_total
  
  se_log_rr <- sqrt((1/treatment_success - 1/n_treatment) + (1/control_success - 1/n_control))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  
  se_rd <- sqrt(control_rate*(1-control_rate)/n_control + treatment_rate*(1-treatment_rate)/n_treatment)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100
  
  chisq_test <- chisq.test(tab)
  fisher_test <- fisher.test(tab)
  
  result <- data.frame(
    Outcome = outcome_name,
    Treatment_Rate = paste0(round(treatment_rate * 100, 2), "% (", treatment_success, "/", treatment_total, ")"),
    Control_Rate = paste0(round(control_rate * 100, 2), "% (", control_success, "/", control_total, ")"),
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0(round(rr_ci_lower, 2), " - ", round(rr_ci_upper, 2)),
    Risk_Difference = paste0(round(rd * 100, 2), "%"),
    RD_CI = paste0(round(rd_ci_lower, 2), "% to ", round(rd_ci_upper, 2), "%"),
    Chisq_pvalue = format.pval(chisq_test$p.value, digits = 3),
    Fisher_pvalue = format.pval(fisher_test$p.value, digits = 3)
  )
  
  return(result)
}
result <- calculate_risk_from_table(
  treatment_success = 17, 
  treatment_total = 39, 
  control_success = 10, 
  control_total = 41,
  outcome_name = "CPR"
)

print(result)

###################Table S10
#####Clinical pregnancy rate in fresh ET-%
calculate_risk_from_table <- function(treatment_success, treatment_total, control_success, control_total, outcome_name = "Outcome") {
  tab <- matrix(c(treatment_success, treatment_total - treatment_success,
                  control_success, control_total - control_success), 
                nrow = 2, byrow = TRUE)
  
  treatment_rate <- treatment_success / treatment_total
  control_rate <- control_success / control_total
  
  rr <- treatment_rate / control_rate
  rd <- treatment_rate - control_rate
  
  n_treatment <- treatment_total
  n_control <- control_total
  
  se_log_rr <- sqrt((1/treatment_success - 1/n_treatment) + (1/control_success - 1/n_control))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  
  se_rd <- sqrt(control_rate*(1-control_rate)/n_control + treatment_rate*(1-treatment_rate)/n_treatment)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100
  
  chisq_test <- chisq.test(tab)
  fisher_test <- fisher.test(tab)
  
  result <- data.frame(
    Outcome = outcome_name,
    Treatment_Rate = paste0(round(treatment_rate * 100, 2), "% (", treatment_success, "/", treatment_total, ")"),
    Control_Rate = paste0(round(control_rate * 100, 2), "% (", control_success, "/", control_total, ")"),
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0(round(rr_ci_lower, 2), " - ", round(rr_ci_upper, 2)),
    Risk_Difference = paste0(round(rd * 100, 2), "%"),
    RD_CI = paste0(round(rd_ci_lower, 2), "% to ", round(rd_ci_upper, 2), "%"),
    Chisq_pvalue = format.pval(chisq_test$p.value, digits = 3),
    Fisher_pvalue = format.pval(fisher_test$p.value, digits = 3)
  )
  
  return(result)
}
result <- calculate_risk_from_table(
  treatment_success = 1, 
  treatment_total = 9, 
  control_success = 0, 
  control_total = 4,
  outcome_name = "CPR"
)

print(result)

#############DOR
calculate_risk_from_table <- function(treatment_success, treatment_total, control_success, control_total, outcome_name = "Outcome") {
  tab <- matrix(c(treatment_success, treatment_total - treatment_success,
                  control_success, control_total - control_success), 
                nrow = 2, byrow = TRUE)
  
  treatment_rate <- treatment_success / treatment_total
  control_rate <- control_success / control_total
  
  rr <- treatment_rate / control_rate
  rd <- treatment_rate - control_rate
  
  n_treatment <- treatment_total
  n_control <- control_total
  
  se_log_rr <- sqrt((1/treatment_success - 1/n_treatment) + (1/control_success - 1/n_control))
  rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
  
  se_rd <- sqrt(control_rate*(1-control_rate)/n_control + treatment_rate*(1-treatment_rate)/n_treatment)
  rd_ci_lower <- (rd - 1.96 * se_rd) * 100
  rd_ci_upper <- (rd + 1.96 * se_rd) * 100
  
  chisq_test <- chisq.test(tab)
  fisher_test <- fisher.test(tab)
  
  result <- data.frame(
    Outcome = outcome_name,
    Treatment_Rate = paste0(round(treatment_rate * 100, 2), "% (", treatment_success, "/", treatment_total, ")"),
    Control_Rate = paste0(round(control_rate * 100, 2), "% (", control_success, "/", control_total, ")"),
    Risk_Ratio = round(rr, 3),
    RR_CI = paste0(round(rr_ci_lower, 2), " - ", round(rr_ci_upper, 2)),
    Risk_Difference = paste0(round(rd * 100, 2), "%"),
    RD_CI = paste0(round(rd_ci_lower, 2), "% to ", round(rd_ci_upper, 2), "%"),
    Chisq_pvalue = format.pval(chisq_test$p.value, digits = 3),
    Fisher_pvalue = format.pval(fisher_test$p.value, digits = 3)
  )
  
  return(result)
}
result <- calculate_risk_from_table(
  treatment_success = 1, 
  treatment_total = 4, 
  control_success = 0, 
  control_total = 3,
  outcome_name = "CPR"
)

print(result)






