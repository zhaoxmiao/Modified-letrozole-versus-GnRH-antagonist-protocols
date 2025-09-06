
###Table 1
# library(tableone)
if(!require(tableone))
  install.packages("tableone")
library(Matrix)
library(lme4)
library(performance) 
library(broom.mixed)
library(readr)
install.packages("readxl")


data1<-read_csv("/data/20250511_RCT_data.csv")
data1$group<-ifelse(data1$group==2,0,data1$group)

data1$group<-as.factor(data1$group)
data1$inferity_type<-as.factor(data1$inferity_type)
data1$fre_trans<-as.factor(data1$fre_trans)

data1$Gndaygroup <- cut(data1$Gn_day,
                        breaks = c(-Inf, 9, 13, Inf),
                        labels = c(1, 2, 3),
                        right = FALSE,
                        include.lowest = TRUE)

data1$Gndaygroup[is.na(data1$Gndaygroup)] <- 1
data1$Gndaygroup<-as.factor(data1$Gndaygroup)


dataGn1 <- data1[!is.na(data1$Gndaygroup) & data1$Gndaygroup == 1, ]
dataGn2 <- data1[!is.na(data1$Gndaygroup) & data1$Gndaygroup == 2, ]
dataGn3 <- data1[data1$Gndaygroup == 3, ]

attach(data1)

datacase<-data1[data1$group == 1, ]
datacon<-data1[data1$group ==0, ]

attach(datacase)
convars1<-c("age","BMI","AMH","FSH","LH","E2","TT","PRL","TSH","GLU","IN","AFC","cycle","inferity")
cate<-c("inferity_type","fre_trans")
var<-c(convars1,cate)
a<-CreateTableOne(vars =var,strata="Gndaygroup",datacase,factorVars = cate)
b<-print(a, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)
write.csv(b,"table1.csv")

convars2<-c("Gn_day","Gn_dose","oocytes","MII","FOR","Endo_th",
           "emb_num","high_embr","No_embr_fre")
a<-CreateTableOne(convars2,strata="group",data1)
b<-print(a, showAllLevels = TRUE,quote = TRUE, noSpaces = TRUE)
write.csv(b,"table1-1.csv")

####logistic  Table3
library(Matrix)
library(lme4)
library(dplyr)
library(performance)

data1$center<-as.factor(data1$center)
attach(data1)

# 拟合模型CPR
model <- glmer(CPR ~ group + age + (1 | center), 
               family = binomial, data = data1)
summary(model)
# 计算OR及95% CI
fixed_effects <- as.data.frame(coef(summary(model)))
fixed_effects$OR <- exp(fixed_effects$Estimate)
fixed_effects$OR_CI <- sprintf(
  "%.2f (%.2f–%.2f)", 
  fixed_effects$OR,
  exp(fixed_effects$Estimate - 1.96 * fixed_effects$`Std. Error`),
  exp(fixed_effects$Estimate + 1.96 * fixed_effects$`Std. Error`)
)

print(fixed_effects[, c("Estimate", "OR", "OR_CI", "Pr(>|z|)")])

# 计算ICC
icc(model)

####Bio_PR
model <- glmer(Bio_PR ~ group + age + (1 | center), 
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

####miscarr
model <- glmer(miscarr ~ group + age + (1 | center), 
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

p_values <- c(0.991,0.283)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)


####排除失访退出
dataex<-subset(data1,data1$exclu==0)
# CPR
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

####Bio_PR
model <- glmer(Bio_PR ~ group + age + (1 | center), 
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

####miscarr
model <- glmer(miscarr ~ group + age + (1 | center), 
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

p_values <- c(0.991, 0.275)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)

##新鲜周期
data2<-subset(data1,data1$fre_trans==1)
####LBR
model <- glmer(LBR ~ group + age + (1 | center), 
               family = binomial, data = data2)
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

# CPR
model <- glmer(CPR ~ group + age + (1 | center), 
               family = binomial, data = data2)
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

####miscarr
model <- glmer(miscarr ~ group + age + (1 | center), 
               family = binomial, data = data2)
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



p_values <- c(0.508, 0.794)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)

###########补Pregnancy loss rate in fresh ET-%-Among biochemical pregnancy
##########BIO=1 中 miscarr1 为1的logistics
library(readxl)
data3<-read_excel("/data/RCT_fresh_data.xlsx")

####miscarr
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

###########补Pregnancy loss rate in fresh ET-%-Among clinical pregnancy
##########CPR=1 中 miscarr 为1的logistics

####miscarr
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

p_values <- c(0.244, 0.665)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)





##交互作用
library(sandwich) 
library(lmtest)

##新鲜周期
data2<-subset(data1,data1$fre_trans==1)

# 拟合Log-Binomial回归（直接估计RR）
fit <- glm(
  LBR ~ group * factor(inferity_type), 
  family = binomial(link = "log"),  
  data = data2
)
# 计算稳健标准误（防止过离散）
coeftest(fit, vcov = vcovHC(fit, type = "HC1"))


#新鲜周期DOR
data3<-subset(data2,data2$inferity_type==2)
fit <- glm(
  CPR ~ group * factor(No_embr_fre), 
  family = binomial(link = "log"),  
  data = data3
)
coeftest(fit, vcov = vcovHC(fit, type = "HC1"))


#新鲜周期old
data4<-subset(data2,data2$inferity_type==1)
fit <- glm(
  CPR ~ group * factor(No_embr_fre), 
  family = binomial(link = "log"),  
  data = data4
)
coeftest(fit, vcov = vcovHC(fit, type = "HC1"))
 

##新鲜周期妊娠
data5<-subset(data2,data2$CPR==1)
fit <- glm(
  miscarr ~ group * factor(inferity_type), 
  family = binomial(link = "log"),  
  data = data5
)
coeftest(fit, vcov = vcovHC(fit, type = "HC1"))




#Table 3 Benjamini-Hochberg校正的P值
# 示例原始P值
p_values <- c(0.991,0.283)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
# 输出结果
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)

p_values <- c(0.991,0.275)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
# 输出结果
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)

p_values <- c(0.508,0.794)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
# 输出结果
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)

p_values <- c(0.164,0.753)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
# 输出结果
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)

p_values <- c(0.244,0.665)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
# 输出结果
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)

p_values <- c(0.132,0.056)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
# 输出结果
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)

p_values <- c(0.132,0.056,1)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
# 输出结果
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)

p_values <- c(1,0.015)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
# 输出结果
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)

p_values <- c(0.408,0.949,0.4)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
# 输出结果
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)

p_values <- c(0.615,0.64)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
# 输出结果
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)


#####补数据：Cumulative live birth rate per woman following one ET
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

# CPR
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

####miscarr
model <- glmer(miscarr ~ group + age + (1 | center), 
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

####Bio_PR
model <- glmer(Bio_PR ~ group + age + (1 | center), 
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

p_values <- c(0.508,0.194, 0.794, 0.109)
# Benjamini-Hochberg校正
adjusted_p <- p.adjust(p_values, method = "BH")  # 或 method = "fdr"（二者等价）
# 输出结果
data.frame(
  Original_P = p_values,
  Adjusted_P_BH = adjusted_p
)

#############center correction
#CPR
# 构建数据框
data <- data.frame(
  center = rep(1:6, each = 2),
  treatment = rep(c("mLP", "GnRH-ant"), 6),
  responder = c(6,11,13,17,2,1,4,1,22,24,4,0),
  total = c(16,18,54,40,11,11,12,9,54,70,6,2)
)

# 拟合GLMM模型
fit <- glmer(cbind(responder, total - responder) ~ treatment + (1 | center),
             family = binomial, data = data)

# --------------- 输出结果 ---------------
# 1. 固定效应（含OR和95%CI）
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

# 2. 随机效应方差
re <- as.data.frame(VarCorr(fit))

# 3. 计算ICC及其95%CI
icc_result <- icc(fit)  # 计算ICC和CI
# 3. 手动计算 ICC 的 95% CI（基于 Delta 方法 + 逻辑回归残差方差 ≈ π²/3）
# ------------------------------------------
# 提取随机效应方差
vc <- as.data.frame(VarCorr(fit))
var_random <- vc$vcov[1]  # 随机截距方差
se_random <- vc$sdcor[1]   # 随机截距的标准误

# 逻辑回归的残差方差（固定值）
var_residual <- pi^2 / 3   # ≈ 3.289868

# 计算 ICC（调整后的）
icc_value <- var_random / (var_random + var_residual)

# 计算 ICC 的标准误（Delta 方法）
se_icc <- se_random * var_residual / (var_random + var_residual)^2

# 计算 95% 置信区间（并限制在 [0, 1] 范围内）
icc_ci_lower <- max(0, icc_value - 1.96 * se_icc)
icc_ci_upper <- min(1, icc_value + 1.96 * se_icc)
# ------------------------------------------

# 4. 输出结果
cat("===================== ICC (Adjusted) =====================\n")
cat(sprintf(
  "ICC = %.3f (95%% CI: [%.3f, %.3f])\n",
  icc_value,
  icc_ci_lower,
  icc_ci_upper
))

# 对比 performance::icc() 的点估计值（验证一致性）
cat("\n对比 performance::icc() 的结果：\n")
print(icc_result)
# --------------- 打印结果 ---------------
cat("===================== Fixed Effects =====================\n")
print(fe_df)

cat("\n=================== Random Effects ===================\n")
print(re)

cat("\n=================== ICC with 95% CI ===================\n")
print(icc_result)



#IR
# 构建数据框
data <- data.frame(
  center = rep(1:6, each = 2),
  treatment = rep(c("mLP", "GnRH-ant"), 6),
  responder = c(6,11,16,18,2,1,4,1,25,25,4,0),
  total = c(16,18,54,40,11,11,12,9,54,70,6,2)
)

# 拟合GLMM模型
fit <- glmer(cbind(responder, total - responder) ~ treatment + (1 | center),
             family = binomial, data = data)

# --------------- 输出结果 ---------------
# 1. 固定效应（含OR和95%CI）
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

# 2. 随机效应方差
re <- as.data.frame(VarCorr(fit))

# 3. 计算ICC及其95%CI
icc_result <- icc(fit)  # 计算ICC和CI
# 3. 手动计算 ICC 的 95% CI（基于 Delta 方法 + 逻辑回归残差方差 ≈ π²/3）
# ------------------------------------------
# 提取随机效应方差
vc <- as.data.frame(VarCorr(fit))
var_random <- vc$vcov[1]  # 随机截距方差
se_random <- vc$sdcor[1]   # 随机截距的标准误

# 逻辑回归的残差方差（固定值）
var_residual <- pi^2 / 3   # ≈ 3.289868

# 计算 ICC（调整后的）
icc_value <- var_random / (var_random + var_residual)

# 计算 ICC 的标准误（Delta 方法）
se_icc <- se_random * var_residual / (var_random + var_residual)^2

# 计算 95% 置信区间（并限制在 [0, 1] 范围内）
icc_ci_lower <- max(0, icc_value - 1.96 * se_icc)
icc_ci_upper <- min(1, icc_value + 1.96 * se_icc)
# ------------------------------------------

# 4. 输出结果
cat("===================== ICC (Adjusted) =====================\n")
cat(sprintf(
  "ICC = %.3f (95%% CI: [%.3f, %.3f])\n",
  icc_value,
  icc_ci_lower,
  icc_ci_upper
))

# 对比 performance::icc() 的点估计值（验证一致性）
cat("\n对比 performance::icc() 的结果：\n")
print(icc_result)
# --------------- 打印结果 ---------------
cat("===================== Fixed Effects =====================\n")
print(fe_df)

cat("\n=================== Random Effects ===================\n")
print(re)

cat("\n=================== ICC with 95% CI ===================\n")
print(icc_result)

#LBR
# 构建数据框
data <- data.frame(
  center = rep(1:6, each = 2),
  treatment = rep(c("mLP", "GnRH-ant"), 6),
  responder = c(6,6,10,14,2,0,1,0,18,16,2,0),
  total = c(16,18,54,40,11,11,12,9,54,70,6,2)
)

# 拟合GLMM模型
fit <- glmer(cbind(responder, total - responder) ~ treatment + (1 | center),
             family = binomial, data = data)

# --------------- 输出结果 ---------------
# 1. 固定效应（含OR和95%CI）
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

# 2. 随机效应方差
re <- as.data.frame(VarCorr(fit))

# 3. 计算ICC及其95%CI
icc_result <- icc(fit)  # 计算ICC和CI
# 3. 手动计算 ICC 的 95% CI（基于 Delta 方法 + 逻辑回归残差方差 ≈ π²/3）
# ------------------------------------------
# 提取随机效应方差
vc <- as.data.frame(VarCorr(fit))
var_random <- vc$vcov[1]  # 随机截距方差
se_random <- vc$sdcor[1]   # 随机截距的标准误

# 逻辑回归的残差方差（固定值）
var_residual <- pi^2 / 3   # ≈ 3.289868

# 计算 ICC（调整后的）
icc_value <- var_random / (var_random + var_residual)

# 计算 ICC 的标准误（Delta 方法）
se_icc <- se_random * var_residual / (var_random + var_residual)^2

# 计算 95% 置信区间（并限制在 [0, 1] 范围内）
icc_ci_lower <- max(0, icc_value - 1.96 * se_icc)
icc_ci_upper <- min(1, icc_value + 1.96 * se_icc)
# ------------------------------------------

# 4. 输出结果
cat("===================== ICC (Adjusted) =====================\n")
cat(sprintf(
  "ICC = %.3f (95%% CI: [%.3f, %.3f])\n",
  icc_value,
  icc_ci_lower,
  icc_ci_upper
))

# 对比 performance::icc() 的点估计值（验证一致性）
cat("\n对比 performance::icc() 的结果：\n")
print(icc_result)
# --------------- 打印结果 ---------------
cat("===================== Fixed Effects =====================\n")
print(fe_df)

cat("\n=================== Random Effects ===================\n")
print(re)

cat("\n=================== ICC with 95% CI ===================\n")
print(icc_result)


##################ITT
#CPR
data <- data.frame(
  center = rep(1:6, each = 2),
  treatment = rep(c("mLP", "GnRH-ant"), 6),
  responder = c(6,11,13,17,2,1,4,1,22,24,4,0),
  total = c(16,18,60,48,11,12,12,9,54,70,6,2)
)

fit <- glmer(cbind(responder, total - responder) ~ treatment + (1 | center),
             family = binomial, data = data)

# --------------- 输出结果 ---------------
# 1. 固定效应（含OR和95%CI）
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

# 2. 随机效应方差
re <- as.data.frame(VarCorr(fit))

# 3. 计算ICC及其95%CI
icc_result <- icc(fit)  # 计算ICC和CI
# 3. 手动计算 ICC 的 95% CI（基于 Delta 方法 + 逻辑回归残差方差 ≈ π²/3）
# ------------------------------------------
# 提取随机效应方差
vc <- as.data.frame(VarCorr(fit))
var_random <- vc$vcov[1]  # 随机截距方差
se_random <- vc$sdcor[1]   # 随机截距的标准误
# 逻辑回归的残差方差（固定值）
var_residual <- pi^2 / 3   # ≈ 3.289868
# 计算 ICC（调整后的）
icc_value <- var_random / (var_random + var_residual)
# 计算 ICC 的标准误（Delta 方法）
se_icc <- se_random * var_residual / (var_random + var_residual)^2
# 计算 95% 置信区间（并限制在 [0, 1] 范围内）
icc_ci_lower <- max(0, icc_value - 1.96 * se_icc)
icc_ci_upper <- min(1, icc_value + 1.96 * se_icc)
# ------------------------------------------
# 4. 输出结果
cat("===================== ICC (Adjusted) =====================\n")
cat(sprintf(
  "ICC = %.3f (95%% CI: [%.3f, %.3f])\n",
  icc_value,
  icc_ci_lower,
  icc_ci_upper
))

# 对比 performance::icc() 的点估计值（验证一致性）
cat("\n对比 performance::icc() 的结果：\n")
print(icc_result)
# --------------- 打印结果 ---------------
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

# --------------- 输出结果 ---------------
# 1. 固定效应（含OR和95%CI）
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
# 2. 随机效应方差
re <- as.data.frame(VarCorr(fit))
# 3. 计算ICC及其95%CI
icc_result <- icc(fit)  # 计算ICC和CI
# 3. 手动计算 ICC 的 95% CI（基于 Delta 方法 + 逻辑回归残差方差 ≈ π²/3）
# ------------------------------------------
# 提取随机效应方差
vc <- as.data.frame(VarCorr(fit))
var_random <- vc$vcov[1]  # 随机截距方差
se_random <- vc$sdcor[1]   # 随机截距的标准误
# 逻辑回归的残差方差（固定值）
var_residual <- pi^2 / 3   # ≈ 3.289868
# 计算 ICC（调整后的）
icc_value <- var_random / (var_random + var_residual)
# 计算 ICC 的标准误（Delta 方法）
se_icc <- se_random * var_residual / (var_random + var_residual)^2
# 计算 95% 置信区间（并限制在 [0, 1] 范围内）
icc_ci_lower <- max(0, icc_value - 1.96 * se_icc)
icc_ci_upper <- min(1, icc_value + 1.96 * se_icc)
# ------------------------------------------
# 4. 输出结果
cat("===================== ICC (Adjusted) =====================\n")
cat(sprintf(
  "ICC = %.3f (95%% CI: [%.3f, %.3f])\n",
  icc_value,
  icc_ci_lower,
  icc_ci_upper
))

# 对比 performance::icc() 的点估计值（验证一致性）
cat("\n对比 performance::icc() 的结果：\n")
print(icc_result)
# --------------- 打印结果 ---------------
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
  total = c(16,18,60,48,11,12,12,9,54,70,6,2)
)
fit <- glmer(cbind(responder, total - responder) ~ treatment + (1 | center),
             family = binomial, data = data)
# --------------- 输出结果 ---------------
# 1. 固定效应（含OR和95%CI）
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
# 2. 随机效应方差
re <- as.data.frame(VarCorr(fit))
# 3. 计算ICC及其95%CI
icc_result <- icc(fit)  # 计算ICC和CI
# 3. 手动计算 ICC 的 95% CI（基于 Delta 方法 + 逻辑回归残差方差 ≈ π²/3）
# ------------------------------------------
# 提取随机效应方差
vc <- as.data.frame(VarCorr(fit))
var_random <- vc$vcov[1]  # 随机截距方差
se_random <- vc$sdcor[1]   # 随机截距的标准误
# 逻辑回归的残差方差（固定值）
var_residual <- pi^2 / 3   # ≈ 3.289868
# 计算 ICC（调整后的）
icc_value <- var_random / (var_random + var_residual)
# 计算 ICC 的标准误（Delta 方法）
se_icc <- se_random * var_residual / (var_random + var_residual)^2
# 计算 95% 置信区间（并限制在 [0, 1] 范围内）
icc_ci_lower <- max(0, icc_value - 1.96 * se_icc)
icc_ci_upper <- min(1, icc_value + 1.96 * se_icc)
# ------------------------------------------
# 4. 输出结果
cat("===================== ICC (Adjusted) =====================\n")
cat(sprintf(
  "ICC = %.3f (95%% CI: [%.3f, %.3f])\n",
  icc_value,
  icc_ci_lower,
  icc_ci_upper
))

# 对比 performance::icc() 的点估计值（验证一致性）
cat("\n对比 performance::icc() 的结果：\n")
print(icc_result)
# --------------- 打印结果 ---------------
cat("===================== Fixed Effects =====================\n")
print(fe_df)

cat("\n=================== Random Effects ===================\n")
print(re)

cat("\n=================== ICC with 95% CI ===================\n")
print(icc_result)

###卡方
library(stats)
a<-matrix(c(51,61,54,77),nrow=2)
chisq.test(a)


b<-matrix(c(6,19,2,14),nrow=2)
fisher.test(b)
chisq.test(b)


c<-matrix(c(9,66,18,67),nrow=2)
chisq.test(c)
fisher.test(c)

d<-matrix(c(23,23,18,35),nrow=2)
chisq.test(d)

e<-matrix(c(70,198,121,273),nrow=2)
fisher.test(e)
chisq.test(e)
