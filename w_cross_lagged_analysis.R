library(lavaan) 
rm(list=ls()) 

data <- read.csv("/Users/loadWXF/IMA_analy/cross_lagged_ana/audit_dep8_FU1FU2.csv") 

regression_model_X1 <- lm(dep8_FU1 ~ Gender_Male + Site1 + Site2 + Site3 + Site4 + Site5 + Site6 + Site7 + Handedness, data = data)
x1 <- regression_model_X1$residuals

regression_model_X2 <- lm(dep8_FU2 ~ Gender_Male + Site1 + Site2 + Site3 + Site4 + Site5 + Site6 + Site7 + Handedness, data = data)
x4 <- regression_model_X2$residuals

regression_model_Y1 <- lm(audit_freq_FU1 ~ Gender_Male +Site1 + Site2 + Site3 + Site4 + Site5 + Site6 + Site7 + Handedness, data = data)
y1 <- regression_model_Y1$residuals

regression_model_Y2 <- lm(audit_freq_FU2 ~  Gender_Male +Site1 + Site2 + Site3 + Site4 + Site5 + Site6 + Site7 + Handedness, data = data)
y4 <- regression_model_Y2$residuals

data_cross <- cbind(x1,x4,y1,y4)

clpmModel <- #yes, "Model" is redundant
  '
kappa =~ 1*x1 + 1*x4
omega =~ 1*y1 + 1*y4

x1 ~ mu1*1 #intercepts 
x4 ~ mu4*1

y1 ~ pi1*1
y4 ~ pi4*1

kappa ~~ 0*kappa #variance nope 
omega ~~ 0*omega #variance nope
kappa ~~ 0*omega #covariance not even 

#laten vars for AR and cross-lagged effects 
p1 =~ 1*x1 #each factor loading set to 1   
p4 =~ 1*x4
q1 =~ 1*y1
q4 =~ 1*y4

p4 ~ alpha2*p1 + beta2*q1

q4 ~ delta2*q1 + gamma2*p1

p1 ~~ p1 #variance 
p4 ~~ u4*p4
q1 ~~ q1 #variance 
q4 ~~ v4*q4

p1 ~~ q1 #p1 and q1 covariance  
p4 ~~ q4 #p4 and q4 covariance'

fitCLPM <- lavaan(clpmModel, data = data_cross,
                  missing = 'ML', #for the missing data!
                  int.ov.free = F,
                  int.lv.free = F,
                  auto.fix.first = F,
                  auto.fix.single = F,
                  auto.cov.lv.x = F,
                  auto.cov.y = F,
                  auto.var = F)

summary(fitCLPM, standardized = T, fit.measures = TRUE)


model_summary<-summary(fitCLPM, standardized = T, fit.measures = TRUE)
View(model_summary[["pe"]])



##########################################  Gender CLPM  ##############################################
rm(list=ls())  

data <- read.csv("/Users/loadWXF/IMA_analy/cross_lagged_ana/audit_dep8_FU1FU2.csv") # 读取数据文件

sex_data <- data[data$Gender_Male == 1, ]


regression_model_X1 <- lm(dep8_FU1 ~ Site1 + Site2 + Site3 + Site4 + Site5 + Site6 + Site7 + Handedness, data = sex_data)
x1 <- regression_model_X1$residuals

regression_model_X2 <- lm(dep8_FU2 ~ Site1 + Site2 + Site3 + Site4 + Site5 + Site6 + Site7 + Handedness, data = sex_data)
x4 <- regression_model_X2$residuals

regression_model_Y1 <- lm(audit_freq_FU1 ~ Site1 + Site2 + Site3 + Site4 + Site5 + Site6 + Site7 + Handedness, data = sex_data)
y1 <- regression_model_Y1$residuals

regression_model_Y2 <- lm(audit_freq_FU2 ~  Site1 + Site2 + Site3 + Site4 + Site5 + Site6 + Site7 + Handedness, data = sex_data)
y4 <- regression_model_Y2$residuals

data_cross <- cbind(x1,x4,y1,y4)

clpmModel <- #yes, "Model" is redundant
  '
kappa =~ 1*x1 + 1*x4
omega =~ 1*y1 + 1*y4

x1 ~ mu1*1 #intercepts 
x4 ~ mu4*1

y1 ~ pi1*1
y4 ~ pi4*1

kappa ~~ 0*kappa #variance nope 
omega ~~ 0*omega #variance nope
kappa ~~ 0*omega #covariance not even 

#laten vars for AR and cross-lagged effects 
p1 =~ 1*x1 #each factor loading set to 1   
p4 =~ 1*x4
q1 =~ 1*y1
q4 =~ 1*y4

p4 ~ alpha2*p1 + beta2*q1

q4 ~ delta2*q1 + gamma2*p1

p1 ~~ p1 #variance 
p4 ~~ u4*p4
q1 ~~ q1 #variance 
q4 ~~ v4*q4

p1 ~~ q1 #p1 and q1 covariance  
p4 ~~ q4 #p4 and q4 covariance'

fitCLPM <- lavaan(clpmModel, data = data_cross,
                  missing = 'ML', 
                  int.ov.free = F,
                  int.lv.free = F,
                  auto.fix.first = F,
                  auto.fix.single = F,
                  auto.cov.lv.x = F,
                  auto.cov.y = F,
                  auto.var = F)

summary(fitCLPM, standardized = T, fit.measures = TRUE)


model_summary<-summary(fitCLPM, standardized = T, fit.measures = TRUE)
View(model_summary[["pe"]])
