rm(list=ls())  

setwd("/Users/loadWXF/IMA_analy/Apaper_result/Mediation/medi_4items")
getwd()

data<- read.csv("medi_FU1.csv") 

regression_model_X1 <- lm(h_mean ~ Gender_Male + Site1 + Site2 + Site3 + Site4 + Site5 + Site6 + Site7 + Handedness, data = data)
M1 <- regression_model_X1$residuals

regression_model_X2 <- lm(imp_mean ~ Gender_Male + Site1 + Site2 + Site3 + Site4 + Site5 + Site6 + Site7 + Handedness, data = data)
M2 <- regression_model_X2$residuals

regression_model_Y1 <- lm(SCORESUM ~ Gender_Male +Site1 + Site2 + Site3 + Site4 + Site5 + Site6 + Site7 + Handedness, data = data)
X <- regression_model_Y1$residuals

regression_model_Y2 <- lm(audit_sum ~  Gender_Male +Site1 + Site2 + Site3 + Site4 + Site5 + Site6 + Site7 + Handedness, data = data)
Y <- regression_model_Y2$residuals

data_cross <- cbind(M1,M2,X,Y)

library(lavaan)


model0=
  "
Y ~ z*X
"
fit0<-sem(model0,data = data_cross, se="boot", bootstrap = 10000)
fitsum0 <-lavaan::parameterEstimates(fit0, boot.ci.type = "bca.simple")
oriexp_dep<-summary(fit0,standardized=T)


model1 <- "
  M1 ~ a1 * X
  M2 ~ a2 * X + d21 * M1
  Y ~  cp * X + b1  * M1 + b2 * M2
  ind_eff := a1 * d21 * b2
  de_eff := cp
"

fit1 <- lavaan::sem(model = model1, data = data_cross, se = "boot", bootstrap = 10000)
fitsum1 <-lavaan::parameterEstimates(fit1, boot.ci.type = "bca.simple")
mediationpath1_dep<-summary(fit1,standardized=T)


indirect_effect <- mediationpath1_dep[["pe"]][["est"]][11]#ind_eff
total_effect <- mediationpath1_dep[["pe"]][["est"]][11] + mediationpath1_dep[["pe"]][["est"]][12]#de_eff
mediation_proportion <- indirect_effect / total_effect



model2=
  "
   M1 ~ a1*X
   Y ~ b1*M1 + b2*X

ie := a1*b1
de := b2
     
"

fit2<-sem(model2,data = data_cross, se="boot", bootstrap = 10000)
fit_summary2<-summary(fit2,standardized=T)

indirect_effect <- fit_summary2[["pe"]][["est"]][7]#ind_eff
total_effect <- fit_summary2[["pe"]][["est"]][7] + fit_summary2[["pe"]][["est"]][8]#de_eff
mediation_proportion <- indirect_effect / total_effect



model3=
  "
   M2 ~ a2*X
   Y ~ b2*M2 + b3*X

ie := a2*b2
de := b3
"

fit3<-sem(model3,data = data_cross, se="boot", bootstrap = 10000)
fit_summary3<-summary(fit3,standardized=T)


indirect_effect <- fit_summary3[["pe"]][["est"]][7]#ind_eff
total_effect <- fit_summary3[["pe"]][["est"]][7] + fit_summary3[["pe"]][["est"]][8]#de_eff
mediation_proportion <- indirect_effect / total_effect



model4=
  "
   M2 ~ a4*M1
   Y ~ b4*M2 + b5*M1

ie := a4*b4
de := b5
"

fit4<-sem(model4,data = data_cross, se="boot", bootstrap = 10000)
fit_summary4<-summary(fit4,standardized=T)
fitsum4 <-lavaan::parameterEstimates(fit4, boot.ci.type = "bca.simple")


indirect_effect <- fit_summary4[["pe"]][["est"]][7]#ind_eff
total_effect <- fit_summary4[["pe"]][["est"]][7] + fit_summary4[["pe"]][["est"]][8]#de_eff
mediation_proportion <- indirect_effect / total_effect



model5=
  "
   M1 ~ a5*X
   M2 ~ b5*M1 + b6*X

ie := a5*b5
de := b6
"

fit5<-sem(model5,data = data_cross, se="boot", bootstrap = 10000)
fit_summary5<-summary(fit5,standardized=T)


indirect_effect <- fit_summary5[["pe"]][["est"]][7]#ind_eff
total_effect <- fit_summary5[["pe"]][["est"]][7] + fit_summary5[["pe"]][["est"]][8]#de_eff
mediation_proportion <- indirect_effect / total_effect
