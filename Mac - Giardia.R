#Generate, plot, and test regression models for Macarthur Giardia

#Install and load necessary packages
install.packages("pscl")
install.packages("forecast")
install.packages("glarma")
install.packages("MASS")
install.packages("statmod")
install.packages("DHARMa")
install.packages("sanwich")
install.packages("mfx")
library(pscl)
library(forecast)
library(glarma)
library(MASS)
library(statmod)
library(DHARMa)
library(sandwich)
library(mfx)

#Generate the recovery rate variable
Macqmra$rateG<-Macqmra$S1_NG/Macqmra$S1_offsetG

#============ Mac - Giardia - neg bin - phase model =================

#Fit neg bin
mac_g_nb_step=glm.nb(Macqmra$S1_NG ~ Macqmra$Week+
                       Macqmra$Phase+
                       offset(log(Macqmra$S1_offsetG)))
  
#Predicted values of the recovery rate
mac_g_nb_step_rate_pred = mac_g_nb_step$fitted.values/Macqmra$S1_offsetG
plot(Macqmra$Week, mac_g_nb_step_rate_pred, type="l", col="blue", lwd=2)

#Test for overall uniformity of the fitted neg bin model 
#(KS test) which can detect heteroskedacticity
testUniformity(simulateResiduals(mac_g_nb_step)) 

# ##Robust standard error - beta and IRR. 
# #Not indicated due to lack of heteroskedacity.
# #Robust standard errors for coefficients (beta) (and 95% CI).
# #Applies a correction factor to achieve heteroskedacticity-consistent 
# #estimate of robust standard error.
# coeftest(mac_g_nb_step, vcov = vcovHC(mac_g_nb_step, type="HC0"))
# coefci(mac_g_nb_step,vcov = vcovHC(mac_g_nb_step, type="HC0"))
# #IRR, IRR SE, z, and p values, 
# #this method is necessary to find these using robust standard error.
# mac_g_nb_step_negbinirr = negbinirr(formula = mac_g_nb_step$call, 
#                                     data = Macqmra, 
#                                     robust = TRUE)
# mac_g_nb_step_negbinirr

#Non-robust standard error for coefficients (beta) (and 95% CI)
#Indicated due to homoskedacticity
coeftest(mac_g_nb_step, vcov = vcov(mac_g_nb_step))
coefci(mac_g_nb_step,vcov = vcov(mac_g_nb_step))
#non-robust IRR, IRR SE, z, and p values:
mac_g_nb_step_negbinirr = negbinirr(formula = mac_g_nb_step$call, 
                                    data = Macqmra, 
                                    robust = FALSE)
mac_g_nb_step_negbinirr

#Plot ACF and PACF to check for AR and MA processes
Acf(qresiduals(mac_g_nb_step))
Pacf(qresiduals(mac_g_nb_step))
hist(qresiduals(mac_g_nb_step))
plot(simulateResiduals(mac_g_nb_step)) 

#Perform tests
odTest(mac_g_nb_step) #Test overdispersion (H0 = Poisson; H1 = overdispersed)
testDispersion(simulateResiduals(mac_g_nb_step)) #Test overdispersion of n.b.
testZeroInflation(simulateResiduals(mac_g_nb_step)) #Test zero-inflation
testTemporalAutocorrelation(simulateResiduals(mac_g_nb_step), 
                            time=Macqmra$Week) #D-W test autocorr. at lag = 1

#Extract model summary information
summary(mac_g_nb_step)
mac_g_nb_step$theta
mac_g_nb_step$SE.theta
mac_g_nb_step$aic

#GLARMA modelling not required due to lack of AR or MA signatures.

#Test additional covariates

# Colour
#Recode "<2" as 2, and convert to numeric
Macqmra$Col400cont = Macqmra$Col400
Macqmra$Col400cont[Macqmra$Col400cont == "<2"] = 2
Macqmra$Col400cont = as.numeric(Macqmra$Col400cont)
mac_g_nb_step_colour=glm.nb(Macqmra$S1_NG ~ Macqmra$Week+
                              Macqmra$Phase+
                              Macqmra$Col400cont+
                              offset(log(Macqmra$S1_offsetG)))
coeftest(mac_g_nb_step_colour, vcov = vcov(mac_g_nb_step_colour))
coefci(mac_g_nb_step_colour,vcov = vcov(mac_g_nb_step_colour))
mac_g_nb_step_negbinirr = negbinirr(formula = mac_g_nb_step_colour$call, 
                                    data = Macqmra, 
                                    robust = FALSE)
mac_g_nb_step_negbinirr
mac_g_nb_step_colour$aic
mac_g_nb_step_colour$theta

#Temperature
mac_g_nb_step_temp=glm.nb(Macqmra$S1_NG ~ Macqmra$Week+
                            Macqmra$Phase+
                            Macqmra$Temp+
                            offset(log(Macqmra$S1_offsetG)))
coeftest(mac_g_nb_step_temp, vcov = vcov(mac_g_nb_step_temp))
coefci(mac_g_nb_step_temp,vcov = vcov(mac_g_nb_step_temp))
mac_g_nb_step_negbinirr = negbinirr(formula = mac_g_nb_step_temp$call, 
                                    data = Macqmra, 
                                    robust = FALSE)
mac_g_nb_step_negbinirr
mac_g_nb_step_temp$aic
mac_g_nb_step_temp$theta

#Turbidity
mac_g_nb_step_turbidity=glm.nb(Macqmra$S1_NG ~ Macqmra$Week+
                                 Macqmra$Phase+
                                 Macqmra$Turbidity+
                                 offset(log(Macqmra$S1_offsetG)))
coeftest(mac_g_nb_step_turbidity, vcov = vcov(mac_g_nb_step_turbidity))
coefci(mac_g_nb_step_turbidity,vcov = vcov(mac_g_nb_step_turbidity))
mac_g_nb_step_negbinirr = negbinirr(formula = mac_g_nb_step_turbidity$call, 
                                    data = Macqmra, 
                                    robust = FALSE)
mac_g_nb_step_negbinirr
mac_g_nb_step_turbidity$aic
mac_g_nb_step_turbidity$theta

#Colour and temperature
mac_g_nb_step_coltemp=glm.nb(Macqmra$S1_NG ~ Macqmra$Week+
                               Macqmra$Phase+
                               Macqmra$Col400cont+
                               Macqmra$Temp+
                               offset(log(Macqmra$S1_offsetG)))
coeftest(mac_g_nb_step_coltemp, vcov = vcov(mac_g_nb_step_coltemp))
coefci(mac_g_nb_step_coltemp,vcov = vcov(mac_g_nb_step_coltemp))
mac_g_nb_step_negbinirr = negbinirr(formula = mac_g_nb_step_coltemp$call, 
                                    data = Macqmra, 
                                    robust = FALSE)
mac_g_nb_step_negbinirr
mac_g_nb_step_coltemp$aic
mac_g_nb_step_coltemp$theta


#============ Mac - Giardia - neg bin - phase and interact model ============

#Fit neg bin
mac_g_nb_steprate=glm.nb(Macqmra$S1_NG ~ Macqmra$Week*Macqmra$Phase+
                           offset(log(Macqmra$S1_offsetG)))

#Predicted values of the recovery rate
mac_g_nb_steprate_rate_pred = 
  mac_g_nb_steprate$fitted.values/Macqmra$S1_offsetG
plot(Macqmra$Week, mac_g_nb_steprate_rate_pred, type="l", col="blue", lwd=2)

#Test for overall uniformity of the fitted neg bin model 
#(KS test) which can detect heteroskedacticity
testUniformity(simulateResiduals(mac_g_nb_steprate)) 

# #Robust standard errors for coefficients (beta) (and 95% CI)
# #Not indicated due to lack of heteroskedacity
# coeftest(mac_g_nb_steprate, vcov = vcovHC(mac_g_nb_steprate, type="HC0"))
# coefci(mac_g_nb_steprate,vcov = vcovHC(mac_g_nb_steprate, type="HC0"))
# #IRR, IRR SE, z, and p values
# #this method is necessary to find these using robust standard error
# mac_g_nb_steprate_negbinirr = negbinirr(formula = mac_g_nb_steprate$call, 
#                                         data = Macqmra, 
#                                         robust = TRUE)
# mac_g_nb_steprate_negbinirr

#Non-robust standard error for coefficients (beta) (and 95% CI)
#Indicated due to homoskedacticity
coeftest(mac_g_nb_steprate, vcov = vcov(mac_g_nb_steprate))
coefci(mac_g_nb_steprate,vcov = vcov(mac_g_nb_steprate))
#Non-robust IRR, IRR SE, z, and p values:
mac_g_nb_steprate_negbinirr = negbinirr(formula = mac_g_nb_steprate$call, 
                                        data = Macqmra, 
                                        robust = FALSE)
mac_g_nb_steprate_negbinirr

#Plot ACF and PACF to check for AR and MA processes
Acf(qresiduals(mac_g_nb_steprate))
Pacf(qresiduals(mac_g_nb_steprate))
hist(qresiduals(mac_g_nb_steprate))
plot(simulateResiduals(mac_g_nb_steprate)) 

#Perform tests
odTest(mac_g_nb_steprate) #Test overdispersion (H0 = Poiss.; H1 = overdisp.)
testDispersion(simulateResiduals(mac_g_nb_steprate)) #Test overdisp. of n.b.
testZeroInflation(simulateResiduals(mac_g_nb_steprate)) #Test zero-inflation
testTemporalAutocorrelation(simulateResiduals(mac_g_nb_steprate), 
                            time=Macqmra$Week) #D-W test autocorr. at lag = 1

summary(mac_g_nb_steprate)
mac_g_nb_steprate$theta
mac_g_nb_steprate$SE.theta
mac_g_nb_steprate$aic


#============ Mac - Giardia - ZINB - phase model =================

#Fit zero-inflated neg bin
mac_g_zinb_step = zeroinfl(Macqmra$S1_NG ~ Macqmra$Week+Macqmra$Phase, 
                           offset=offset(log(Macqmra$S1_offsetG)), 
                           dist="negbin")

#Plot results
plot(mac_g_zinb_step$fitted.values/Macqmra$S1_offsetG,type="l")

#Extract model summary information
summary(mac_g_zinb_step)
extractAIC(mac_g_zinb_step)

#plot ACF and PACF to check for AR and MA processes
par(mfrow=c(1,2)) 
Acf(mac_g_zinb_step$residuals, lag.max=NULL, type=c("correlation"), 
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(mac_g_zinb_step$residuals, lag.max=NULL, plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)
hist(qresiduals(mac_g_zinb_step))


#============ Mac - Giardia - ZINB - phase and interact model ============

#Fit zero-inflated neg bin
mac_g_zinb_steprate = zeroinfl(Macqmra$S1_NG ~ Macqmra$Week*Macqmra$Phase, 
                               offset=offset(log(Macqmra$S1_offsetG)), 
                               dist="negbin")

#Plot results
plot(mac_g_zinb_steprate$fitted.values/Macqmra$S1_offsetG,type="l")

#Extract model summary information
summary(mac_g_zinb_steprate)
extractAIC(mac_g_zinb_steprate)

#Plot ACF and PACF to check for AR and MA processes
Acf(mac_g_zinb_steprate$residuals, lag.max=NULL, type=c("correlation"),
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(mac_g_zinb_steprate$residuals, lag.max=NULL, plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)
hist(qresiduals(mac_g_zinb_steprate))


#============ Mac - Giardia - Poisson - phase model ===========

#Fit Poisson
mac_g_p_step=glm(Macqmra$S1_NG ~ Macqmra$Week+Macqmra$Phase, 
                 offset=offset(log(Macqmra$S1_offsetG)),
                 family=poisson)

#Plot results
plot(mac_g_p_step$fitted.values/Macqmra$S1_offsetG, type="l")

#Extract model summary information
summary(mac_g_p_step)

#Plot ACF and PACF to check for AR and MA processes
Acf(qresiduals(mac_g_p_step))
Pacf(qresiduals(mac_g_p_step))
hist(qresiduals(mac_g_p_step))
plot(simulateResiduals(mac_g_p_step)) 

#Perform tests
testOverdispersion(simulateResiduals(mac_g_p_step)) #Test overdisp. of Poiss.
testZeroInflation(simulateResiduals(mac_g_p_step)) #Test zero-inflation
testUniformity(simulateResiduals(mac_g_p_step)) #Detect heteroskedacticity
testTemporalAutocorrelation(simulateResiduals(mac_g_p_step), 
                            time=Macqmra$Week) #D-W test autocorr. at lag = 1

#Overdispersion was detected. The Poisson approach would be abandoned.


#============ Mac - Giardia - Poisson - phase and interact mdoel =========

#Fit Poisson
mac_g_p_steprate=glm(Macqmra$S1_NG ~ Macqmra$Week*Macqmra$Phase, 
                     offset=offset(log(Macqmra$S1_offsetG)),
                     family=poisson)

#Plot results
plot(mac_g_p_steprate$fitted.values/Macqmra$S1_offsetG, type="l")

#Extract model summary information
summary(mac_g_p_steprate)

#Plot ACF and PACF to check for AR and MA processes
Acf(qresiduals(mac_g_p_steprate))
Pacf(qresiduals(mac_g_p_steprate))
hist(qresiduals(mac_g_p_steprate))
plot(simulateResiduals(mac_g_p_step)) 

#Perform tests
testDispersion(simulateResiduals(mac_g_p_steprate)) #Test overdisp. of Poiss.
testZeroInflation(simulateResiduals(mac_g_p_steprate)) #Test zero-inflation
testTemporalAutocorrelation(simulateResiduals(mac_g_p_steprate), 
                            time=Macqmra$Week) #D-W test autocorr. at lag = 1

#Overdispersion was detected. The Poisson approach would be abandoned.


#============ Mac - Giardia - ZIP - phase model ================

#Fit ZIP
mac_g_zip_step = zeroinfl(Macqmra$S1_NG ~ Macqmra$Week+Macqmra$Phase, 
                          offset=offset(log(Macqmra$S1_offsetG)), 
                          dist="poisson")

#Fit a simple linear regression
mac_g_zip_step_glm = glm(formula = mac_g_zip_step$fitted.values ~ 
                           Macqmra$Week+Macqmra$Phase)

#Plot results
plot(mac_g_zip_step$fitted.values/Macqmra$S1_offsetG,type="l")

#Extract model summary information
summary(mac_g_zip_step)
extractAIC(mac_g_zip_step)

#Plot ACF and PACF to check for AR and MA processes
Acf(mac_g_zip_step$residuals, lag.max=NULL, type=c("correlation"), 
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(mac_g_zip_step$residuals, lag.max=NULL, plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)
hist(qresiduals(mac_g_zip_step))


#============ Mac - Giardia - ZIP - phase and interact model ===========

#Fit ZIP
mac_g_zip_steprate = zeroinfl(Macqmra$S1_NG ~ Macqmra$Week*Macqmra$Phase, 
                              offset=offset(log(Macqmra$S1_offsetG)), 
                              dist="poisson")

#Plot results
plot(mac_g_zip_steprate$fitted.values/Macqmra$S1_offsetG,type="l")

#Extract model summary information
summary(mac_g_zip_steprate)
extractAIC(mac_g_zip_steprate)

#Plot ACF and PACF to check for AR and MA processes
Acf(mac_g_zip_steprate$residuals, lag.max=NULL,type=c("correlation"), 
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(mac_g_zip_steprate$residuals, lag.max=NULL, plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)
hist(qresiduals(mac_g_zip_steprate))
