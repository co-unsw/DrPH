#Generate, plot, and test regression models for Macarthur Cryptosporidium

#Install and load necessary packages
install.packages("pscl")
install.packages("forecast")
install.packages("MASS")
install.packages("statmod")
install.packages("glarma")
install.packages("DHARMa")
install.packages("mfx")
library(pscl)
library(forecast)
library(MASS)
library(statmod)
library(glarma)
library(DHARMa)
library(mfx)

#Generate the recovery rate variable
Macqmra$rate<-Macqmra$S1_N/Macqmra$S1_offset

#============ Mac - Crypto - neg bin - phase model =================
#Fit neg bin
mac_c_nb_step=glm.nb(Macqmra$S1_N ~ Macqmra$Week+
                       Macqmra$Phase+
                       offset(log(Macqmra$S1_offset)))
plot(Macqmra$Week, mac_c_nb_step_rate_pred, type="l", col="blue", lwd=2)

#Test for overall uniformity of the fitted neg bin model 
#(KS test) which can detect heteroskedacticity
testUniformity(simulateResiduals(mac_c_nb_step))

# ##Robust standard error - beta and IRR. 
# #Not indicated due to lack of heteroskedacity.
# #Robust standard errors for coefficients (beta) (and 95% CI).
# #Applies a correction factor to achieve heteroskedacticity-consistent 
# #estimate of robust standard error.
# coeftest(mac_c_nb_step, vcov = vcovHC(mac_c_nb_step, type="HC0"))
# coefci(mac_c_nb_step,vcov = vcovHC(mac_c_nb_step, type="HC0"))
# #IRR, IRR SE, z, and p values, 
# #this method is necessary to find these using robust standard error.
# mac_c_nb_steprate_negbinirr = negbinirr(formula = mac_c_nb_step$call,
#                                         data = Macqmra, 
#                                         robust = TRUE)
# mac_c_nb_steprate_negbinirr

#Non-robust standard error for coefficients (beta) (and 95% CI)
#Indicated due to homoskedacticity
coeftest(mac_c_nb_step, vcov = vcov(mac_c_nb_step))
coefci(mac_c_nb_step,vcov = vcov(mac_c_nb_step))
#Non-robust IRR, IRR SE, z, and p values:
mac_c_nb_step_negbinirr = negbinirr(formula = mac_c_nb_step$call, 
                                    data = Macqmra, 
                                    robust = FALSE)
mac_c_nb_step_negbinirr

#Plot ACF and PACF to check for AR and MA processes
Acf(qresiduals(mac_c_nb_step))
Pacf(qresiduals(mac_c_nb_step))
hist(qresiduals(mac_c_nb_step))
plot(simulateResiduals(mac_c_nb_step)) 

#Perform tests
odTest(mac_c_nb_step) #Test overdispersion (H0 = Poisson; H1 = overdispersed)
testDispersion(simulateResiduals(mac_c_nb_step)) #Test overdispersion of n.b.
testZeroInflation(simulateResiduals(mac_c_nb_step)) #Test zero-inflation
testTemporalAutocorrelation(simulateResiduals(mac_c_nb_step), 
                            time=Macqmra$Week) #D-W test autocorr. at lag = 1

#Extract model summary information
summary(mac_c_nb_step)
mac_c_nb_step$theta
mac_c_nb_step$SE.theta
mac_c_nb_step$aic

#GLARMA modelling not required due to lack of AR or MA signatures.

#Test additional covariates

#Colour
#Recode "<2" as 2, and convert to numeric
Macqmra$Col400cont = Macqmra$Col400
Macqmra$Col400cont[Macqmra$Col400cont == "<2"] = 2
Macqmra$Col400cont = as.numeric(Macqmra$Col400cont)
mac_c_nb_step_colour=glm.nb(Macqmra$S1_N ~ Macqmra$Week+
                              Macqmra$Phase+
                              Macqmra$Col400cont+
                              offset(log(Macqmra$S1_offset)))
coeftest(mac_c_nb_step_colour, vcov = vcov(mac_c_nb_step_colour))
coefci(mac_c_nb_step_colour,vcov = vcov(mac_c_nb_step_colour))
mac_c_nb_step_negbinirr = negbinirr(formula = mac_c_nb_step_colour$call, 
                                    data = Macqmra, 
                                    robust = FALSE)
mac_c_nb_step_negbinirr
mac_c_nb_step_colour$aic
mac_c_nb_step_colour$theta
mac_c_nb_step_colour$SE.theta

#Temperature
mac_c_nb_step_temp=glm.nb(Macqmra$S1_N ~ Macqmra$Week+
                            Macqmra$Phase+
                            Macqmra$Temp+
                            offset(log(Macqmra$S1_offset)))
coeftest(mac_c_nb_step_temp, vcov = vcov(mac_c_nb_step_temp))
coefci(mac_c_nb_step_temp,vcov = vcov(mac_c_nb_step_temp))
mac_c_nb_step_negbinirr = negbinirr(formula = mac_c_nb_step_temp$call, 
                                    data = Macqmra, robust = FALSE)
mac_c_nb_step_negbinirr
mac_c_nb_step_temp$aic
mac_c_nb_step_temp$theta
mac_c_nb_step_temp$SE.theta


#Turbidity
mac_c_nb_step_turbidity=glm.nb(Macqmra$S1_N ~ Macqmra$Week+
                                 Macqmra$Phase+
                                 Macqmra$Turbidity+
                                 offset(log(Macqmra$S1_offset)))
coeftest(mac_c_nb_step_turbidity, vcov = vcov(mac_c_nb_step_turbidity))
coefci(mac_c_nb_step_turbidity,vcov = vcov(mac_c_nb_step_turbidity))
mac_c_nb_step_negbinirr = negbinirr(formula = mac_c_nb_step_turbidity$call, 
                                    data = Macqmra, robust = FALSE)
mac_c_nb_step_negbinirr
mac_c_nb_step_turbidity$aic
mac_c_nb_step_turbidity$theta
mac_c_nb_step_turbidity$SE.theta


#============ Mac - Crypto - neg bin - step and rate ============

#Fit neg bin
mac_c_nb_steprate=glm.nb(Macqmra$S1_N ~ Macqmra$Week*Macqmra$Phase+
                           offset(log(Macqmra$S1_offset)))

#Predicted values of the recovery rate
mac_c_nb_steprate_rate_pred = 
  ac_c_nb_steprate$fitted.values/Macqmra$S1_offset
plot(Macqmra$Week, mac_c_nb_steprate_rate_pred, type="l", col="blue", lwd=2)

# #Robust standard errors for coefficients (beta) (and 95% CI)
# #Not indicated due to lack of heteroskedacity
# coeftest(mac_c_nb_steprate, vcov = vcovHC(mac_c_nb_steprate, type="HC0"))
# coefci(mac_c_nb_steprate,vcov = vcovHC(mac_c_nb_steprate, type="HC0"))
# #IRR, IRR SE, z, and p values
# #this method is necessary to find these using robust standard error
# mac_c_nb_steprate_negbinirr = negbinirr(formula = mac_c_nb_steprate$call, 
#                                         data = Macqmra, 
#                                         robust = TRUE)
# mac_c_nb_steprate_negbinirr

#Non-robust standard error for coefficients (beta) (and 95% CI)
#Indicated due to homoskedacticity
coeftest(mac_c_nb_steprate, vcov = vcov(mac_c_nb_steprate))
coefci(mac_c_nb_steprate,vcov = vcovHC(mac_c_nb_steprate))
#non-robust IRR, IRR SE, z, and p values:
mac_c_nb_steprate_negbinirr = negbinirr(formula = mac_c_nb_steprate$call, 
                                        data = Macqmra, 
                                        robust = FALSE)
mac_c_nb_steprate_negbinirr

#Plot ACF and PACF to check for AR and MA processes
Acf(qresiduals(mac_c_nb_steprate))
Pacf(qresiduals(mac_c_nb_steprate))
hist(qresiduals(mac_c_nb_steprate))
plot(simulateResiduals(mac_c_nb_steprate)) 

#Perform tests
odTest(mac_c_nb_steprate) #Test overdispersion (H0 = Poiss.; H1 = overdisp.)
testDispersion(simulateResiduals(mac_c_nb_steprate)) #Test overdisp. of n.b.
testZeroInflation(simulateResiduals(mac_c_nb_steprate)) #Test zero-inflation
testUniformity(simulateResiduals(mac_c_nb_steprate)) #Indicates heterosked.
testTemporalAutocorrelation(simulateResiduals(mac_c_nb_steprate), 
                            time=Macqmra$Week) #D-W test autocorr. at lag = 1

#Extract summary model information
summary(mac_c_nb_steprate)
mac_c_nb_steprate$theta
mac_c_nb_steprate$SE.theta
mac_c_nb_steprate$aic

#GLARMA modelling not required due to lack of AR or MA signatures.

#Test additional covariates

#Colour
#Recode "<2" as 2, and convert to numeric
Macqmra$Col400cont = Macqmra$Col400
Macqmra$Col400cont[Macqmra$Col400cont == "<2"] = 2
Macqmra$Col400cont = as.numeric(Macqmra$Col400cont)
mac_c_nb_steprate_colour=glm.nb(Macqmra$S1_N ~ Macqmra$Week*Macqmra$Phase+
                                  Macqmra$Col400cont+
                                  offset(log(Macqmra$S1_offset)))
coeftest(mac_c_nb_steprate_colour, vcov = vcov(mac_c_nb_steprate_colour))
coefci(mac_c_nb_steprate_colour,vcov = vcov(mac_c_nb_steprate_colour))
mac_c_nb_steprate_negbinirr = negbinirr(formula = 
                                          mac_c_nb_steprate_colour$call, 
                                        data = Macqmra, 
                                        robust = FALSE)
mac_c_nb_steprate_negbinirr
mac_c_nb_steprate_colour$aic
mac_c_nb_steprate_colour$theta
mac_c_nb_steprate_colour$SE.theta

#Temperature
mac_c_nb_steprate_temp=glm.nb(Macqmra$S1_N ~ Macqmra$Week*Macqmra$Phase+
                                Macqmra$Temp+
                                offset(log(Macqmra$S1_offset)))
coeftest(mac_c_nb_steprate_temp, vcov = vcov(mac_c_nb_steprate_temp))
coefci(mac_c_nb_steprate_temp,vcov = vcov(mac_c_nb_steprate_temp))
mac_c_nb_steprate_negbinirr = negbinirr(formula = mac_c_nb_steprate_temp$call,
                                        data = Macqmra, 
                                        robust = FALSE)
mac_c_nb_steprate_negbinirr
mac_c_nb_steprate_temp$aic
mac_c_nb_steprate_temp$theta
mac_c_nb_steprate_temp$SE.theta

#Turbidity
mac_c_nb_steprate_turbidity=glm.nb(Macqmra$S1_N ~ Macqmra$Week*Macqmra$Phase+
                                     Macqmra$Turbidity+
                                     offset(log(Macqmra$S1_offset)))
coeftest(mac_c_nb_steprate_turbidity, 
         vcov = vcov(mac_c_nb_steprate_turbidity))
coefci(mac_c_nb_steprate_turbidity,
       vcov = vcov(mac_c_nb_steprate_turbidity))
mac_c_nb_steprate_negbinirr = negbinirr(formula = 
                                          mac_c_nb_steprate_turbidity$call, 
                                        data = Macqmra, 
                                        robust = FALSE)
mac_c_nb_steprate_negbinirr
mac_c_nb_steprate_turbidity$aic
mac_c_nb_steprate_turbidity$theta
mac_c_nb_steprate_turbidity$SE.theta


#============ Mac - Crypto - ZINB - phase model =================

#Fit zero-inflated neg bin
mac_c_zinb_step = zeroinfl(Macqmra$S1_N ~ Macqmra$Week+Macqmra$Phase, 
                           offset=offset(log(Macqmra$S1_offset)), 
                           dist="negbin")

#Plot results
plot(mac_c_zinb_step$fitted.values/Macqmra$S1_offset,type="l")

#Extract model summary information
summary(mac_c_zinb_step)
extractAIC(mac_c_zinb_step)

#plot ACF and PACF to check for AR and MA processes
Acf(mac_c_zinb_step$residuals, lag.max=NULL, type=c("correlation"),
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(mac_c_zinb_step$residuals, lag.max=NULL, plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)


#============ Mac - Crypto - ZINB - phase and interact model ============

#Fit zero-inflated neg bin
mac_c_zinb_steprate = zeroinfl(Macqmra$S1_N ~ Macqmra$Week*Macqmra$Phase, 
                               offset=offset(log(Macqmra$S1_offset)), 
                               dist="negbin")

#Plot results
plot(mac_c_zinb_steprate$fitted.values/Macqmra$S1_offset,type="l")

#Extract model summary information
summary(mac_c_zinb_steprate)
extractAIC(mac_c_zinb_steprate)

#Plot ACF and PACF to check for ARMA processes
Acf(mac_c_zinb_steprate$residuals, lag.max=NULL, type=c("correlation"),
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(mac_c_zinb_steprate$residuals, lag.max=NULL,plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)


#============ Mac - Crypto - Poisson - phase model ===========

#Fit Poisson
mac_c_p_step=glm(Macqmra$S1_N ~ Macqmra$Week+Macqmra$Phase, 
                 offset=offset(log(Macqmra$S1_offset)),
                 family=poisson)

#Plot results
plot(mac_c_p_step$fitted.values/Macqmra$S1_offset, type="l")

#Extract model summary information
summary(mac_c_p_step)

#Plot ACF and PACF to check for ARMA processes
Acf(qresiduals(mac_c_p_step))
Pacf(qresiduals(mac_c_p_step))
hist(qresiduals(mac_c_p_step))
plot(simulateResiduals(mac_c_p_step))

#Perform tests
testDispersion(simulateResiduals(mac_c_p_step)) #Test overdisp. of Poiss.
testZeroInflation(simulateResiduals(mac_c_p_step)) #Test zero-inflation
testUniformity(simulateResiduals(mac_c_p_step)) #Detect heteroskedacticity
testTemporalAutocorrelation(simulateResiduals(mac_c_p_step), 
                            time=Macqmra$Week) #D-W test autocorr. at lag = 1


#============ Mac - Crypto - Poisson - phase and interact model =========

#Fit Poisson
mac_c_p_steprate=glm(Macqmra$S1_N ~ Macqmra$Week*Macqmra$Phase, 
                     offset=offset(log(Macqmra$S1_offset)),
                     family=poisson)

#Plot results
plot(mac_c_p_steprate$fitted.values/Macqmra$S1_offset, type="l")

#Extract model summary information
summary(mac_c_p_steprate)

#Plot ACF and PACF to check for ARMA processes
Acf(qresiduals(mac_c_p_steprate))
Pacf(qresiduals(mac_c_p_steprate))
hist(qresiduals(mac_c_p_steprate))
plot(simulateResiduals(mac_c_p_steprate)) 
testDispersion(simulateResiduals(mac_c_p_steprate)) #Test overdisp. of Poiss.
testZeroInflation(simulateResiduals(mac_c_p_steprate)) #Test zero-inflation
testUniformity(simulateResiduals(mac_c_p_steprate)) #Detect heteroskedacticity
testTemporalAutocorrelation(simulateResiduals(mac_c_p_steprate), 
                            time=Macqmra$Week) #D-W test autocorr. at lag = 1


#============ Mac - Crypto - ZIP - phase model ================
#Fit ZIP
mac_c_zip_step = zeroinfl(Macqmra$S1_N ~ Macqmra$Week+Macqmra$Phase, 
                          offset=offset(log(Macqmra$S1_offset)), 
                          dist="poisson")

#Plot results
plot(mac_c_zip_step$fitted.values/Macqmra$S1_offset,type="l")

#Extract model summary information
summary(mac_c_zip_step)
extractAIC(mac_c_zip_step)

#Plot ACF and PACF to check for AR and MA processes
Acf(mac_c_zip_step$residuals, lag.max=NULL, type=c("correlation"),
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(mac_c_zip_step$residuals, lag.max=NULL, plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)


#============ Mac - Crypto - ZIP - phase and interact model ===========

#Fit ZIP
mac_c_zip_steprate = zeroinfl(Macqmra$S1_N ~ Macqmra$Week*Macqmra$Phase, 
                              offset=offset(log(Macqmra$S1_offset)), 
                              dist="poisson")

#Plot results
plot(mac_c_zip_steprate$fitted.values/Macqmra$S1_offset,type="l")

#Extract model summary information
summary(mac_c_zip_steprate)
extractAIC(mac_c_zip_steprate)

#Plot ACF and PACF to check for AR and MA processes
Acf(mac_c_zip_steprate$residuals, lag.max=NULL, type=c("correlation"),
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(mac_c_zip_steprate$residuals, lag.max=NULL, plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)
