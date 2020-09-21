#Generate, plot, and test regression models for North Richmond Cryptosporidium

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
NRqmra$rate<-NRqmra$S1_N/NRqmra$S1_offset

#============ NR - Crypto - neg bin - phase model =================

#Fit neg bin
nr_c_nb_step=glm.nb(NRqmra$S1_N ~ NRqmra$Week+NRqmra$Phase+
                      offset(log(NRqmra$S1_offset)))

#Predicted values of the recovery rate
nr_c_nb_step_rate_pred = nr_c_nb_step$fitted.values/NRqmra$S1_offset
plot(NRqmra$Week, nr_c_nb_step_rate_pred, type="l", col="blue", lwd=2)

#Test for overall uniformity of the fitted neg bin model 
#(KS test) which can detect heteroskedacticity
testUniformity(simulateResiduals(nr_c_nb_step))

# ##Robust standard error - beta and IRR. 
# #Not indicated due to lack of heteroskedacity.
# #Robust standard errors for coefficients (beta) (and 95% CI).
# #Applies a correction factor to achieve heteroskedacticity-consistent 
# #estimate of robust standard error.
# coeftest(nr_c_nb_step, vcov = vcovHC(nr_c_nb_step, type="HC0"))
# coefci(nr_c_nb_step,vcov = vcovHC(nr_c_nb_step, type="HC0"))
# #IRR, IRR SE, z, and p values, 
# #this method is necessary to find these using robust standard error.
# nr_c_nb_step_negbinirr = negbinirr(formula = nr_c_nb_step$call,
#                                    data = NRqmra, 
#                                    robust = TRUE)
# nr_c_nb_step_negbinirr

#Non-robust standard error for coefficients (beta) (and 95% CI)
#Indicated due to homoskedacticity
coeftest(nr_c_nb_step, vcov = vcov(nr_c_nb_step))
coefci(nr_c_nb_step,vcov = vcov(nr_c_nb_step))
#Non-robust IRR, IRR SE, z, and p values:
nr_c_nb_step_negbinirr = negbinirr(formula = nr_c_nb_step$call, 
                                   data = NRqmra, 
                                   robust = FALSE)
nr_c_nb_step_negbinirr

#Plot ACF and PACF to check for AR and MA processes
Acf(qresiduals(nr_c_nb_step))
Pacf(qresiduals(nr_c_nb_step))
hist(qresiduals(nr_c_nb_step))
plot(simulateResiduals(nr_c_nb_step)) 

#Perform tests
odTest(nr_c_nb_step) #Test overdispersion (H0 = Poisson; H1 = overdispersed)
testDispersion(simulateResiduals(nr_c_nb_step)) #Test overdispersion of n.b.
testZeroInflation(simulateResiduals(nr_c_nb_step)) #Test zero-inflation
testTemporalAutocorrelation(simulateResiduals(nr_c_nb_step), 
                            time=NRqmra$Week) #D-W test autocorr. at lag = 1

#Extract model summary information
summary(nr_c_nb_step)
nr_c_nb_step$theta
nr_c_nb_step$SE.theta
nr_c_nb_step$aic

#Test additional covariates

#Colour
#Recode "<2" as 2, and convert to numeric
NRqmra$Col400cont = NRqmra$Col400
NRqmra$Col400cont[NRqmra$Col400cont == "<2"] = 2
NRqmra$Col400cont = as.numeric(NRqmra$Col400cont)
nr_c_nb_step_colour=glm.nb(NRqmra$S1_N ~ 
                             NRqmra$Week+
                             NRqmra$Phase+
                             NRqmra$Col400cont+
                             offset(log(NRqmra$S1_offset)))
coeftest(nr_c_nb_step_colour, vcov = vcov(nr_c_nb_step_colour))
coefci(nr_c_nb_step_colour,vcov = vcov(nr_c_nb_step_colour))
nr_c_nb_step_negbinirr = negbinirr(formula = nr_c_nb_step_colour$call, 
                                   data = NRqmra, 
                                   robust = FALSE)
nr_c_nb_step_negbinirr
nr_c_nb_step_colour$aic
nr_c_nb_step_colour$theta

#Turbidity
nr_c_nb_step_turbidity=glm.nb(NRqmra$S1_N ~ NRqmra$Week+
                                NRqmra$Phase+
                                NRqmra$Turbidity+
                                offset(log(NRqmra$S1_offset)))
coeftest(nr_c_nb_step_turbidity, vcov = vcov(nr_c_nb_step_turbidity))
coefci(nr_c_nb_step_turbidity,vcov = vcov(nr_c_nb_step_turbidity))
nr_c_nb_step_negbinirr = negbinirr(formula = nr_c_nb_step_turbidity$call, 
                                   data = NRqmra, 
                                   robust = FALSE)
nr_c_nb_step_negbinirr
nr_c_nb_step_turbidity$aic
nr_c_nb_step_turbidity$theta


#============ NR - Crypto - neg bin - phase and interact model ============

#Fit neg bin
nr_c_nb_steprate=glm.nb(NRqmra$S1_N ~ NRqmra$Week*NRqmra$Phase+
                          offset(log(NRqmra$S1_offset)))

#Predicted values of the recovery rate
nr_c_nb_steprate_rate_pred = nr_c_nb_steprate$fitted.values/NRqmra$S1_offset
plot(NRqmra$Week, nr_c_nb_steprate_rate_pred, type="l", col="blue", lwd=2)

#Test for overall uniformity of the fitted neg bin model 
#(KS test) which can detect heteroskedacticity
testUniformity(simulateResiduals(nr_c_nb_steprate))

# #Robust standard errors for coefficients (beta) (and 95% CI)
# #Not indicated due to lack of heteroskedacity
# coeftest(nr_c_nb_steprate, vcov = vcovHC(nr_c_nb_steprate, type="HC0"))
# coefci(nr_c_nb_steprate,vcov = vcovHC(nr_c_nb_steprate, type="HC0"))
# #IRR, IRR SE, z, and p values
# #this method is necessary to find these using robust standard error
# nr_c_nb_steprate_negbinirr = negbinirr(formula = nr_c_nb_steprate$call, 
#                                        data = NRqmra, 
#                                        robust = TRUE)
# nr_c_nb_steprate_negbinirr

#Non-robust standard error for coefficients (beta) (and 95% CI)
#Indicated due to homoskedacticity
coeftest(nr_c_nb_steprate, vcov = vcov(nr_c_nb_steprate))
coefci(nr_c_nb_steprate,vcov = vcov(nr_c_nb_steprate))
#Non-robust IRR, IRR SE, z, and p values:
nr_c_nb_step_negbinirr = negbinirr(formula = nr_c_nb_steprate$call, 
                                   data = NRqmra, 
                                   robust = FALSE)
nr_c_nb_step_negbinirr

#Plot ACF and PACF to check for AR and MA processes
Acf(qresiduals(nr_c_nb_steprate))
Pacf(qresiduals(nr_c_nb_steprate))
hist(qresiduals(nr_c_nb_steprate))
plot(simulateResiduals(nr_c_nb_steprate)) 
odTest(nr_c_nb_steprate) #Test overdispersion (H0 = Poiss.; H1 = overdisp.)
testDispersion(simulateResiduals(nr_c_nb_steprate)) #Test overdisp. of n.b.
testZeroInflation(simulateResiduals(nr_c_nb_steprate)) #Test zero-inflation
testTemporalAutocorrelation(simulateResiduals(nr_c_nb_steprate), 
                            time=NRqmra$Week) #D-W test autocorr. at lag = 1

summary(nr_c_nb_steprate)
nr_c_nb_steprate$theta
nr_c_nb_steprate$SE.theta
nr_c_nb_steprate$aic


#============ NR - Crypto - ZINB - phase model =================

#Fit zero-inflated neg bin
nr_c_zinb_step = zeroinfl(NRqmra$S1_N ~ NRqmra$Week+NRqmra$Phase, 
                          offset=offset(log(NRqmra$S1_offset)), 
                          dist="negbin")

#Plot results
plot(nr_c_zinb_step$fitted.values/NRqmra$S1_offset,type="l")
summary(nr_c_zinb_step)
extractAIC(nr_c_zinb_step)

#Plot ACF and PACF to check for AR and MA processes
Acf(nr_c_zinb_step$residuals, lag.max=NULL, type=c("correlation"), 
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(nr_c_zinb_step$residuals,lag.max=NULL, plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)


#============ NR - Crypto - ZINB - phase and interact model ============

#Fit zero-inflated neg bin
nr_c_zinb_steprate = zeroinfl(NRqmra$S1_N ~ NRqmra$Week*NRqmra$Phase, 
                              offset=offset(log(NRqmra$S1_offset)), 
                              dist="negbin")

#Plot results
plot(nr_c_zinb_steprate$fitted.values/NRqmra$S1_offset,type="l")

#Extract model summary information
summary(nr_c_zinb_steprate)
extractAIC(nr_c_zinb_steprate)

#Plot ACF and PACF to check for AR and MA processes
Acf(nr_c_zinb_steprate$residuals, lag.max=NULL,type=c("correlation"),
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(nr_c_zinb_steprate$residuals, lag.max=NULL,
     plot = TRUE, na.action = na.contiguous, demean = TRUE)


#============ NR - Crypto - Poisson - phase model ===========

#Fit Poisson
nr_c_p_step=glm(NRqmra$S1_N ~ NRqmra$Week+NRqmra$Phase, 
                offset=offset(log(NRqmra$S1_offset)),
                family=poisson)

#Plot the results
plot(nr_c_p_step$fitted.values/NRqmra$S1_offset, type="l")

#Extract model summary information
summary(nr_c_p_step)

#Plot ACF and PACF to check for AR and MA processes
Acf(nr_c_p_step$residuals, lag.max=NULL, type=c("correlation"), 
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(nr_c_p_step$residuals, lag.max=NULL, plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)


#============ NR - Crypto - Poisson - phase and interact model =========

#Fit Poisson
nr_c_p_steprate=glm(NRqmra$S1_N ~ NRqmra$Week*NRqmra$Phase, 
                    offset=offset(log(NRqmra$S1_offset)),
                    family=poisson)

#Plot the results
plot(nr_c_p_steprate$fitted.values/NRqmra$S1_offset, type="l")

#Extract model summary information
summary(nr_c_p_steprate)

#Plot ACF and PACF to check for AR and MA processes
Acf(nr_c_p_steprate$residuals, lag.max=NULL, type=c("correlation"),
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(nr_c_p_steprate$residuals, lag.max=NULL, plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)


#============ NR - Crypto - ZIP - phase model ================

#Fit ZIP
nr_c_zip_step = zeroinfl(NRqmra$S1_N ~ NRqmra$Week+NRqmra$Phase, 
                         offset=offset(log(NRqmra$S1_offset)),
                         dist="poisson")

#Plot results
plot(nr_c_zip_step$fitted.values/NRqmra$S1_offset,type="l")

#Extract model summary information
summary(nr_c_zip_step)
extractAIC(nr_c_zip_step)

#Plot ACF and PACF to check for AR and MA processes
Acf(nr_c_zip_step$residuals, lag.max=NULL, type=c("correlation"),
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(nr_c_zip_step$residuals, lag.max=NULL, plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)


#============ NR - Crypto - ZIP - phase and interact model ===========

#Fit ZIP
nr_c_zip_steprate = zeroinfl(NRqmra$S1_N ~ NRqmra$Week*NRqmra$Phase, 
                             offset=offset(log(NRqmra$S1_offset)), 
                             dist="poisson")

#Plot results
plot(nr_c_zip_steprate$fitted.values,type="l")

#Extract model summary information
summary(nr_c_zip_steprate)
extractAIC(nr_c_zip_steprate)

#Plot ACF and PACF to check for AR and MA processes
Acf(nr_c_zip_steprate$residuals, lag.max=NULL, type=c("correlation"),
    plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(nr_c_zip_steprate$residuals, lag.max=NULL, plot = TRUE, 
     na.action = na.contiguous, demean = TRUE)
