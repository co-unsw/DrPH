#Generate, plot, and test regression models for North Richmond Giardia

#Install and load necessary packages

install.packages("pscl")
install.packages("forecast")
install.packages("glarma")
install.packages("statmod")
install.packages("DHARMa")
install.packages("lmtest")
install.packages("sandwich")
install.packages("aod")
install.packages("mfx")
library(pscl)
library(forecast)
library(glarma)
library(statmod)
library(DHARMa)
library(lmtest)
library(sandwich)
library(aod)
library(mfx)

#Generate the recovery rate variable 
NRqmra$rateG<-NRqmra$S1_NG/NRqmra$S1_offsetG

#============ NR - Giardia - neg bin - phase model =================

#Fit neg bin model
nr_g_nb_step=glm.nb(NRqmra$S1_NG ~
                      NRqmra$Week+NRqmra$Phase+
                      offset(log(NRqmra$S1_offsetG)))

#Predicted values of the recovery rate
nr_g_nb_step_rate_pred=nr_g_nb_step$fitted.values/NRqmra$S1_offsetG
plot(NRqmra$Week, nr_g_nb_step_rate_pred, type="l", col="blue", lwd=2)

#Test for overall uniformity of the fitted neg bin model 
#(KS test) which can detect heteroskedacticity
testUniformity(simulateResiduals(nr_g_nb_step))

# ##Robust standard error - beta and IRR. 
# #Not indicated due to lack of heteroskedacity.
# #Robust standard errors for coefficients (beta) (and 95% CI).
# #Applies a correction factor to achieve heteroskedacticity-consistent 
# #estimate of robust standard error.
# coeftest(nr_g_nb_step, vcov=vcovHC(nr_g_nb_step, type="HC0"))
# coefci(nr_g_nb_step,vcov=vcovHC(nr_g_nb_step, type="HC0"))
# #IRR, IRR SE, z, and p values, 
# #this method is necessary to find these using robust standard error.
# nr_g_nb_step_negbinirr=negbinirr(formula=nr_g_nb_step$call, 
#                                    data=Macqmra, 
#                                    robust=TRUE)
# nr_g_nb_step_negbinirr

###Non-robust standard error - beta and IRR:
#Non-robust standard error for coefficients (beta) (and 95% CI):
#Indicated due to homoskedacticity
coeftest(nr_g_nb_step, vcov=vcov(nr_g_nb_step))
coefci(nr_g_nb_step,vcov=vcov(nr_g_nb_step))
#Non-robust IRR, IRR SE, z, and p values:
nr_g_nb_step_negbinirr=negbinirr(formula=nr_g_nb_step$call,
                                 data=NRqmra, 
                                 robust=FALSE)
nr_g_nb_step_negbinirr

#plot ACF and PACF to check for AR and MA processes
Acf(qresiduals(nr_g_nb_step))
Pacf(qresiduals(nr_g_nb_step))
hist(qresiduals(nr_g_nb_step))
plot(simulateResiduals(nr_g_nb_step)) 

#Perform tests
odTest(nr_g_nb_step) #Test overdispersion (H0=Poisson; H1=overdispersed)
testDispersion(simulateResiduals(nr_g_nb_step)) #Test overdispersion of n.b.
testZeroInflation(simulateResiduals(nr_g_nb_step)) #Test zero-inflation
testTemporalAutocorrelation(simulateResiduals(nr_g_nb_step), 
                            time=NRqmra$Week) #D-W test autocorr. at lag=1

#Extract model summary information
summary(nr_g_nb_step)
nr_g_nb_step$theta
nr_g_nb_step$SE.theta
nr_g_nb_step$aic

#GLARMA modelling needed - AR or MA present
Xstep=as.matrix(cbind("Intercept"=1, 
                      "Week"=NRqmra$Week, 
                      "Phase"=NRqmra$Phase))
nr_g_nb_step_glarma=glarma(NRqmra$S1_NG, 
                           X=Xstep, 
                           offset=log(NRqmra$S1_offsetG), 
                           type="NegBin", 
                           thetaLags=1) #Most parsimonious model shown
plot(nr_g_nb_step_glarma)

# #Robust standard errors for coefficients (beta) (and 95% CI)
# #Not indicated due to lack of heteroskedacity
# coeftest(nr_g_nb_step, vcov=vcovHC(nr_g_nb_step, type="HC0"))
# coefci(nr_g_nb_step,vcov=vcovHC(nr_g_nb_step, type="HC0"))
# #IRR, IRR SE, z, and p values
# #this method is necessary to find these using robust standard error
# nr_g_nb_steprate_negbinirr=negbinirr(formula=nr_g_nb_step$call, 
#                                        data=Macqmra, 
#                                        robust=TRUE)
# nr_g_nb_steprate_negbinirr

#Plot ACF and PACF to check ARMA processes are controlled
Acf(normRandPIT(nr_g_nb_step_glarma)$rt)
Pacf(normRandPIT(nr_g_nb_step_glarma)$rt)
hist(normRandPIT(nr_g_nb_step_glarma)$rt)

#Extract model summary information
summary.glarma(nr_g_nb_step_glarma, tests=FALSE)

#Test additional covariates

#Colour
#Recode "<2" as 2, and convert to numeric
NRqmra$Col400cont=NRqmra$Col400
NRqmra$Col400cont[NRqmra$Col400cont=="<2"]=2
NRqmra$Col400cont=as.numeric(NRqmra$Col400cont)

Xstep=as.matrix(cbind("Intercept"=1, 
                      "Week"=NRqmra$Week, 
                      "Phase"=NRqmra$Phase, 
                      "Colour"=NRqmra$Col400cont))
nr_g_nb_step_glarma_colour=glarma(NRqmra$S1_NG, 
                                  X=Xstep, 
                                  offset=log(NRqmra$S1_offsetG),
                                  type="NegBin", 
                                  thetaLags=1)
summary.glarma(nr_g_nb_step_glarma_colour, tests=FALSE)

#Temperature
Xstep=as.matrix(cbind("Intercept"=1, 
                      "Week"=NRqmra$Week, 
                      "Phase"=NRqmra$Phase, 
                      "Temp"=NRqmra$Temp))
nr_g_nb_step_glarma_temp=glarma(NRqmra$S1_NG, 
                                X=Xstep, 
                                offset=log(NRqmra$S1_offsetG), 
                                type="NegBin", 
                                thetaLags=1) #No reasonable convergence
summary.glarma(nr_g_nb_step_glarma_temp, tests=FALSE)

#Turbidity
Xstep=as.matrix(cbind("Intercept"=1,
                      "Week"=NRqmra$Week, 
                      "Phase"=NRqmra$Phase, 
                      "Turbidity"=NRqmra$Turbidity))
nr_g_nb_step_glarma_turbidity=glarma(NRqmra$S1_NG, 
                                     X=Xstep, 
                                     offset=log(NRqmra$S1_offsetG), 
                                     type="NegBin", 
                                     thetaLags=1)
summary.glarma(nr_g_nb_step_glarma_turbidity, tests=FALSE)

#Turbidity and colour
Xstep=as.matrix(cbind("Intercept"=1,
                      "Week"=NRqmra$Week, 
                      "Phase"=NRqmra$Phase, 
                      "Colour"=NRqmra$Col400cont, 
                      "Turbidity"=NRqmra$Turbidity))
nr_g_nb_step_glarma_colturb=glarma(NRqmra$S1_NG,
                                   X=Xstep, 
                                   offset=log(NRqmra$S1_offsetG), 
                                   type="NegBin", 
                                   thetaLags=1)
summary.glarma(nr_g_nb_step_glarma_colturb, tests=FALSE)


#============ NR - Giardia - neg bin - phase and interact model ============

#Fit neg bin
nr_g_nb_steprate=glm.nb(NRqmra$S1_NG ~ 
                          NRqmra$Week*NRqmra$Phase+
                          offset(log(NRqmra$S1_offsetG)))

#Predicted values of the recovery rate
nr_g_nb_steprate_rate_pred=nr_g_nb_step$fitted.values/NRqmra$S1_offsetG
plot(NRqmra$Week, NRqmra$rateG, type="l")
lines(NRqmra$Week, nr_g_nb_steprate_rate_pred, type="l", col="blue", lwd=2)

#Test for overall uniformity of the fitted neg bin model 
#(KS test) which can detect heteroskedacticity
testUniformity(simulateResiduals(nr_g_nb_steprate))

# #Robust standard errors for coefficients (beta) (and 95% CI)
# #Not indicated due to lack of heteroskedacity
# coeftest(nr_g_nb_steprate, vcov=vcovHC(nr_g_nb_steprate, type="HC0"))
# coefci(nr_g_nb_steprate,vcov=vcovHC(nr_g_nb_steprate, type="HC0"))
# #IRR, IRR SE, z, and p values
# #this method is necessary to find these using robust standard error
# nr_g_nb_steprate_negbinirr=negbinirr(formula=nr_g_nb_steprate$call, 
#                                      data=Macqmra, 
#                                      robust=TRUE)
# nr_g_nb_steprate_negbinirr

#Non-robust standard error for coefficients (beta) (and 95% CI)
coeftest(nr_g_nb_steprate, vcov=vcov(nr_g_nb_steprate))
coefci(nr_g_nb_steprate,vcov=vcov(nr_g_nb_steprate))
#Non-robust IRR, IRR SE, z, and p values:
nr_g_nb_steprate_negbinirr=negbinirr(formula=nr_g_nb_steprate$call, 
                                     data=NRqmra, 
                                     robust=FALSE)
nr_g_nb_steprate_negbinirr

#Plot ACF and PACF to check for AR and MA processes
Acf(qresiduals(nr_g_nb_steprate))
Pacf(qresiduals(nr_g_nb_steprate))
hist(qresiduals(nr_g_nb_steprate))
plot(simulateResiduals(nr_g_nb_steprate)) 
odTest(nr_g_nb_steprate) #Test overdispersion (H0=Poiss.; H1=overdisp.)
testDispersion(simulateResiduals(nr_g_nb_steprate)) #Test overdisp. of n.b.
testZeroInflation(simulateResiduals(nr_g_nb_steprate)) #Test zero-inflation
testTemporalAutocorrelation(simulateResiduals(nr_g_nb_steprate), 
                            time=NRqmra$Week) #D-W test autocorr. at lag=1

#Extract model summary information
summary(nr_g_nb_steprate)
nr_g_nb_steprate$theta
nr_g_nb_steprate$SE.theta
nr_g_nb_steprate$aic

#GLARMA modelling needed - AR or MA present
Xsteprate=as.matrix(cbind("Intercept"=1, 
                          "Week"=NRqmra$Week, 
                          "Phase"=NRqmra$Phase, 
                          "Interact"=NRqmra$Interact))
nr_g_nb_steprate_glarma=glarma(y=NRqmra$S1_NG, 
                               X=Xsteprate, 
                               offset=log(NRqmra$S1_offsetG), 
                               type="NegBin", 
                               phiLags=1) #Most parsimonious
plot(nr_g_nb_steprate_glarma)

#Plot ACF and PACF to check ARMA processes are controlled
Acf(normRandPIT(nr_g_nb_steprate_glarma)$rt)
Pacf(normRandPIT(nr_g_nb_steprate_glarma)$rt)
summary.glarma(nr_g_nb_steprate_glarma, tests=FALSE)

#Extract model summary information
extractAIC(nr_g_nb_steprate_glarma)
likTests(nr_g_nb_steprate_glarma)#Throws exception related to offset
glarmaModelEstimates(nr_g_nb_steprate_glarma)

#Need to do Wald or L-R test to compare ARMA model with GLM,
#however, likTests() function throwing unknown exception relating to offset.
#The following code is extracted from the likTests function, glarma package.
#Because the H0 model is already fitted (i.e. nr_g_nb_steprate),
#the same function is achieved with the following code.

#Check errors
nr_g_nb_steprate_glarma$errCode
nr_g_nb_steprate_glarma$WError

####START of code adapted from likTests(), glarma package####
ll.null <- nr_g_nb_steprate_glarma$r - nr_g_nb_steprate$aic/2
ll.alt <- nr_g_nb_steprate_glarma$logLik
pq <- nr_g_nb_steprate_glarma$pq

#L-R test for whether glarma is more suitable than glm
LRT <- 2 * (ll.alt - ll.null)
LRT.P <- 1 - pchisq(LRT, pq)
#Print the LR test statistic: H0=neg bin glm; H1 indicates neg bin GLARMA
LRT 
#Print the LR test p-value
#to confirm stat. sig. difference of neg bin GLARMA over neg bin glm
LRT.P 

#Wald test for whether glarma is more suitable than glm
thetahat <- nr_g_nb_steprate_glarma$delta[nr_g_nb_steprate_glarma$r + 1:pq]
Wald <- (thetahat) %*%
  solve(nr_g_nb_steprate_glarma$cov[nr_g_nb_steprate_glarma$r + 1:pq, 
                                    nr_g_nb_steprate_glarma$r + 1:pq]) %*%
  thetahat
Wald.P <- 1 - pchisq(Wald, pq)
#Print the Wald test statistic: H0=neg bin glm; H1 indicates neg bin GLARMA
Wald 
#Print the Wald test p-value
#to confirm stat. sig. difference of neg bin GLARMA over neg bin glm
Wald.P 
####END of code adapted from likTests(), glarma package####

#Wald test for coefficients (beta)
#This is different Wald test to the above. 
#It is for the p-value for the stat. sig. of the regressors
#(2=Week; 3=Phase; 4=Interact)
print(wald.test(Sigma=nr_g_nb_steprate_glarma$cov, 
                b=glarmaModelEstimates(nr_g_nb_steprate_glarma)[,1], 
                Terms=2, 
                verbose=TRUE), digits=7)
print(wald.test(Sigma=nr_g_nb_steprate_glarma$cov,
                b=glarmaModelEstimates(nr_g_nb_steprate_glarma)[,1], 
                Terms=3, 
                verbose=TRUE), digits=7)
print(wald.test(Sigma=nr_g_nb_steprate_glarma$cov, 
                b=glarmaModelEstimates(nr_g_nb_steprate_glarma)[,1], 
                Terms=4, 
                verbose=TRUE), digits=7)


#============ NR - Giardia - ZINB - phase model =================

#Fit zero-inflated neg bin
nr_g_zinb_step=zeroinfl(NRqmra$S1_NG ~ NRqmra$Week+NRqmra$Phase, 
                        offset=offset(log(NRqmra$S1_offsetG)), 
                        dist="negbin")

#Plot results
plot(nr_g_zinb_step$fitted.values/NRqmra$S1_offsetG,type="l")

#Extract model summary information
summary(nr_g_zinb_step)
extractAIC(nr_g_zinb_step)

#plot ACF and PACF to check for AR and MA processes
qqrplot(nr_g_zinb_step)
Acf(nr_g_zinb_step$residuals,
    lag.max=NULL, type=c("correlation"),
    plot=TRUE, na.action=na.contiguous, demean=TRUE)
Pacf(nr_g_zinb_step$residuals,
     lag.max=NULL, plot=TRUE, 
     na.action=na.contiguous, demean=TRUE)


#============ NR - Giardia - ZINB - phase and interact model ============
#Fit zero-inflated neg bin
nr_g_zinb_steprate=zeroinfl(NRqmra$S1_NG ~ NRqmra$Week*NRqmra$Phase, 
                            offset=offset(log(NRqmra$S1_offsetG)), 
                            dist="negbin")
#Plot results
plot(nr_g_zinb_steprate$fitted.values/NRqmra$S1_offsetG,type="l")

#Extract model summary information
summary(nr_g_zinb_steprate)
extractAIC(nr_g_zinb_steprate)

#Plot ACF and PACF to check for AR and MA processes
Acf(nr_g_zinb_steprate$residuals,
    lag.max=NULL, type=c("correlation"),
    plot=TRUE, na.action=na.contiguous, demean=TRUE)
Pacf(nr_g_zinb_steprate$residuals,
     lag.max=NULL, plot=TRUE, na.action=na.contiguous, demean=TRUE)

#Zero-inflation model component non-significant.
#Glarma modelling of ZINB not available.

#============ NR - Giardia - Poisson - phase model ===========

#Fit Poisson
nr_g_p_step=glm(NRqmra$S1_NG ~ NRqmra$Week+NRqmra$Phase, 
                offset=offset(log(NRqmra$S1_offsetG)),
                family=poisson)

#Plot results
plot(nr_g_p_step$fitted.values/NRqmra$S1_offsetG, type="l")

#Extract model summary information
summary(nr_g_p_step)

#Plot ACF and PACF to check for AR and MA processes
Acf(qresiduals(nr_g_p_step))
Pacf(qresiduals(nr_g_p_step))
hist(qresiduals(nr_g_p_step))
plot(simulateResiduals(nr_g_p_step)) 

#Perform tests
testDispersion(simulateResiduals(nr_g_p_step)) #Test overdisp. of Poisson
testZeroInflation(simulateResiduals(nr_g_p_step)) #Test zero-inflation
testUniformity(simulateResiduals(nr_g_p_step)) #Detect heteroskedacticity
testTemporalAutocorrelation(simulateResiduals(nr_g_p_step), 
                            time=Macqmra$Week) #D-W test autocorr. at lag=1

#Overdispersion was detected. The Poisson approach would be abandoned.
#If GLARMA modelling was needed, the following would be used (as an example)
Xstep=as.matrix(cbind("Intercept"=1, 
                      "Week"=NRqmra$Week, 
                      "Phase"=NRqmra$Phase))
nr_g_p_step_glarma=glarma(NRqmra$S1_NG, 
                          X=Xstep, 
                          offset=log(NRqmra$S1_offsetG), 
                          type="Poi", 
                          phiLags=1, thetaLags=2) #Most parsimonious

#Plot ACF and PACF to check ARMA processes are controlled
Acf(normRandPIT(nr_g_p_step_glarma)$rt)
Pacf(normRandPIT(nr_g_p_step_glarma)$rt)

#Extract model summary information
summary.glarma(nr_g_p_step_glarma, tests=FALSE)


#============ NR - Giardia - Poisson - phase and interact model =========

#Fit Poisson
nr_g_p_steprate=glm(NRqmra$S1_NG ~ NRqmra$Week*NRqmra$Phase, 
                    offset=offset(log(NRqmra$S1_offsetG)),
                    family=poisson)

#Plot the results
plot(nr_g_p_steprate$fitted.values/NRqmra$S1_offsetG, type="l")

summary(nr_g_p_steprate)

#Plot ACF and PACF to check for AR and MA processes
Acf(qresiduals(nr_g_p_steprate))
Pacf(qresiduals(nr_g_p_steprate))
hist(qresiduals(nr_g_p_steprate))
plot(simulateResiduals(nr_g_p_steprate))

#Perform tests
testDispersion(simulateResiduals(nr_g_p_steprate)) #Test overdisp. of Poiss.
testZeroInflation(simulateResiduals(nr_g_p_steprate)) #Test zero-inflation
testUniformity(simulateResiduals(nr_g_p_steprate)) #Detect heteroskedacticity
testTemporalAutocorrelation(simulateResiduals(nr_g_p_steprate), 
                            time=Macqmra$Week) #D-W test autocorr. at lag=1

#Overdispersion was detected. The Poisson approach would be abandoned.


#============ NR - Giardia - ZIP - phase model ================

#Fit Zero-inflated Poisson
nr_g_zip_step=zeroinfl(NRqmra$S1_NG ~ NRqmra$Week+NRqmra$Phase, 
                       offset=offset(log(NRqmra$S1_offsetG)), 
                       dist="poisson")

#Plot results
plot(nr_g_zip_step$fitted.values/NRqmra$S1_offsetG, type="l")

#Extract summary model information
summary(nr_g_zip_step)
extractAIC(nr_g_zip_step)

#Plot ACF and PACF to check for AR and MA processes
Acf(nr_g_zip_step$residuals,
    lag.max=NULL, type=c("correlation"),
    plot=TRUE, na.action=na.contiguous, demean=TRUE)
Pacf(nr_g_zip_step$residuals, lag.max=NULL,
     plot=TRUE, na.action=na.contiguous, demean=TRUE)

#Glarma modelling of ZINB not available.


#============ NR - Giardia - ZIP - phase and interact model ===========

#Fit Zero-inflated Poisson
nr_g_zip_steprate=zeroinfl(NRqmra$S1_NG ~ NRqmra$Week*NRqmra$Phase, 
                           offset=offset(log(NRqmra$S1_offsetG)), 
                           dist="poisson")

#Plot results
plot(nr_g_zip_steprate$fitted.values/NRqmra$S1_offsetG,type="l")

#Extract model summary information
summary(nr_g_zip_steprate)
extractAIC(nr_g_zip_steprate)

#plot ACF and PACF to check for AR and MA processes
Acf(nr_g_zip_steprate$residuals, 
    lag.max=NULL,type=c("correlation"), plot=TRUE, 
    na.action=na.contiguous, demean=TRUE)
Pacf(nr_g_zip_steprate$residuals, lag.max=NULL,
     plot=TRUE, na.action=na.contiguous, demean=TRUE)

#Glarma modelling of ZINB not available.