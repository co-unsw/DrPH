#Performs the deterministic QMRA on the preferred ITS models
#and plots the results against the tolerable (negligible) limits 
#for Pinf,a and DALY/p/y

#Install and load necessary packages
install.packages("ggplot2")
install.packages("viridis")
install.packages("scales")
install.packages("gtable")
install.packages("grid")
library(ggplot2)
library(viridis)
library(scales)
library(gtable)
library(grid)

#Set up global variables - QMRA assumptions
macLRV=3 #LRV achieved by direct filtration under given CCPs
nrLRV=3.5 #LRV achieved by conventional treatment under given CCPs
Ev=1 #Exposure of 1 L unheated water per person per day
Ef=365.25 #Exposure frequency per year
Hi=0.3 #WSAA assumption for human-infective fraction of total oocysts
P=1 #WSAA assumption for human-pathogenic fraction of infective oocysts
r_c=0.2 #dose-response for Cryptosporidium
r_g=0.01982 #dose-response for Giardia
Pill_c=0.7 #Pill|inf for Cryptosporidium
Pill_g=0.24 #Pill|inf for Giardia
Bc_c=2.7e-3 #DALY/case for Cryptosporidium
Bc_g=1.7e-3 #DALY/case for Giardia
Sf=1 #Susceptible fraction

#Trend code below is replicated for convenience
#Prepare the fitted trend vectors
#Macarthur Cryptosporidium
B0=mac_c_nb_step$coefficients[[1]] #Intercept
B1=mac_c_nb_step$coefficients[[2]] #Week
B2=mac_c_nb_step$coefficients[[3]] #Phase
mactrendc=exp(B0+(B1*Macqmra$Week)+(B2*Macqmra$Phase))

#Macarthur Giardia
B0=mac_g_nb_step_coltemp$coefficients[[1]] #Intercept
B1=mac_g_nb_step_coltemp$coefficients[[2]] #Week
B2=mac_g_nb_step_coltemp$coefficients[[3]] #Phase
mactrendg=exp(B0+(B1*Macqmra$Week)+(B2*Macqmra$Phase))

#North Richmond Cryptosporidium
B0=nr_c_nb_step$coefficients[[1]] #Intercept
B1=nr_c_nb_step$coefficients[[2]] #Week
B2=nr_c_nb_step$coefficients[[3]] #Phase
nrtrendc=exp(B0+(B1*NRqmra$Week)+(B2*NRqmra$Phase))

#North Richmond Giardia
glarmaB0=nr_g_nb_step_glarma$delta[1,1] #Intercept
glarmaB1=nr_g_nb_step_glarma$delta[2,1] #Week
glarmaB2=nr_g_nb_step_glarma$delta[3,1] #Phase
nrtrendg=exp(glarmaB0+(glarmaB1*NRqmra$Week)+(glarmaB2*NRqmra$Phase))

#Print to console the log-valued effect estimates for protozoal level
log10(mactrendc[366])-log10(mactrendc[365]) #Macarthur Cryptosporidium
log10(mactrendg[366])-log10(mactrendg[365]) #Macarthur Giardia
log10(nrtrendc[366])-log10(nrtrendc[365]) #North Richmond Cryptosporidium
log10(nrtrendg[366])-log10(nrtrendg[365]) #North Richmond Giardia

#QMRA
#Annual probability of infection (Pinf,a)
#Macarthur Cryptosporidium
mac_qmra_Pinfatrend_c=
  1-((1-(1-exp(-r_c*(10^(log10(mactrendc)-macLRV))*Hi*Ev)))^Ef) 
#Macarthur Giardia
mac_qmra_Pinfatrend_g=
  1-((1-(1-exp(-r_g*(10^(log10(mactrendg)-macLRV))*Hi*Ev)))^Ef) 
#North Richmond Cryptosporidium
nr_qmra_Pinfatrend_c=
  1-((1-(1-exp(-r_c*(10^(log10(nrtrendc)-nrLRV))*Hi*Ev)))^Ef) 
#North Richmond Giardia
nr_qmra_Pinfatrend_g=
  1-((1-(1-exp(-r_g*(10^(log10(nrtrendg)-nrLRV))*Hi*Ev)))^Ef) 

#Effect estimates, Pinf,a
#Macarthur Cryptosporidium
mac_pinfa_c_effect=toString(signif(mac_qmra_Pinfatrend_c[366]-
                                     mac_qmra_Pinfatrend_c[365], digits=3)) 
#Macarthur Giardia
mac_pinfa_g_effect=toString(signif(mac_qmra_Pinfatrend_g[366]-
                                       mac_qmra_Pinfatrend_g[365], digits=3))
#North Richmond Cryptosporidium
nr_pinfa_c_effect=toString(signif(nr_qmra_Pinfatrend_c[366]-
                                      nr_qmra_Pinfatrend_c[365], digits=3)) 
#North Richmond Giardia
nr_pinfa_g_effect=toString(signif(nr_qmra_Pinfatrend_g[366]-
                                      nr_qmra_Pinfatrend_g[365], digits=3)) 

#Disease burden (Ba)
#Macarthur Cryptosporidium
mac_qmra_Batrend_c=mac_qmra_Pinfatrend_c*Pill_c*Sf*Bc_c 
#Macarthur Giardia
mac_qmra_Batrend_g=mac_qmra_Pinfatrend_g*Pill_g*Sf*Bc_g 
#North Richmond Cryptosporidium
nr_qmra_Batrend_c=nr_qmra_Pinfatrend_c*Pill_c*Sf*Bc_c 
#North Richmond Giardia
nr_qmra_Batrend_g=nr_qmra_Pinfatrend_g*Pill_g*Sf*Bc_g 

#Effect estimates, Ba
#Macarthur Cryptosporidium
mac_ba_c_effect=toString(signif(mac_qmra_Batrend_c[366]-
                                    mac_qmra_Batrend_c[365], digits=3))
#Macarthur Giardia
mac_ba_g_effect=toString(signif(mac_qmra_Batrend_g[366]-
                                    mac_qmra_Batrend_g[365], digits=3)) 
#North Richmond Cryptosporidium
nr_ba_c_effect=toString(signif(nr_qmra_Batrend_c[366]-
                                   nr_qmra_Batrend_c[365], digits=3)) 
#North Richmond Giardia
nr_ba_g_effect=toString(signif(nr_qmra_Batrend_g[366]-
                                   nr_qmra_Batrend_g[365], digits=3)) 

#Plot the figures
#Construct the figure for annual probability of infection
figqmra1=ggplot(Macqmra, aes(y=S1_V, x=Week))+
  
  #Trends
  #Pinf,a North Richmond Giardia
  geom_line(aes(y=nr_qmra_Pinfatrend_g,
                colour="North Richmond Giardia"), size=1)+
  #Pinf,a North Richmond Cryptosporidium
  geom_line(aes(y=nr_qmra_Pinfatrend_c, 
                colour="North Richmond Cryptosporidium"), size=1)+
  #Pinf,a Macarthur Cryptosporidium
  geom_line(aes(y=mac_qmra_Pinfatrend_c, 
                colour="Macarthur Cryptosporidium"), 
            linetype="dashed", size=1)+
  #Pinf,a Macarthur Giardia
  geom_line(aes(y=mac_qmra_Pinfatrend_g, 
                colour="Macarthur Giardia"), 
            linetype="dashed", size=1)+
  
  #Intervention and tolerable risk lines
  geom_vline(aes(xintercept=366), linetype="dotted", size=0.75)+
  geom_hline(aes(yintercept=10^-4),  size=1, colour="red")+#Tolerable Pinf,a
  
  #Effect estimates (x,y are nudged to improve readability)
  geom_text(x=366-4,
            y=mac_qmra_Pinfatrend_c[366]*1.05, 
            aes(colour="Macarthur Cryptosporidium", 
                label=paste("EE: ",mac_pinfa_c_effect)), 
            size=3, hjust=1, show.legend=FALSE)+
  geom_text(x=366-4,
            y=mac_qmra_Pinfatrend_g[366]*0.85, 
            aes(colour="Macarthur Giardia", 
                label=paste("EE: ",mac_pinfa_g_effect)), 
            size=3, hjust=1, show.legend=FALSE)+
  geom_text(x=366-4,
            y=nr_qmra_Pinfatrend_c[366]*1.05, 
            aes(colour="North Richmond Cryptosporidium", 
                label=paste("EE: ",nr_pinfa_c_effect)), 
            size=3, hjust=1, show.legend=FALSE)+
  geom_text(x=366-4,
            y=nr_qmra_Pinfatrend_g[366]*1.05, 
            aes(colour="North Richmond Giardia", 
                label=paste("EE: ",nr_pinfa_g_effect)), 
            size=3, hjust=1, show.legend=FALSE)+
  
  #Set up logarithmic y-axis
  scale_y_continuous(expand=c(0,0), 
                     breaks=trans_breaks("log10", function(x) 10^x), 
                     labels=trans_format("log10", math_format(10^.x)))+
  coord_trans(y="log10")+
  
  #Set up x-axis
  scale_x_continuous(expand=c(0,0))+
  expand_limits(x=c(0,494), y=c(10^-3.5,10^-5))+
  labs(x="Time (week)", 
       y="Probability of infection\n(person\U207B\U00B9 year\U207B\u00B9)")+
  
  #Labels for intervention and tolerable risk lines
  annotate("text", x=370, y=10^-3.55, 
           label="Intervention", hjust=0)+
  annotate("text", x=10, y=10^-3.945, 
           label="Reference level of risk", hjust=0)+

  #Apply aesthetics
  scale_color_viridis_d(labels=c(
    expression(paste("Macarthur ", italic("Cryptosporidium"), " risk")), 
    expression(paste("Macarthur ", italic("Giardia"), " risk")), 
    expression(paste("North Richmond ", italic("Cryptosporidium"), " risk")), 
    expression(paste("North Richmond ", italic("Giardia"), " risk"))
  ))+scale_fill_viridis_d()+theme_bw()+
  guides(colour=guide_legend(reverse=TRUE, nrow=2, 
                               override.aes=list(linetype=c("solid", 
                                                                "solid",
                                                                "dashed", 
                                                                "dashed"))))+
  theme(legend.title=element_blank(), 
        legend.position="bottom", 
        legend.text.align=0)

figqmra1 #Print the figure

#Construct the figure for annual disease burden
figqmra2=ggplot(Macqmra, aes(y=S1_V, x=Week))+
  
  #Trends
  #Ba North Richmond Cryptosporidium
  geom_line(aes(y=nr_qmra_Batrend_c, 
                colour="North Richmond Cryptosporidium"), size=1)+
  #Ba North Richmond Giardia  
  geom_line(aes(y=nr_qmra_Batrend_g, 
                colour="North Richmond Giardia"), size=1)+
  #Ba Macarthur Cryptosporidium
  geom_line(aes(y=mac_qmra_Batrend_c, 
                colour="Macarthur Cryptosporidium"), 
            linetype="dashed", size=1)+
  #Ba Macarthur Giardia
  geom_line(aes(y=mac_qmra_Batrend_g, colour="Macarthur Giardia"), 
            linetype="dashed", size=1)+

  #Intervention and tolerable risk lines
  geom_vline(aes(xintercept=366), linetype="dotted", size=0.75)+
  geom_hline(aes(yintercept=10^-6), size=1, colour="red")+#Tolerable Ba
  
  #Effect estimates (x,y are nudged to improve readability)
  geom_text(x=366-2,
            y=mac_qmra_Batrend_c[366]*1.15, 
            aes(colour="Macarthur Cryptosporidium", 
                label=paste("EE: ",mac_ba_c_effect)), 
            size=3, hjust=1, show.legend=FALSE)+
  geom_text(x=366-2,
            y=mac_qmra_Batrend_g[366]*0.8, 
            aes(colour="Macarthur Giardia", 
                label=paste("EE: ",mac_ba_g_effect)), 
            size=3, hjust=1, show.legend=FALSE)+
  geom_text(x=366-2,
            y=nr_qmra_Batrend_c[366]*1.15, 
            aes(colour="North Richmond Cryptosporidium", 
                label=paste("EE: ",nr_ba_c_effect)), 
            size=3, hjust=1, show.legend=FALSE)+
  geom_text(x=366-2,
            y=nr_qmra_Batrend_g[366]*1.15, 
            aes(colour="North Richmond Giardia", 
                label=paste("EE: ",nr_ba_g_effect)), 
            size=3, hjust=1, show.legend=FALSE)+
  
  #Set up logarithmic y-axis
  scale_y_continuous(expand=c(0,0), 
                     breaks=trans_breaks("log10", function(x) 10^x), 
                     labels=trans_format("log10", math_format(10^.x)))+
  coord_trans(y="log10")+
  
  #Set up x-axis
  scale_x_continuous(expand=c(0,0))+
  expand_limits(x=c(0,494), y=c(10^-5.5,10^-8.5))+
  labs(x="Time (week)", 
       y="Annual disease burden\n(disability-adjusted life years 
       person\U207B\U00B9 year\U207B\u00B9)")+
  
  #Labels for intervention and tolerable risk lines
  annotate("text", x=370, y=10^-5.6, label="Intervention", hjust=0)+
  annotate("text", x=10, y=10^-5.9, label="Reference level of risk", hjust=0)+

  #Apply aesthetics
  scale_color_viridis_d(labels=c(
    expression(paste("Macarthur ", italic("Cryptosporidium"), " risk")), 
    expression(paste("Macarthur ", italic("Giardia"), " risk")), 
    expression(paste("North Richmond ", italic("Cryptosporidium"), " risk")), 
    expression(paste("North Richmond ", italic("Giardia"), " risk"))))+
  scale_fill_viridis_d()+
  theme_bw()+
  guides(colour=guide_legend(reverse=TRUE, nrow=2, 
                               override.aes=list(linetype=c("solid", 
                                                                "solid",
                                                                "dashed", 
                                                                "dashed"))))+
  theme(legend.title=element_blank(), 
        legend.text.align=0, 
        legend.position="bottom")

figqmra2 #Print the figure