#########################################################################################

### Multi-species occupancy models for Nepal - data from 2019  ###

#### PART 2 - EFFECT OF AGRICULTURE  ####
#(interaction between management and agriculture)
# model psi ~manage*agric, p ~mount + ForCov_50m
# estimates effect of proportion of agriculture on occupancy in the three management regimes assessed 

##########################################################################################
setwd("###/country_files/nepal/our_papers/management_effect_2019CTsurvey")
# loading spp matrices and covs
load("###/country_files/nepal/our_papers/management_effect_2019CTsurvey/40NativeSppMatrix&Covs_MSOM_ENV.RData")
# X is the 3D array 
# dimensions of array are: 148 CT sites, 6 5-day survey occasions, and 40 species (25 real spp + 15 all-0 augmented spp)
# preparing data
X <- type.convert(X) # from numeric to integer

# just checking covs; numerical covs are scaled
head(covs) 
Management <- as.character(covs$Management)  
# transforming Management in indicator variables
BZ <- Management
BZ[BZ=="BZ"] <- 1
BZ[BZ=="NP" | BZ=="OBZ"] <- 0
BZ <- as.numeric(BZ)
NP <- Management
NP[NP=="NP"] <- 1
NP[NP=="BZ" | NP=="OBZ"] <- 0
NP <- as.numeric(NP)
OBZ <- Management
OBZ[OBZ=="OBZ"] <- 1
OBZ[OBZ=="BZ" | OBZ=="NP"] <- 0
OBZ <- as.numeric(OBZ)
sum(OBZ)

# proportion of agriculture in 500-m buffer
propAgric500 <- covs$propAgric500
range(propAgric500)
#type of camera mmount; either tree (T) or pole (P)
covs$Mount <- sub(pattern=" ", replacement="", covs$Mount) #just excluding 'space' from two sites
Mount <- covs$Mount 
unique(Mount)
# need to trasnform in number for analysis
# one category must be 0 - this will be intercept for detection
Mount[Mount =="T"] <- 1
Mount[Mount =="P"] <- 0
Mount <- as.numeric(Mount)

#forest cover in a 50m buffer around each CT site
ForCov <- covs$propForest50
range(ForCov)
# addition info for the model
n <- dim(X)[3]-15 # number of observed spp (total n of species minus the number of augmented species)
nzeroes <- dim(X)[3]-n # number of augmented (all 0) species in the data set
J <- as.integer(dim(X)[1]) # number of CT sites
# number of survey occasions per site (each 5-day survey occasion)
Ktemp <- X[,,1]
Ktemp[Ktemp==0] <- 1  # 1 indicates CT was working, therefore a valid survey occasion
Ktemp[is.na(Ktemp)] <- 0 # 0 indicates CT was NOT working, therefore not a valid survey occasion and represented as NA in the species presence-abs matrix
Knamed <- rowSums(Ktemp) # sum of rows (i.e. of 1s) indicate numberr of valid survey occasions at each site
K <- unname(Knamed) # number of surveys at each CT site
Knamed-K # fine, all 0
K <- as.integer(K)
unique(K)
length(K) # must match number of CT sites - ok

# Bundle and summarize data
str(sp.data <- list(n = n, nzeroes = nzeroes, J = J, K = K, X = X,
                    BZ=BZ, NP=NP, OBZ=OBZ,
                    PropAgric=propAgric500,
                    Mount = Mount, ForCov = ForCov))

# Defining initial values for the MCMC
wst <- rep(1, n+nzeroes)             # Simply set everybody at occurring
zst <- array(1, dim = c(J, n+nzeroes)) # ditto
range(wst)
range(zst)

sp.inits <- function() {
  omegaGuess = runif(1, n/(n+nzeroes), 1)
  psi.meanGuess = runif(1, .25,1)
  list(omega=omegaGuess, Z = zst, w = wst,
       psiBZ = rnorm(n = n+nzeroes), 
       psiNP = rnorm(n = n+nzeroes), psiOBZ = rnorm(n = n+nzeroes), 
       alphaBZAgric = rnorm(n = n+nzeroes), alphaNPAgric = rnorm(n = n+nzeroes), 
       alphaOBZAgric = rnorm(n = n+nzeroes), 
       betaPmount = rnorm(n = n+nzeroes), betaTmount = rnorm(n = n+nzeroes), 
       betaForCov = rnorm(n = n+nzeroes))
    }


# model
sink("model_managAgric.txt")
cat("
    model {
    # model adapted from Ferreira et al 2020 Biological Conservation
    # Prior distribution for community-level params - hyperpriors
    
    omega ~ dunif(0,1)   # 'hyper-psi'
    mu.BZ  ~ dnorm(0, 0.001)
    mu.NP ~ dnorm(0, 0.001)
    mu.OBZ ~ dnorm(0, 0.001)
    mu.alphaBZAgric  ~ dnorm(0, 0.001)
    mu.alphaNPAgric ~ dnorm(0, 0.001)
    mu.alphaOBZAgric ~ dnorm(0, 0.001)
    mu.betaPmount ~ dnorm(0, 0.001)
    mu.betaTmount ~ dnorm(0, 0.001)
    mu.betaForCov ~ dnorm(0, 0.001)
    
    # tau and sd based on Kery & Royle 2016; 
    tau.BZ <- pow(sd.BZ,-2) 
    tau.NP <- pow(sd.NP,-2) 
    tau.OBZ <- pow(sd.OBZ,-2)
    tau.alphaBZAgric  <- pow(sd.alphaBZAgric,-2)
    tau.alphaNPAgric  <- pow(sd.alphaNPAgric,-2)
    tau.alphaOBZAgric <- pow(sd.alphaOBZAgric,-2)
    tau.betaPmount <- pow(sd.betaPmount,-2) 
    tau.betaTmount <- pow(sd.betaTmount,-2) 
    tau.betaForCov <- pow(sd.betaForCov,-2) 
    
    sd.BZ ~ dunif(0,2) 
    sd.NP ~ dunif(0,2) 
    sd.OBZ ~ dunif(0,2) 
    sd.alphaBZAgric ~ dunif(0,2)
    sd.alphaNPAgric ~ dunif(0,2) 
    sd.alphaOBZAgric ~ dunif(0,2) 
    sd.betaPmount ~ dunif(0,2) 
    sd.betaTmount ~ dunif(0,2) 
    sd.betaForCov ~ dunif(0,2) 
    
    # specify species-level priors for species i (out of 40) influenced by comm-level priors
    # this is exactly the same in Kery & Royle and in Zipkin
    for (i in 1:(n+nzeroes)) {
    
    w[i] ~ dbern(omega)   # Superpopulation process: Ntotal species sampled out of all spp available
    psiBZ[i] ~ dnorm(mu.BZ, tau.BZ)
    psiNP[i] ~ dnorm(mu.NP, tau.NP)
    psiOBZ[i] ~ dnorm(mu.OBZ, tau.OBZ)
    alphaBZAgric[i] ~ dnorm(mu.alphaBZAgric, tau.alphaBZAgric)
    alphaNPAgric[i] ~ dnorm(mu.alphaNPAgric, tau.alphaNPAgric)
    alphaOBZAgric[i] ~ dnorm(mu.alphaOBZAgric, tau.alphaOBZAgric)
    betaPmount[i] ~ dnorm(mu.betaPmount, tau.betaPmount)
    betaTmount[i] ~ dnorm(mu.betaTmount, tau.betaTmount)
    betaForCov[i] ~ dnorm(mu.betaForCov, tau.betaForCov)
    
    # Ecological model for true occurrence (process model) of sp i at site j
    # loop to define Z-matrix ('true' matrix of 1-0), obs: spp loop above still open
    for (j in 1:J) {
    logit(psi[j,i]) <- psiBZ[i]*BZ[j] + psiNP[i]*NP[j] + psiOBZ[i]*OBZ[j] + 
                        alphaBZAgric[i]*BZ[j]*PropAgric[j] + alphaNPAgric[i]*NP[j]*PropAgric[j] +
                        alphaOBZAgric[i]*OBZ[j]*PropAgric[j]
    
    mu.psi[j,i] <- psi[j,i] * w[i]
    Z[j,i] ~ dbern(mu.psi[j,i])
    

    # Observation model for replicated detection/nondetection observations; observed 1-0 matrix
    # detetection of species i at site j for survey occasion k; obs spp and site loops still open
    for (k in 1:K[j]) {
    logit(p[j,k,i]) <- betaPmount[i] + betaTmount[i] * Mount[j] + 
                       betaForCov[i] * ForCov[j]
    
    mu.p[j,k,i] <- p[j,k,i] * Z[j,i]
    X[j,k,i] ~ dbern(mu.p[j,k,i])
    
    Xnew[j,k,i] ~ dbern(mu.p[j,k,i])  # replicate data
    
    # model assessment
    # Observed dataset
    chi2.actual[j,k,i] <- pow(X[j,k,i] - mu.p[j,k,i], 2)/ (mu.p[j,k,i] + 0.0001)  # Add small value to denominator to prevent division by zero
    
    # Expected dataset  
    chi2.sim[j,k,i] <- pow(Xnew[j,k,i] - mu.p[j,k,i], 2)/ (mu.p[j,k,i] + 0.0001)  # Add small value to denominator to prevent division by zero
    
    } #  temporal replicate loop
    
    chi2.actual.sum[j,i] <- sum(chi2.actual[j,1:K[j],i])
    chi2.sim.sum[j,i] <- sum(chi2.sim[j,1:K[j],i])  

    } #  site loop
    } # close species loop
    
    
    # Calculate the discrepancy measure (Pearsons chi-squared) for each spp, which is then defined as the mean(p.fit > p.fitnew)
    # note that this loop does not include the all-zero spp
    for(g in 1:n){
    fit.sp.actual[g] <- sum(chi2.actual.sum[,g])      # Species-specific fit statistic for actual dataset
    fit.sp.sim[g] <- sum(chi2.sim.sum[,g])            # Species-specific fit statistic for simulated dataset
    
    c.hat.sp[g] <- fit.sp.actual[g]/fit.sp.sim[g]
    bpv.sp[g] <- step(fit.sp.sim[g] - fit.sp.actual[g])
    }
    
    # Calculate overall Pearson's chi-squared discrepency measure, defined post-hoc as: mean(p.fit>p.fit.new)
    fit.actual <- sum(chi2.actual.sum[1:J, 1:n])
    fit.sim <- sum(chi2.sim.sum[1:J, 1:n])
    
    c.hat <- fit.actual/fit.sim
    
    bpv <- step(fit.sim - fit.actual)    
    
    ###### Derived quantities ########
    ## Bayesian p-values - model assessment
    for(o in 1:n){
    spp.Pvalue[o] <- mean(fit.sp.actual[o] > fit.sp.sim[o])
    }
    overall.Pvalue <- mean(fit.actual > fit.sim)
    
    ## porportion of sites occupied (finite sample - fs)
    for (q in 1:(n+nzeroes)){
    psi.fs[q] <- sum(Z[,q])/148       # proportion of sites occupied overall and per manag type
    psi.fs.BZ[q] <- sum(Z[1:50,q])/50  
    psi.fs.NP[q] <- sum(Z[51:100,q])/50  
    psi.fs.OBZ[q] <- sum(Z[101:148,q])/48
    }
    
    ## Site species richness 
    for (s in 1:J){
    Nsite[s] <- sum(Z[s,])          # Number of occurring species at each site 
    }
    
    ## Mean site spp richness for each type of management
    # All species
    Nsite.BZ <- Nsite[1:50]    # sprich per site for BZ; no need to monitor
    Nsite.NP <-Nsite[51:100]          # sprich per site for NP; no need to monitor 
    Nsite.OBZ <-Nsite[101:148]  # sprich per site for OBZ; no need to monitor 
   
    mean.Nsite.BZ <- mean(Nsite.BZ) # mean SPrich per site for BZ, param to be monitored
    mean.Nsite.NP <- mean(Nsite.NP)   # mean SPrich per site for NP, param to be monitored
    mean.Nsite.OBZ <- mean(Nsite.OBZ)   # mean SPrich per site for OBZ, param to be monitored
    
    diff.Nsite.NP.BZ <- mean.Nsite.NP - mean.Nsite.BZ # difference in spp rich NP vs BZ
    diff.Nsite.NP.OBZ <- mean.Nsite.NP - mean.Nsite.OBZ
    diff.Nsite.BZ.OBZ <- mean.Nsite.BZ - mean.Nsite.OBZ
    
    Ntotal <- sum(w[])                 # Total metacommunity size
    
    
    }
    ",fill = TRUE)
sink()

params1 <- c("omega", "mu.BZ", "mu.NP", "mu.OBZ", "mu.alphaBZAgric", "mu.alphaNPAgric", "mu.alphaOBZAgric",
             "mu.betaPmount", "mu.betaTmount", "mu.betaForCov", 
             "psiBZ", "psiNP", "psiOBZ", "alphaBZAgric", "alphaNPAgric", "alphaOBZAgric",
             "betaPmount", "betaTmount", "betaForCov",
             "psi.fs", "psi.fs.BZ", "psi.fs.NP", "psi.fs.OBZ",
             "Nsite",
             "mean.Nsite.BZ", "mean.Nsite.NP", "mean.Nsite.OBZ",
             "diff.Nsite.NP.BZ","diff.Nsite.NP.OBZ","diff.Nsite.BZ.OBZ", 
             "Ntotal","c.hat.sp", "c.hat", "bpv.sp", "bpv", 
             "spp.Pvalue", "overall.Pvalue")

# MCMC settings
ni <- 150000   ;   nt <- 10   ;   nb <- 50000   ;   nc <- 3

# Run JAGS, check convergence and summarize posteriors
library(jagsUI)
outmodel_managAgric<- jags(sp.data, sp.inits, params1, "model_managAgric.txt", 
                                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

save.image("###/country_files/nepal/our_papers/management_effect_2019CTsurvey/output_managAgric_ENV.RData")
summary_managAgric <- as.data.frame(outmodel_managAgric$summary)
write.csv(summary_managAgric, file="summarymodel_managAgric.csv")

summary_managAgric[1:15,c(1:3,7,8)]

# model ends here

#### EXPLORING AND PLOTTING MODEL RESULTS ####
# getting model coeffs
agric.coeffs.comm.spp<- summary_managAgric[131:250,] # spp coeffs
tail(agric.coeffs.comm.spp,30)
head(agric.coeffs.comm.spp)
temp.comm.agric <-summary_managAgric[5:7,] # community coeffs

agric.coeffs.comm.spp <- rbind(temp.comm.agric,agric.coeffs.comm.spp) # joining comm & spp coeffs
agric.coeffs.comm.spp$parameter <- rownames(agric.coeffs.comm.spp)

library(tidyverse)
agric.coeffs.comm.spp <-separate(agric.coeffs.comm.spp, parameter, c("parameter", "spIndex"), sep="\\[")
agric.coeffs.comm.spp$spIndex <- sub(pattern="\\]", replacement="", agric.coeffs.comm.spp$spIndex)
agric.coeffs.comm.spp[1:3,13] <- c(97,98,99) # adding random Index number for comm parameters
agric.coeffs.comm.spp$spIndex <- as.numeric(agric.coeffs.comm.spp$spIndex)
agric.coeffs.comm.spp <- agric.coeffs.comm.spp[,c(1,3,7,12:13)]
agric.coeffs.comm.spp <- agric.coeffs.comm.spp %>% 
  rename(
    'LCI'='2.5%',
    'UCI'='97.5%'
  )

agric.coeffs.comm.spp$management <- c(rep("BZ",43), rep("NP",40), rep("OBZ",40))

# resetting the first 3 rows of management
agric.coeffs.comm.spp[1:3,6]=c("BZ","NP", "OBZ")

#eliminating all params from NP as only want to show effect of agric in BZ and OBZ

agric.coeffs.comm.spp <- agric.coeffs.comm.spp[agric.coeffs.comm.spp$management!="NP",]
str(agric.coeffs.comm.spp)
# using sppIndex to keep spp with 5+recs and community params
agric.coeffs.comm.spp <- filter(agric.coeffs.comm.spp, spIndex==97| spIndex==99| spIndex==1|spIndex==2|spIndex==6|spIndex==7|spIndex==8
                                |spIndex==10|spIndex==11|spIndex==12|spIndex==15|spIndex==16|spIndex==17
                                |spIndex==20|spIndex==22|spIndex==25)

agric.coeffs.comm.spp$species <- c("community", "community", "barking deer", "chital", "grey langur", "hare",
                                   "hog deer", "grey mongoose", "jackal", "jungle cat", "macaque", "masked palm civet",
                                   "nilgai", "sambar", "small indian civet", "wild boar", 
                                   "barking deer", "chital", "grey langur", "hare",
                                   "hog deer", "grey mongoose", "jackal", "jungle cat", "macaque", "masked palm civet",
                                   "nilgai", "sambar", "small indian civet", "wild boar")
agric.coeffs.comm.spp$row <- 1:30

# need to eliminate params for spp with <5 recs in a management regime
# b deer , sambar, g langur in OBZ [19, 22, 28]
# h deer and palm civet in BZ [7,12]
agric.coeffs.comm.spp.red <- agric.coeffs.comm.spp[-c(7,12,17,19,28),]

agric.coeffs.comm.spp.red$species <- 
  factor(agric.coeffs.comm.spp.red$species, #ordering by type of response
         levels = c("community",
                    "nilgai", "wild boar", "macaque","hare","jackal",
                    "grey mongoose",
                    "chital","small indian civet", "jungle cat",
                    "barking deer", "sambar","grey langur",
                    "hog deer", "masked palm civet"))


#### Fig 4A - plotting model coefficients ####
ggplot(agric.coeffs.comm.spp.red , aes(x=species, #x=reorder(species, sortSpp), 
                                       y=mean, ymax = UCI, ymin = LCI, colour= management),  
       fatten = 4, lwd=2)+#, shape= size)) + 
  geom_pointrange(position=position_dodge2(width=.6))+
  geom_hline(yintercept=0, linetype="dashed") +
  #scale_y_continuous(breaks = seq(-1,1,by=0.25))+
  theme_bw()+ # white background (as opposed to the default grey)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #coord_flip()+
  theme(axis.title.y = element_text(size=11, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x  = element_text(size=10, angle=45, hjust = 1),
        axis.text.y  = element_text(size=10),
        axis.title.x = element_blank(),
        legend.position="bottom",
        #legend.position=c(0.80,0.85),
        legend.title = element_blank()) + 
  scale_color_manual(values=c("#FC8D62","#8DA0CB"))+
  ylab("Model coefficients")

### predicting hyper-parameters management*Agric; community response ####
pred.propAgric <- seq(from = min(propAgric500), to = max(propAgric500),length.out=250) # scaled values of Dist_water to plug into the prediction formula

# now actually predicting
str(tmp <- (outmodel_managAgric$sims.list))
nsamp <- length (tmp[[1]]) # number of MCM samplings
predmanagAgric <- array(NA,dim=c(250, nsamp, 3))

for(i in 1:nsamp){
  predmanagAgric[,i,1] <- plogis(tmp$mu.BZ[i] + tmp$mu.alphaBZAgric[i] * pred.propAgric) # interaction between Agric and manag PA
  predmanagAgric[,i,2] <- plogis(tmp$mu.NP[i] + tmp$mu.alphaNPAgric[i] * pred.propAgric) # 
  predmanagAgric[,i,3] <- plogis(tmp$mu.OBZ[i] + tmp$mu.alphaOBZAgric[i] * pred.propAgric) 
 }
str(predmanagAgric)

predmanagAgric.mean <- apply(predmanagAgric, c(1,3), mean)  # getting the posterior mean of the predicted value 

# getting 95% credible interval of predictions
cri.predmanagAgric <- apply(predmanagAgric, c(1,3), function (x) quantile (x, probs = c(0.025, 0.975)))
lower.predmanagAgric <- cri.predmanagAgric[1,,] # just spliting array into DF
upper.predmanagAgric <- cri.predmanagAgric[2,,] # ditto
full.predmanagAgric <- as.data.frame(cbind(predmanagAgric.mean, lower.predmanagAgric, upper.predmanagAgric))
head(full.predmanagAgric)
names(full.predmanagAgric) <- c("BZ_Agric", 
                       "NP_Agric", 
                       "OBZ_Agric", 
                       "LCI_BZ_Agric", 
                       "LCI_NP_Agric", 
                       "LCI_OBZ_Agric",
                       "UCI_BZ_Agric", 
                       "UCI_NP_Agric", 
                       "UCI_OBZ_Agric")

library(reshape2)
tmp2 <- melt(full.predmanagAgric) # will stack all collumns and create a new col with variable name
tail(tmp2,20)  #checking
tmp2[740:760,]#checking
mean <- tmp2[1:750 ,]
LCI <- tmp2[751:1500,]
UCI <- tmp2[1501:2250,]
tail(UCI)

new.predmanagAgric <- cbind(mean,LCI,UCI)
head(new.predmanagAgric)
names(new.predmanagAgric) <- c("var","mean","var","LCI", "var","UCI")
         
new.predmanagAgric$management <- c(rep("BZ",250), rep("NP",250), rep("OBZ",250))
#just checking
new.predmanagAgric[245:255,] #ok
new.predmanagAgric[495:505,]#ok

# adding cols with real value for variables; will be used in x-axis of graphs
#range propAgric (%): 0-100
real.PropAgric <- as.data.frame(seq(from = 0, to = 100,length.out=250)) 
real.PropAgric <- rbind(real.PropAgric, real.PropAgric, real.PropAgric)
names(real.PropAgric) <- "propAgric"
new.predmanagAgric$propAgric <- unlist(real.PropAgric)
head(new.predmanagAgric)
new.predmanagAgric <- new.predmanagAgric[,c(2,4,6:8)] # just eliminating repeated cols
head(new.predmanagAgric)


#### Fig. 4B - plotting community prediction for BZ and OBZ ####
# eliminating NP as agric is negligible in this mgmnt
plot.managAgric <- new.predmanagAgric[new.predmanagAgric$management != "NP",] 
plot.managAgric$facet <- rep("community",500)

comm.agric.plot <- ggplot(plot.managAgric, aes(x=propAgric, y=mean))+
  facet_wrap(~facet, scales="free_y", nrow=1) +
  theme_bw()+ # white background (as opposed to the default grey)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=management),alpha=0.2)+
  geom_line(aes(colour=management, linetype=management), size=1.5) + 
  scale_fill_manual(values=c("#FC8D62", "#8DA0CB"))+
  scale_colour_manual(values=c("#FC8D62", "#8DA0CB"))+
  scale_linetype_manual(values=c("solid", "longdash"))+
  scale_size_manual(values=c(1, 5))+
  theme(legend.position=c(.85,.85))+
  theme(legend.title=element_blank())+
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text.x  = element_text(size=10),
        axis.text.y  = element_text(size=10)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  ylab("Occupancy probability")+
  xlab("Proportion of agricultural land")

#### getting coeffs and predictions of agric effect for all spp with at least 5 recs in either BZ and OBZ ####
# first finding spp with at least 5 recs in either BZ and OBZ
# splitting the 3D spp matrices into BZ and OBZ
apply(X,c(3),sum, na.rm=T) # sum of records overall
BZ.recs <- X[1:50,,]
OBZ.recs <- X[101:148,,]
BZ.OBZ.recs <- X[c(1:50,101:148),,] # sum of records in BZ and OBZ


BZ.sum <- apply(BZ.recs,c(3),sum, na.rm=T)
OBZ.sum <- apply(OBZ.recs,c(3),sum, na.rm=T)
BZ.OBZ.sum <- cbind(BZ.sum,OBZ.sum)
index <- 1:40
BZ.OBZ.sum <- cbind(BZ.OBZ.sum,index)
# spp w/ 5+ recs in either BZ or OBZ:
#                       BZ.sum OBZ.sum index
# BarkingDeer             55       0     1
# Chital                 100      17     2
# GreyLangur              10       0     6
# Hare                    15      11     7
# HogDeer                  0       6     8
# IndianGreyMongoose       7      10    10
# Jackal                  50      35    11
# JungleCat               20      44    12
# Macaque                 54      19    15
# MaskedPalmCivet          3       5    16
# Nilgai                  17      14    17
# Sambar                  14       0    20
# SmallIndianCivet         9      11    22
# WildBoar                83      35    25


# predicting spp responses
pred.propAgric <- seq(from = min(propAgric500), to = max(propAgric500),length.out=100) # scaled values of prop Agric to plug into the prediction formula

str(tmp <- (outmodel_managAgric$sims.list))
nsamp <- length (tmp[[1]]) # number of MCM samplings
spp.predmanagAgric.5plus <- array(NA,dim=c(100, nsamp, 23)) #3rd dimension is for the spp coeffs, only spp with 5+ recs in agiven mgmnt
str(spp.predmanagAgric.5plus)
for(i in 1:nsamp){
  spp.predmanagAgric.5plus[,i,1] <- plogis(tmp$psiBZ[i,1] + tmp$alphaBZAgric[i,1] * pred.propAgric) # interaction between Agric and BZ for spp 1
  spp.predmanagAgric.5plus[,i,2] <- plogis(tmp$psiBZ[i,2] + tmp$alphaBZAgric[i,2] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,3] <- plogis(tmp$psiBZ[i,6] + tmp$alphaBZAgric[i,6] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,4] <- plogis(tmp$psiBZ[i,7] + tmp$alphaBZAgric[i,7] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,5] <- plogis(tmp$psiBZ[i,10] + tmp$alphaBZAgric[i,10] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,6] <- plogis(tmp$psiBZ[i,11] + tmp$alphaBZAgric[i,11] * pred.propAgric) 
  spp.predmanagAgric.5plus[,i,7] <- plogis(tmp$psiBZ[i,12] + tmp$alphaBZAgric[i,12] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,8] <- plogis(tmp$psiBZ[i,15] + tmp$alphaBZAgric[i,15] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,9] <- plogis(tmp$psiBZ[i,17] + tmp$alphaBZAgric[i,17] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,10] <- plogis(tmp$psiBZ[i,20] + tmp$alphaBZAgric[i,20] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,11] <- plogis(tmp$psiBZ[i,22] + tmp$alphaBZAgric[i,22] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,12] <- plogis(tmp$psiBZ[i,25] + tmp$alphaBZAgric[i,25] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,13] <- plogis(tmp$psiOBZ[i,2] + tmp$alphaOBZAgric[i,2] * pred.propAgric) # interaction between Agric and OBZ for spp 2
  spp.predmanagAgric.5plus[,i,14] <- plogis(tmp$psiOBZ[i,7] + tmp$alphaOBZAgric[i,7] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,15] <- plogis(tmp$psiOBZ[i,8] + tmp$alphaOBZAgric[i,8] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,16] <- plogis(tmp$psiOBZ[i,10] + tmp$alphaOBZAgric[i,10] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,17] <- plogis(tmp$psiOBZ[i,11] + tmp$alphaOBZAgric[i,11] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,18] <- plogis(tmp$psiOBZ[i,12] + tmp$alphaOBZAgric[i,12] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,19] <- plogis(tmp$psiOBZ[i,15] + tmp$alphaOBZAgric[i,15] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,20] <- plogis(tmp$psiOBZ[i,16] + tmp$alphaOBZAgric[i,16] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,21] <- plogis(tmp$psiOBZ[i,17] + tmp$alphaOBZAgric[i,17] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,22] <- plogis(tmp$psiOBZ[i,22] + tmp$alphaOBZAgric[i,22] * pred.propAgric)
  spp.predmanagAgric.5plus[,i,23] <- plogis(tmp$psiOBZ[i,25] + tmp$alphaOBZAgric[i,25] * pred.propAgric)
}


spp.predmanagAgric.5plus.mean <- apply(spp.predmanagAgric.5plus, c(1,3), mean)  # getting the posterior mean of the predicted value 

# getting 95% credible interval of predictions
cri.spp.predmanagAgric.5plus <- apply(spp.predmanagAgric.5plus, c(1,3), function (x) quantile (x, probs = c(0.025, 0.975)))
lower.spp.predmanagAgric.5plus <- cri.spp.predmanagAgric.5plus[1,,] # just spliting array into DF
upper.spp.predmanagAgric.5plus <- cri.spp.predmanagAgric.5plus[2,,] # ditto
full.spp.predmanagAgric.5plus <- as.data.frame(cbind(spp.predmanagAgric.5plus.mean, lower.spp.predmanagAgric.5plus, upper.spp.predmanagAgric.5plus))
head(full.spp.predmanagAgric.5plus)
names(full.spp.predmanagAgric.5plus) <- c("BZ_barking deer", "BZ_chital", "BZ_grey langur", "BZ_hare","BZ_grey mongoose",
                                                 "BZ_jackal","BZ_jungle cat", "BZ_macaque","BZ_nilgai",  "BZ_sambar", "BZ_small indian civet",
                                                 "BZ_wild boar",
                                                 "OBZ_chital", "OBZ_hare", "OBZ_hog deer", "OBZ_grey mongoose", "OBZ_jackal","OBZ_jungle cat", 
                                                 "OBZ_macaque", "OBZ_masked palm civet","OBZ_nilgai", "OBZ_small indian civet","OBZ_wild boar",
                                                 "LCI_BZ_barking deer", "LCI_BZ_chital", "LCI_BZ_grey langur", "LCI_BZ_hare","LCI_BZ_grey mongoose",
                                                 "LCI_BZ_jackal","LCI_BZ_jungle cat", "LCI_BZ_macaque","LCI_BZ_nilgai", "LCI_BZ_sambar", "LCI_BZ_small indian civet",
                                                 "LCI_BZ_wild boar",
                                                 "LCI_OBZ_chital", "LCI_OBZ_hare", "LCI_OBZ_hog deer", "LCI_OBZ_grey mongoose", "LCI_OBZ_jackal","LCI_OBZ_jungle cat", 
                                                 "LCI_OBZ_macaque", "LCI_OBZ_masked palm civet","LCI_OBZ_nilgai", "LCI_OBZ_small indian civet","LCI_OBZ_wild boar",
                                                 "UCI_BZ_barking deer", "UCI_BZ_chital", "UCI_BZ_grey langur", "UCI_BZ_hare","UCI_BZ_grey mongoose",
                                                 "UCI_BZ_jackal","UCI_BZ_jungle cat", "UCI_BZ_macaque","UCI_BZ_nilgai", "UCI_BZ_sambar", "UCI_BZ_small indian civet",
                                                 "UCI_BZ_wild boar",
                                                 "UCI_OBZ_chital", "UCI_OBZ_hare", "UCI_OBZ_hog deer", "UCI_OBZ_grey mongoose", "UCI_OBZ_jackal","UCI_OBZ_jungle cat", 
                                                 "UCI_OBZ_macaque", "UCI_OBZ_masked palm civet","UCI_OBZ_nilgai", "UCI_OBZ_small indian civet","UCI_OBZ_wild boar")


str(full.spp.predmanagAgric.5plus)

library(reshape2)
tmp3 <- melt(full.spp.predmanagAgric.5plus) # will stack all columns and create a new col with variable name
tail(tmp3,101)  #checking
tmp3[995:1005,]#checking
mean <- tmp3[1:2300 ,]
LCI <- tmp3[2301:4600,]
UCI <- tmp3[4601:6900,]
tail(UCI)
head(mean)

new.spp.5plus.pred <- cbind(mean,LCI,UCI)
head(new.spp.5plus.pred)
names(new.spp.5plus.pred) <- c("var","mean","var","LCI", "var","UCI")
new.spp.5plus.pred <- new.spp.5plus.pred[,-c(3,5)] # just eliminating unecessary cols

library(dplyr)
library(tidyr)
new.spp.5plus.pred <- new.spp.5plus.pred %>%
  separate(var, c("management", "species"), "_")

new.spp.5plus.pred[245:255,] #ok
new.spp.5plus.pred[495:505,]#ok
new.spp.5plus.pred[1195:1205,]#ok


# adding cols with real value for variables; will be used in x-axis of graphs
#range propAgric (%): 0-100
real.PropAgric <- as.data.frame(seq(from = 0, to = 100,length.out=100)) 
real.PropAgric <- rbind(real.PropAgric, real.PropAgric, real.PropAgric,
                        real.PropAgric, real.PropAgric, real.PropAgric,
                        real.PropAgric, real.PropAgric, real.PropAgric,
                        real.PropAgric, real.PropAgric, real.PropAgric, real.PropAgric,
                        real.PropAgric, real.PropAgric, real.PropAgric,
                        real.PropAgric, real.PropAgric, real.PropAgric,
                        real.PropAgric, real.PropAgric, real.PropAgric,
                        real.PropAgric) # stacking it 23 times, once per coeff predicted
names(real.PropAgric) <- "propAgric"
new.spp.5plus.pred$propAgric <- unlist(real.PropAgric)
head(new.spp.5plus.pred)

new.spp.5plus.pred$management <- factor(new.spp.5plus.pred$management, levels = c("BZ", "OBZ"))

str(new.spp.5plus.pred)
# need to create blank species for layout purposes and reduce spp name to fit facet
new.spp.5plus.pred.blank <- new.spp.5plus.pred
temp.blank <- data.frame(management= rep("BZ",100),
                         species= rep("blank",100),
                         mean= rep(0.1,100),
                         LCI=rep(0.1,100),
                         UCI=rep(0.1,100),
                         propAgric=pred.propAgric)
new.spp.5plus.pred.blank <- rbind(new.spp.5plus.pred.blank,temp.blank)
new.spp.5plus.pred.blank$spp <- new.spp.5plus.pred.blank$species
new.spp.5plus.pred.blank$spp <- sub("small indian civet", "s. indian civet",new.spp.5plus.pred.blank$spp)
new.spp.5plus.pred.blank$spp <- sub("masked palm civet", "m. palm civet",new.spp.5plus.pred.blank$spp)
new.spp.5plus.pred.blank$spp <- sub("grey mongoose", "g. mongoose",new.spp.5plus.pred.blank$spp)
new.spp.5plus.pred.blank$spp <- factor(new.spp.5plus.pred.blank$spp, #ordering by threat
                                              levels = c("barking deer", "sambar","chital", "nilgai","hog deer", "blank",
                                                         "wild boar","s. indian civet", "g. mongoose","grey langur","m. palm civet",
                                                         "macaque","hare","jackal", "jungle cat"))

#### Fig 4B - plotting spp-level predictions #### 

ggplot(new.spp.5plus.pred.blank, aes(x=propAgric, y=mean))+
  facet_wrap(~spp, scales="free_y", nrow=5) +
  theme_bw()+ # white background (as opposed to the default grey)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=management),alpha=0.2)+
  geom_line(aes(colour=management, linetype=management), size=1) + 
  scale_fill_manual(values=c("#FC8D62","#8DA0CB"))+
  scale_colour_manual(values=c("#FC8D62","#8DA0CB"))+
  scale_linetype_manual(values=c("solid", "longdash"))+
  scale_size_manual(values=c(1, 5))+
  #theme(legend.position="top")+
  theme(legend.position=c(.85,.75))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size = 12))+
  theme(#strip.text = element_text(size=10,lineheight=5.0),
    strip.background = element_rect(fill="gray90", colour="black"))+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x  = element_text(size=10),
        axis.text.y  = element_text(size=10)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  ylab("Occupancy probability")+
  xlab("Proportion of agricultural land")




#plotting                     
ggplot(agric.coeffs.comm.spp.red , aes(x=species, #x=reorder(species, sortSpp), 
                                   y=mean, ymax = UCI, ymin = LCI, colour= management),  
       fatten = 4, lwd=2)+#, shape= size)) + 
  geom_pointrange(position=position_dodge2(width=.6))+
  geom_hline(yintercept=0, linetype="dashed") +
  #scale_y_continuous(breaks = seq(-1,1,by=0.25))+
  theme_bw()+ # white background (as opposed to the default grey)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #coord_flip()+
  theme(axis.title.y = element_text(size=11, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x  = element_text(size=10, angle=45, hjust = 1),
        axis.text.y  = element_text(size=10),
        axis.title.x = element_blank(),
        legend.position="bottom",
        #legend.position=c(0.80,0.85),
        legend.title = element_blank()) + 
  scale_color_manual(values=c("#FC8D62","#8DA0CB"))+
  #ylab("Effect of agriculture \n(model coefficient)")
  ylab("Model coefficients")

save.image("###/Nepal2019_occuAnalysis/FinalModels/Jan2022/model_managAgric_vJan2022_ENV.RData")
