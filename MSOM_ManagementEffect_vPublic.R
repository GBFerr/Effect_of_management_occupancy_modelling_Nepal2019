#########################################################################################

### Multi-species occupancy models for Nepal - data from 2019  ###

#### PART 1 - EFFECT OF MANAGEMENT ####
# model psi ~management+DistRiver, p~Mount+ForCov
# estimates effect of management on occupancy
# estimates spp richness in each management: overall, threat status, and functional group

##########################################################################################
## Loading data
setwd("###/country_files/nepal/our_papers/management_effect_2019CTsurvey")
# loading environment with 3D array with presence-absence matrices for 40 spp
load("###/country_files/nepal/our_papers/management_effect_2019CTsurvey/40NativeSpp_matrix_MSOM_ENV.RData")
# X is the 3D array 
# dimensions of array are: 148 CT sites, 6 5-day survey occasions, and 40 species (25 real spp + 15 all-0 augmented spp)
# spp names can be checked in the list that generated array
names(matrix_native_5d)
# just checking array - records for spp 2 (Chital)
str(X[,,2]) 
X <- type.convert(X) # from numeric to integer

# loading covariates for analysis
covs <- read.csv("scaledCovs_148site_Nepal2019.csv", stringsAsFactors = F) # numerical covs have been scaled
Management <- covs$Management  
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
tail(OBZ)

# distance from rivers
DistRiver <- covs$DistRiver

#type of camera mmount; either tree (T) or pole (P)
Mount <- covs$Mount 
# need to trasnform in number for analysis
# one category must be 0 - this will be intercept for detection
Mount[Mount =="T"] <- 1
Mount[Mount =="P"] <- 0
Mount <- as.numeric(Mount)
unique(Mount)
#forest cover in a 50m buffer around each CT site
ForCov <- covs$propForest50
range(ForCov)
# additional info for the model
n <- dim(X)[3]-15 # number of observed spp (total n of species minus the number of augmented species)
nzeroes <- dim(X)[3]-n # number of augmented (all 0) species in the data set
J <- as.integer(dim(X)[1]) # number of CT sites

# number of survey occasions per site (each 5-day survey occasion)
Ktemp <- X[,,1]
Ktemp[Ktemp==0] <- 1  # 1 indicates CT was working, therefore a valid survey occasion
Ktemp[is.na(Ktemp)] <- 0 # 0 indicates CT was NOT working, therefore not a valid survey occasion and represented as NA in the species presence-abs matrix
Knamed <- rowSums(Ktemp) # sum of rows (i.e. of 1s) indicate number of valid survey occasions at each site
K <- unname(Knamed) # number of surveys at each CT site
Knamed-K # fine, all 0
K <- as.integer(K)
unique(K)
length(K) # must match number of CT sites - ok

# Model used in the paper was originally named 'model 2', which has management and river as covs for occupancy; and forest cover and camera mount as covs for detection
#### model2 psi ~ manag+River p ~ mount+50mForCov, with spprich ####
# Bundle and summarize data to be used in the model (loaded in the previous step)
str(sp.data <- list(n = n, nzeroes = nzeroes, J = J, K = K, X = X,
                    BZ=BZ, NP=NP, OBZ=OBZ, DistRiver=DistRiver,
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
       alphaRiver = rnorm(n = n+nzeroes),
       betaPmount = rnorm(n = n+nzeroes), betaTmount = rnorm(n = n+nzeroes), 
       betaForCov = rnorm(n = n+nzeroes))
}

# changing WD before running to save results in another folder
setwd("###/Nepal2019_occuAnalysis/FinalModels")
# model
sink("m2_managRiver_Feb2021.txt")
cat("
    model {
    # model adapted from Ferreira et al 2020 Biological Conservation
    # Prior distribution for community-level params - hyperpriors
    
    omega ~ dunif(0,1)   # 'hyper-psi'
    mu.BZ  ~ dnorm(0, 0.001)
    mu.NP ~ dnorm(0, 0.001)
    mu.OBZ ~ dnorm(0, 0.001)
    mu.alphaRiver ~ dnorm(0, 0.001)
    mu.betaPmount ~ dnorm(0, 0.001)
    mu.betaTmount ~ dnorm(0, 0.001)
    mu.betaForCov ~ dnorm(0, 0.001)
    
    # tau and sd based on Kery & Royle 2016; 
    tau.BZ <- pow(sd.BZ,-2) 
    tau.NP <- pow(sd.NP,-2) 
    tau.OBZ <- pow(sd.OBZ,-2)
    tau.alphaRiver <- pow(sd.alphaRiver,-2)
    tau.betaPmount <- pow(sd.betaPmount,-2) 
    tau.betaTmount <- pow(sd.betaTmount,-2) 
    tau.betaForCov <- pow(sd.betaForCov,-2) 
    
    sd.BZ ~ dunif(0,2) 
    sd.NP ~ dunif(0,2) 
    sd.OBZ ~ dunif(0,2) 
    sd.alphaRiver ~ dunif(0,2) 
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
    alphaRiver[i] ~ dnorm(mu.alphaRiver, tau.alphaRiver)
    betaPmount[i] ~ dnorm(mu.betaPmount, tau.betaPmount)
    betaTmount[i] ~ dnorm(mu.betaTmount, tau.betaTmount)
    betaForCov[i] ~ dnorm(mu.betaForCov, tau.betaForCov)
    
    # Ecological model for true occurrence (process model) of sp i at site j
    # loop to define Z-matrix ('true' matrix of 1-0), obs: spp loop above still open
    for (j in 1:J) {
    logit(psi[j,i]) <- psiBZ[i]*BZ[j] + psiNP[i]*NP[j] + psiOBZ[i]*OBZ[j] + alphaRiver[i]*DistRiver[j]
    
    mu.psi[j,i] <- psi[j,i] * w[i]
    Z[j,i] ~ dbern(mu.psi[j,i])
    

    # Observation model for replicated detection/nondetection observations; observed 1-0 matrix
    # detetection of species i at site j for survey occasion k; obs spp and site loops still open
    for (k in 1:K[j]) {
    logit(p[j,k,i]) <- betaPmount[i] + betaTmount[i] * Mount[j] + 
                       betaForCov[i] * ForCov[j]
    
    mu.p[j,k,i] <- p[j,k,i] * Z[j,i]
    X[j,k,i] ~ dbern(mu.p[j,k,i])
    
    
    } #  temporal replicate loop
    } #  site loop
    } # close species loop
    
    
    ###### Derived quantities ########
    
    ## porportion of sites occupied (finite sample - fs)
    for (q in 1:(n+nzeroes)){
    psi.fs[q] <- sum(Z[,q])/148       # proportion of sites occupied overall and per manag type
    psi.fs.BZ[q] <- sum(Z[1:50,q])/50  
    psi.fs.NP[q] <- sum(Z[51:100,q])/50  
    psi.fs.OBZ[q] <- sum(Z[101:148,q])/48
    }
    
    ## Site species richness, spp features from NepalSpecies_SizeThreatGuild_info.csv
    for (s in 1:J){
    Nsite[s] <- sum(Z[s,])
    Nsite.LgAnim[s] <- sum(Z[s,c(14, 21, 24)])  # Number of large animalivores in each site 
    Nsite.LgHerbiv[s] <- sum(Z[s,c(2, 4, 8, 17, 18, 20, 23, 25)]) # Number of large herbiv in each site 
    Nsite.SmAnim[s] <- sum(Z[s,c(3, 9, 10, 11, 12, 13, 16, 22)])  # Number of small animaliv in each site 
    Nsite.SmHerbiv[s] <- sum(Z[s,c(1, 5, 6, 7, 15, 19)])  # Number of small herbivores in each site 
    Nsite.reg.threat[s] <- sum(Z[s,c(1,2,3,4,8,9,14,17,18,20,21,23,24)])  # Number regionally threatened spp in each site 
    Nsite.reg.Nthreat[s] <- sum(Z[s,c(5,6,7,10,11,12,13,15,16,19,22,25)])
    Nsite.glob.threat[s] <- sum(Z[s,c(4,5,8,18,20,21,23,24)]) # Number of globally threatened spp in each site 
    Nsite.glob.Nthreat[s] <- sum(Z[s,c(1,2,3,6,7,9,10,11,12,13,14,15,16,17,19,22,25)])
    }
    ## Mean site spp richness for each type of management
    # All species
    Nsite.BZ <- Nsite[1:50]    # sprich per site for BZ; no need to monitor
    Nsite.NP <-Nsite[51:100]          # sprich per site for NP; no need to monitor 
    Nsite.OBZ <-Nsite[101:148]  # sprich per site for OBZ; no need to monitor 
   
    mean.Nsite.BZ <- mean(Nsite.BZ) # mean SPrich per site for BZ, param to be monitored
    mean.Nsite.NP <- mean(Nsite.NP)   # mean SPrich per site for NP, param to be monitored
    mean.Nsite.OBZ <- mean(Nsite.OBZ)   # mean SPrich per site for OBZ, param to be monitored
    
    Ntotal <- sum(w[])                 # Total metacommunity size
    
    }
    ",fill = TRUE)
sink()


params1 <- c("omega", "mu.BZ", "mu.NP", "mu.OBZ", "mu.alphaRiver",
             "mu.betaPmount", "mu.betaTmount", "mu.betaForCov", 
             "psiBZ", "psiNP", "psiOBZ", "alphaRiver",
             "betaPmount", "betaTmount", "betaForCov",
             "psi.fs", "psi.fs.BZ", "psi.fs.NP", "psi.fs.OBZ",
             "Nsite","Nsite.LgAnim","Nsite.SmAnim","Nsite.LgHerbiv","Nsite.SmHerbiv",
             "Nsite.reg.threat", "Nsite.reg.Nthreat", "Nsite.glob.threat", "Nsite.glob.Nthreat",
             "mean.Nsite.BZ", "mean.Nsite.NP", "mean.Nsite.OBZ",
             "Ntotal")



# MCMC settings
ni <- 150000   ;   nt <- 10   ;   nb <- 50000   ;   nc <- 3

# Run JAGS, check convergence and summarize posteriors
library(jagsUI)
#Tiny diff when compared with other versions of this model
# only increased 10k iterations (in ni and nb) and included mu.alphaRiver in the parameterrs to be monitored
outm2_managRiver_Feb<- jags(sp.data, sp.inits, params1, "m2_managRiver_Feb2021.txt", 
                                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# ran for 7hs20m

save.image("###/Nepal2019_occuAnalysis/FinalModels/m2_managRiver_vFeb21_ENV.RData")
# to avoid running model again can load inputs and output from model 2 using the comand below

m2.summary <- as.data.frame(outm2_managRiver_Feb$summary) # summary of model output
write.csv(m2.summary, file="summarymodel2_managRiver_vFeb2021.csv") #saving m2 summary as csv
head(m2.summary,8)
# checking convergence
m2.summary[m2.summary$Rhat>1.1,] 
mean(m2.summary$Rhat)
# end of the model

###### EXPLORING AND PLOTTING MODEL OUTPUT ######

#### SPECIES OCCUPANCY RESULTS  ####
setwd("C:/Users/Guilherme/Dropbox/!BHP_storage/Nepal2019_occuAnalysis/FinalModels/Feb2021") # crashed, starting new session
spp.info.2 <- read.csv("NepalSpecies_SizeThreatGuild_info_v2.csv") # spp traits

library(cowplot)
library(ggthemes)
library(viridis)
library(cowplot)
library(gridExtra)
library(grid)

#### calculating difference in predicted occupancy (probability scale) between management regimes ####
str(tmp <- (outm2_managRiver_Feb$sims.list)) # getting  posterior values of  MCMC chains for all parameters monitored
nsamp <- length (tmp[[1]]) # number of MCM sampling
diff.NP.BZ.Prob<- array(NA,dim=c(nsamp, 25)) # it's 2 dimensions only because it's not predicting over different values of a cov (that would be the 3rd dim)
#25 in the array above represents the 25 spp assessed
str(diff.NP.BZ.Prob)
#getting difference between occu in NP and occu in BZ in each MCMC sample for 25 spp
for(i in 1:nsamp){
  for(s in 1:25){
    diff.NP.BZ.Prob[i,s] <- plogis(tmp$psiNP[i,s]) - plogis(tmp$psiBZ[i,s]) # diff NP-BZ
  }}
head(diff.NP.BZ.Prob)    
dim(diff.NP.BZ.Prob) #checking
mean.NP.BZ.Prob <- apply(diff.NP.BZ.Prob, 2, mean)  # getting the posterior mean of the difference; 2 is dimensions of the array representing spp
# getting 95% credible interval of differences in occu
cri.NP.BZ.Prob <- apply(diff.NP.BZ.Prob, 2, function (x) quantile (x, probs = c(0.025, 0.975)))
NP.BZ.effect.Prob <- rbind(mean.NP.BZ.Prob, cri.NP.BZ.Prob)

# rearranging DF
NP.BZ.effect.Prob <- as.data.frame(t(NP.BZ.effect.Prob))
NP.BZ.effect.Prob$pair <- rep("NP-BZ",25)
names(NP.BZ.effect.Prob) <- c("mean", "lower", "upper", "pair")

# NP - OBZ
diff.NP.OBZ.Prob<- array(NA,dim=c(nsamp, 25)) 
#25 in the array above represents the 25 spp assessed
str(diff.NP.OBZ.Prob)
#getting difference between occu in NP and occu in OBZ (coefProb)
for(i in 1:nsamp){
  for(s in 1:25){
    diff.NP.OBZ.Prob[i,s] <- plogis(tmp$psiNP[i,s]) - plogis(tmp$psiOBZ[i,s]) # diff NP-OBZ
  }}
head(diff.NP.OBZ.Prob)    

mean.NP.OBZ.Prob <- apply(diff.NP.OBZ.Prob, 2, mean)  # getting the posterior mean of the differences in occu, 2 is dimensions of the array representing spp
# getting 95% credible interval of differences in occu
cri.NP.OBZ.Prob <- apply(diff.NP.OBZ.Prob, 2, function (x) quantile (x, probs = c(0.025, 0.975)))
NP.OBZ.effect.Prob <- rbind(mean.NP.OBZ.Prob, cri.NP.OBZ.Prob)
# rearranging
NP.OBZ.effect.Prob <- as.data.frame(t(NP.OBZ.effect.Prob))
NP.OBZ.effect.Prob$pair <- rep("NP-OBZ",25)
names(NP.OBZ.effect.Prob) <- c("mean", "lower", "upper", "pair")

# BZ - OBZ
diff.BZ.OBZ.Prob<- array(NA,dim=c(nsamp, 25))  #25 spp assessed
str(diff.BZ.OBZ.Prob)
#getting difference between occu in BZ and occu in OBZ (coefProb)
for(i in 1:nsamp){
  for(s in 1:25){
    diff.BZ.OBZ.Prob[i,s] <- plogis(tmp$psiBZ[i,s]) - plogis(tmp$psiOBZ[i,s]) # diff BZ-OBZ
  }}
head(diff.BZ.OBZ.Prob)    
dim(diff.BZ.OBZ.Prob)   

mean.BZ.OBZ.Prob <- apply(diff.BZ.OBZ.Prob, 2, mean)  # getting the posterior mean of the differences
# getting 95% credible interval of differences
cri.BZ.OBZ.Prob <- apply(diff.BZ.OBZ.Prob, 2, function (x) quantile (x, probs = c(0.025, 0.975)))
BZ.OBZ.effect.Prob <- rbind(mean.BZ.OBZ.Prob, cri.BZ.OBZ.Prob)
# rearranging
BZ.OBZ.effect.Prob <- as.data.frame(t(BZ.OBZ.effect.Prob))
BZ.OBZ.effect.Prob$pair <- rep("BZ-OBZ",25)
names(BZ.OBZ.effect.Prob) <- c("mean", "lower", "upper", "pair")

# joining results for the 3 pair-wise comparisos of occupancy estimate (NP-BZ; NP-OBZ; BZ-OBZ)
manag.effect.occu.Prob <- rbind(NP.BZ.effect.Prob, NP.OBZ.effect.Prob, BZ.OBZ.effect.Prob)
# col indicating CRI values that overlap 0
manag.effect.occu.Prob$ovelap0 <- manag.effect.occu.Prob$lower * manag.effect.occu.Prob$upper
manag.effect.occu.Prob[manag.effect.occu.Prob$ovelap0>0, 5] <- 0  # positive value indicates no overlap
manag.effect.occu.Prob[manag.effect.occu.Prob$ovelap0<0, 5] <- 1  # negative output of multiplication means 95% CRI overlap 0 
range(manag.effect.occu.Prob$ovelap0) #checking
#adding spp names
manag.effect.occu.Prob$spp_names <- spp.info.ltd$Tag  # using 'spp.info' which has the 25 spp names repeated 3 times in the correct order
sum(manag.effect.occu.Prob$ovelap0) # number of cases in which difference in occu overlaps 0
head(manag.effect.occu.Prob)
write.csv(manag.effect.occu.Prob,"managEffect_25spp_occupancy_Probability.csv")


#### plotting  difference in probabilities (relative to the species-specific estimate in OBZ)
manag.effect.Prob <- cbind(manag.effect.occu.Prob, pred.occu.m2)
head(manag.effect.Prob)
manag.effect.Prob <- manag.effect.Prob[manag.effect.Prob$n_recs>4,-c(7:10,13)] # keeping only spp with >4 recs and eliminating unecessary cols
# eliminating NP-BZ comparison, this way OBZ will be the baseline in the graph
manag.effect.Prob <- manag.effect.Prob[manag.effect.Prob$pair!="NP-BZ",]

# replacing Indian grey mongoose for grey mongoose only to make name smaller
manag.effect.Prob$Tag <- sub("indian grey", "grey", manag.effect.Prob$Tag)

manag.effect.Prob <- manag.effect.Prob %>% #renaming 'size' col to avoid confusion with ggplot
  rename("Size_categ" = "size")

OBZ.est <- as.matrix(plogis(m2.summary[89:113,1]))
nrecs <- as.matrix(pred.occu.m2$n_recs[1:25,1])
OBZ.est <- as.data.frame(cbind(OBZ.est,nrecs))
names(OBZ.est) <- c("meanOBZ", "nrecs")
OBZ.est <- OBZ.est[OBZ.est$nrecs>4,] # keeping only spp with >4 records
OBZ.est <- rbind(OBZ.est,OBZ.est)

manag.effect.Prob2 <-cbind(manag.effect.Prob,OBZ.est)


manag.effect.Prob2$Tag <- factor(manag.effect.Prob2$Tag, 
                                 levels =c("jungle cat","jackal", "grey mongoose","nilgai","hare","leopard",
                                           "small indian civet","wild boar", "swamp deer", "four-horned antelope",
                                           "masked palm civet","hog deer", "one-horned rhino", "sloth bear", "macaque",
                                           "porcupine","tiger",
                                           "barking deer","sambar","grey langur","chital"))

# Fig. 2 - plotting effect of management on spp occupancy (OBZ as baseline)
library(ggplot2)
ggplot(manag.effect.Prob2, aes(x=Tag, y=mean, ymax = upper, ymin = lower, colour= pair),  #aes(x=reorder(Tag, mean)
       fatten = 3, lwd=1)+#, shape= size)) + 
  geom_pointrange(position=position_dodge2(width=.6))+
  geom_hline(yintercept=0, linetype="dashed") +
  scale_y_continuous(breaks = seq(-1,1,by=0.25))+
  theme_bw()+ # white background (as opposed to the default grey)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_flip()+
  theme(axis.title.y = element_blank(),
        axis.text.x  = element_text(size=10),
        axis.text.y  = element_text(size=10),
        axis.title.x = element_text(size=11, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.position=c(0.8,0.15),
        legend.title = element_blank()) + 
  scale_color_manual(values=c("#FC8D62", "#66C2A5"))+
  ylab("Difference in occupancy probability")

save.image("###/Nepal2019_occuAnalysis/FinalModels/Feb2021/m2mangRiver_Feb2021_ENV.RData")



#### SPECIES RICHNESS RESULTS ####
spprich.m2 <- as.data.frame(m2.summary[449:1780,c(1:3,7,8)]) # getting only the spp rich parameters
tail(spprich.m2)
head(spprich.m2)
spprich.m2$rich_categ <- c(rep("A) overall", 148),rep("G) large animaliv.", 148),rep("F) small animaliv.", 148),
                           rep("E) large herbiv.", 148),rep("D) small herbiv.", 148),rep("C) threatened", 148),
                           rep("B) non-threatened", 148),rep("threatened global", 148),
                           rep("non-threatened global", 148))


library(dplyr)
spprich.m2 <- spprich.m2 %>% 
  rename("lower" = "2.5%",
         "upper" = "97.5%")
head(spprich.m2)
spprich.m2$manag <- c(rep("BZ", 50), rep("NP",50), rep("OBZ",48),
                      rep("BZ", 50), rep("NP",50), rep("OBZ",48),
                      rep("BZ", 50), rep("NP",50), rep("OBZ",48),
                      rep("BZ", 50), rep("NP",50), rep("OBZ",48),
                      rep("BZ", 50), rep("NP",50), rep("OBZ",48),
                      rep("BZ", 50), rep("NP",50), rep("OBZ",48),
                      rep("BZ", 50), rep("NP",50), rep("OBZ",48),
                      rep("BZ", 50), rep("NP",50), rep("OBZ",48),
                      rep("BZ", 50), rep("NP",50), rep("OBZ",48))
spprich.m2$manag <- factor(spprich.m2$manag, levels = c("NP", "BZ", "OBZ"))

# now plotting

#### calculating difference in spp richness between management regimes for m2 ####
str(tmp <- (outm2_managRiver_Feb$sims.list))
tmp <- tmp[20:28] # keeping only spprich elements 
str(tmp) #checking
tmp[1] #checking

# calculating the difference in spp rich between management regimes
diff.spprich<- array(NA,dim=c(3,3,9))
#rows 1:50 = BZ; 51:100 = NP; 101-148 = OBZ
for(n in 1:9){
  NP.BZ.temp <- (rowMeans(as.data.frame(tmp[[n]][,51:100])))-(rowMeans(as.data.frame(tmp[[n]][,1:50]))) 
  diff.spprich[1,1,n] <- mean(NP.BZ.temp)
  diff.spprich[1,2:3,n] <- quantile(NP.BZ.temp, probs= c(0.025,0.975))
  NP.OBZ.temp <- (rowMeans(as.data.frame(tmp[[n]][,51:100])))-(rowMeans(as.data.frame(tmp[[n]][,101:148]))) 
  diff.spprich[2,1,n] <- mean(NP.OBZ.temp)
  diff.spprich[2,2:3,n] <- quantile(NP.OBZ.temp, probs= c(0.025,0.975))
  BZ.OBZ.temp <- (rowMeans(as.data.frame(tmp[[n]][,1:50])))-(rowMeans(as.data.frame(tmp[[n]][,101:148]))) 
  diff.spprich[3,1,n] <- mean(BZ.OBZ.temp)
  diff.spprich[3,2:3,n] <- quantile(BZ.OBZ.temp, probs= c(0.025,0.975))
}

diff.spprich[,,9] # checking result

# arranging resulting DF with diff in spp rich per manag #
colnames(diff.spprich) <- c("mean", "LCI", "UCI")
rownames(diff.spprich) <- c("NP-BZ", "NP-OBZ", "BZ-OBZ")
subset.names <- names(tmp)
# naming the 3rd dimension of the array w/ subset of spp rich (e.g.: small, medium , large spp rich)
dimnames(diff.spprich)[[3]] <- subset.names 
str(diff.spprich)

# rearranging the resulting 3d array
library(geomorph)
mtx.diff.spprich <- two.d.array(diff.spprich) # transforming 3d arry in 2d matrix
#not quite what I wanted it still needs some manual handling...
temp1 <- mtx.diff.spprich[,1:3]
temp2 <- mtx.diff.spprich[,4:6]
temp3 <- mtx.diff.spprich[,7:9]
new.mtx.spprich <- as.data.frame(rbind(temp1,temp2,temp3))
new.mtx.spprich$subset <- c(rep("NP-BZ",9), rep("NP-OBZ",9), rep("BZ-OBZ",9)) #col identifying pairwise comparisons

mtx.sorted <- new.mtx.spprich[order(rownames(new.mtx.spprich)),] #rearranging row order
colnames(mtx.sorted) <- c("mean", "LCI", "UCI", "manag.pair") # more appropriate col names
mtx.sorted$subset <- rownames(mtx.sorted)
row.names(mtx.sorted) <- 1:27

# col indicating CRI values that overlap 0
mtx.sorted$overlap0 <- mtx.sorted$LCI * mtx.sorted$UCI
mtx.sorted[mtx.sorted$overlap0>0, 6] <- 0  # positive value indicates CRI dos not overlap 0
mtx.sorted[mtx.sorted$overlap0<0, 6] <- 1  # negative output of multiplication means 95% CRI overlap 0 


#reordering cols in DF
col_order <- c("subset", "manag.pair", "mean",
               "LCI", "UCI", "overlap0")
mean.diff.spprich <- mtx.sorted[,col_order]
mean.diff.spprich

#final reordering of rows
mean.diff.spprich$index <- c(1,1,1,9,9,9,8,8,8,5,5,5,4,4,4,6,6,6,7,7,7,3,3,3,2,2,2)
mean.diff.spprich <- mean.diff.spprich[order(mean.diff.spprich$index),]
mean.diff.spprich$subset <- sub(".1","",mean.diff.spprich$subset)
mean.diff.spprich$subset <- sub(".2","",mean.diff.spprich$subset)
write.csv(mean.diff.spprich,"mean_diff_spprich_manag_comparison_m2.csv")

# mean spprich per management 
spprich.manag<- array(NA,dim=c(3,3,9))
for(n in 1:9){
  spprich.manag[1,1,n] <- mean(tmp[[n]][,51:100])
  temp.NP <- rowMeans(as.data.frame(tmp[[n]][,51:100]))
  spprich.manag[1,2:3,n] <- quantile(temp.NP, probs= c(0.025,0.975))
  spprich.manag[2,1,n] <- mean(tmp[[n]][,1:50])
  temp.BZ <- rowMeans(as.data.frame(tmp[[n]][,1:50]))
  spprich.manag[2,2:3,n] <- quantile(temp.BZ, probs= c(0.025,0.975))
  spprich.manag[3,1,n] <- mean(tmp[[n]][,101:148])
  temp.OBZ <- rowMeans(as.data.frame(tmp[[n]][,101:148]))
  spprich.manag[3,2:3,n] <- quantile(temp.OBZ, probs= c(0.025,0.975))
}
spprich.manag[,,1]
spprich.manag.2 <- two.d.array(spprich.manag)

row.names(spprich.manag.2) <- names(tmp)
colnames(spprich.manag.2) <- c("NP", "LCI", "UCI", "BZ", "LCI", "UCI", "OBZ", "LCI", "UCI")

np <- spprich.manag.2[,1:3]
bz <- spprich.manag.2[,4:6]
obz <- spprich.manag.2[,7:9]
spprich.manag.3 <- as.data.frame(rbind(np,bz,obz))
spprich.manag.3$manag <- c(rep("NP",9), rep("BZ",9), rep("OBZ",9))
spprich.manag.3 <- spprich.manag.3[order(rownames(spprich.manag.3)),]
spprich.manag.3$subset <- row.names(spprich.manag.3)
row.names(spprich.manag.3) <- 1:27
#final reordering of rows
spprich.manag.3$index <- c(1,1,1,9,9,9,8,8,8,5,5,5,4,4,4,6,6,6,7,7,7,3,3,3,2,2,2)
spprich.manag.3 <- spprich.manag.3[order(spprich.manag.3$index),]
colnames(spprich.manag.3) <- c("mean", "LCI", "UCI", "manag", "subset", "index")
col_order.2 <- c("subset", "manag", "mean", "LCI", "UCI", "index")
spprich.manag.3 <- spprich.manag.3[,col_order.2]
spprich.manag.3$subset <- sub(".1","",spprich.manag.3$subset)
spprich.manag.3$subset <- sub(".2","",spprich.manag.3$subset)
spprich.manag.3
write.csv(spprich.manag.3,"mean_spprich_permanag_m2.csv")

##### Fig.3: boxplots of spp rich in each management ####
# overall, threat level, and functional groups

overall <- spprich.m2[1:148,]
fgs <- spprich.m2[149:740,]
threat<- spprich.m2[741:1332,]

library(ggplot2)
library(ggsignif)
library(scales)

ggplot(overall, aes(x = manag, y = mean)) +
  labs(x = "Management", y = "Site species richness")+
  facet_wrap(~rich_categ, ncol=1, scales="free_y")+
  scale_y_continuous(breaks = seq(2,13,by=3))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = manag), alpha = 0.85,
              position = position_jitter(width = 0.3))+
  #scale_y_continuous(labels = accuracy = 1) +
  scale_colour_brewer(palette = "Set2")+
  theme_bw()+
  expand_limits(y=c(2,12))+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  #ggtitle("Overall")+ 
  geom_signif(y_position = c(10.3, 11.3), xmin = c(2, 1), xmax = c(3, 3),
  annotation = c("*", "*"), tip_length = 0.01)

#### THREAT
threat.2 <- threat[1:296,]
#threat.2$rich_categ <- sub(" Nepal","", threat.2$rich_categ)
threat.boxplot <- ggplot(threat.2, aes(x = manag, y = mean)) +
  labs(x = "Management", y = "Site species richness")+
  facet_wrap(~rich_categ, ncol=2, scales="free_y")+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = manag), alpha = 0.85,
              position = position_jitter(width = 0.3))+
  scale_colour_brewer(palette = "Set2")+
  theme_bw()+
  #expand_limits(y=c(2,12))+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        #axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
        axis.title.y = element_blank())

# adding annotations for pair-wise comparisons
annotation_df <- data.frame(
  rich_categ = c("B) non-threatened","B) non-threatened","C) threatened","C) threatened","C) threatened"),
  start = c(1,2,1,1,2),
  end = c(2,3,2,3,3),
  y = c(7.3,8,7.5,8.3,6),
  label = c("*", "* ", " *", " * ", " *  "))

threat.boxplot+
  geom_signif(data = annotation_df,aes(
              y_position = y, xmin = start, xmax = end,
              annotations = label, tip_length = 0.01), manual = TRUE)+
  ylim(NA,8.5)

#### functional group
fgs$rich_categ <- factor(fgs$rich_categ , levels = c("D) small herbiv.","E) large herbiv.",
                                                     "F) small animaliv.", "G) large animaliv."))

# changing DF to increase range of y-axis for individual facets
unique(fgs$rich_categ)
fgs$ymin <- rep(c(0,0.5,0,0), each=148)
fgs$ymax <- rep(c(2.55,5,5.5,4.7), each=148)

FG.boxplot <- ggplot(fgs, aes(x = manag, y = mean)) +
  labs(x = "Management", y = "Site species richness")+
  facet_wrap(~rich_categ, ncol=4, scales="free")+
  geom_blank(aes(y = ymin)) + geom_blank(aes(y = ymax))+ # change range of y
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = manag), alpha = 0.85,
              position = position_jitter(width = 0.3))+
  scale_colour_brewer(palette = "Set2")+
  theme_bw()+
  #expand_limits(y=c(2,12))+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

annotationFG_df <- data.frame(
  rich_categ = c("F) small animaliv.","D) small herbiv.","D) small herbiv.","G) large animaliv.",
                 "G) large animaliv.", "E) large herbiv.", "E) large herbiv."),
  start = c(1,1,2,1,1,1,2),
  end = c(3,3,3,2,3,3,3),
  y = c(4,4.4,3.7,2.2,2.4,5,4.4),
  label = c("*", "* ", " *", " * ", " *  ", "  *  ", "  *   "))

#annotationFG_df$rich_categ <- factor(annotationFG_df$rich_categ , levels = c("small herbiv.","large herbiv.", 
#                                                                             "small animaliv.","large animaliv."))

FG.boxplot+
  geom_signif(data = annotationFG_df,aes(
    y_position = y, xmin = start, xmax = end,
    annotations = label, tip_length = 0.01), manual = TRUE)

save.image("C:/Users/Guilherme/Dropbox/!BHP_storage/Nepal2019_occuAnalysis/FinalModels/Feb2021/m2mangRiver_Feb2021_ENV.RData")


