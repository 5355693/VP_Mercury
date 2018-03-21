rm(list = ls())
data.df <- read.csv(url("https://raw.githubusercontent.com/5355693/VP_Mercury/master/hg_data.csv"))
str(data.df)
library(ggplot2)
library(plyr)
library(lubridate)
data.df$Sample_Date2 <- as.POSIXct(data.df$Sample_Date, format = "%m/%d/%y") #create a true date variable.
data.df$Year <- year(data.df$Sample_Date2) #Create a year variable to allow us to separate the 2016 observation.
data.df$individual <- seq(from = 1,to = nrow(data.df), by = 1) #create a unique identifier for each individual.

summary(data.df)
library(R2jags)
data_subset2 <- subset(data.df, Life_Stage != "Adult")
data_subset2 <- subset(data_subset2, (!is.na(data_subset2$Amphib_MeHg)))

#Simplest model: Methylmercury in amphibians ~ Species + Water MeHg + random effect of pool
spp <- data_subset2$Spp
hg <- log(data_subset2$Amphib_MeHg)
pool <- data_subset2$Pool
waterhg <- log(data_subset2$Water_MeHg) #could consider scaling this to get smaller Betas
n <- nrow(data_subset2)
n.groups <- length(levels(data_subset2$Spp))
n.pools <- length(levels(data_subset2$Pool))
jags.params <- c("alpha","alpha0","beta.waterhg", "sigma", "mu","bpvalue") #add mu to monitor predictions!
jags.inits <- function(){
  list(sigma=rlnorm(1))
}
#Model 1: Amphib_MeHg ~ Species + waterhg + 1|Pool
m1 <- function () {
  for(i in 1:n){
    hg[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha0[pool[i]] + alpha[spp[i]] + beta.waterhg*waterhg[i]
  }
  for(i in 1:n.groups){
    alpha[i] ~ dnorm(0,0.001) 
  }
  for(i in 1:n.pools){#random pool effect
    alpha0[i] ~ dnorm(mu.alpha, tau.alpha)
  }
  for(i in 1:n){
    residual[i] <- hg[i]-mu[i]
    predicted[i] <- mu[i]
    sq[i] <- pow(residual[i],2)
    y.new[i] ~ dnorm(mu[i], tau)
    sq.new[i] <- pow(y.new[i] - predicted[i],2)
  }
  fit<- sum(sq[])
  fit.new <- sum(sq.new[])
  test <- step(fit.new - fit)
  bpvalue <- mean(test)
  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)
  beta.waterhg ~ dnorm (0,0.001)
  mu.alpha ~ dnorm(0,0.001)
  tau.alpha <- 1/(sigma.alpha*sigma.alpha)
  sigma.alpha ~ dunif(0,100)
}

jagsfitm1 <- jags.parallel(data = c("spp","pool","hg","waterhg","n","n.groups", "n.pools"), inits = jags.inits, jags.params,
                     n.iter = 500000, n.thin = 100,model.file = m1)
print(jagsfitm1,digits = 5)
traceplot(jagsfitm1)
dev.off()
jagsfitm1$BUGSoutput$mean$bpvalue

## Chains do not look to be mixing well and convergence is poor. It may be necessary to use
## pool as a fixed effect.
spp <- data_subset2$Spp
hg <- log(data_subset2$Amphib_MeHg)
pool <- data_subset2$Pool
waterhg <- log(data_subset2$Water_MeHg) #could consider scaling this to get smaller Betas
n <- nrow(data_subset2)
n.groups <- length(levels(data_subset2$Spp))
n.pools <- length(levels(data_subset2$Pool))
jags.params <- c("alpha","beta","beta.waterhg", "sigma","bpvalue") #add mu to monitor predictions!
jags.inits <- function(){
  list(sigma=rlnorm(1))
}

m1.FE <- function () {
  for(i in 1:n){
    hg[i] ~ dnorm(mu[i], tau)
    mu[i] <-alpha[spp[i]] + beta[pool[i]] + beta.waterhg*waterhg[i]
  }
  for(i in 1:n.groups){
    alpha[i] ~ dnorm(0, 0.001) 
  }
  for (i in 1:n.pools){
    beta[i] ~ dnorm(0, 0.001)
  }
  for(i in 1:n){
    residual[i] <- hg[i]-mu[i]
    predicted[i] <- mu[i]
    sq[i] <- pow(residual[i],2)
    y.new[i] ~ dnorm(mu[i], tau)
    sq.new[i] <- pow(y.new[i] - predicted[i],2)
  }
  fit<- sum(sq[])
  fit.new <- sum(sq.new[])
  test <- step(fit.new - fit)
  bpvalue <- mean(test)
  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)
  beta.waterhg ~ dnorm (0, 0.001)
}

jagsfitm1.FE <- jags.parallel(data = c("spp","hg","waterhg","n","n.groups", "n.pools", "pool"), inits = jags.inits, jags.params,
                  n.iter = 50000, model.file = m1.FE)
print(jagsfitm1.FE,digits = 5)
traceplot(jagsfitm1.FE)
dev.off()
jagsfitm1.FE$BUGSoutput$mean$bpvalue

##Compare expected to observed values:
plot(jagsfitm1$BUGSoutput$median$mu,log(data_subset2$Amphib_MeHg))
abline(0,1)

##Confirming that the JAGS and ML models are the same:
(m1_ML <- lm(log(Amphib_MeHg) ~ Spp-1 + log(Water_MeHg) + Pool, data = data_subset2))
summary(m1_ML)

plot(jagsfitm1$BUGSoutput$median$mu,predict(m1_ML, newdata = data_subset2, type = "response"),
     xlab = "JAGS predictions", 
     ylab = "ML predictions")
abline(0,1)

library(lme4)
REmod <- lmer(log(Amphib_MeHg) ~ (1|Pool) + Spp-1 + log(Water_MeHg), 
              data = data_subset2)
plot(REmod)
summary(REmod)
plot(predict(REmod), log(data_subset2$Amphib_MeHg))
abline(0,1)
FEmod <- glm(log(Amphib_MeHg) ~ Spp-1 + log(Water_MeHg), data = data_subset2)
summary(FEmod)
plot(predict(FEmod), log(data_subset2$Amphib_MeHg))
abline(0,1)
plot(predict(FEmod), predict(REmod))
abline(0,1)

#More complex model: Methylmercury in amphibians ~ Species + Water MeHg + Life_Stage
spp <- data_subset2$Spp
hg <- log(data_subset2$Amphib_MeHg)
stage <- data_subset2$Life_Stage
stage <- droplevels(stage)
waterhg <- log(data_subset2$Water_MeHg) 
n <- nrow(data_subset2)
n.groups <- length(levels(data_subset2$Spp))
n.stages <- 3
jags.params <- c("alpha", "beta.waterhg", "waterhg","sigma",
                 "fit","fit.new","bpvalue","residual","predicted", "beta.stage") #add mu to monitor predictions!
jags.inits <- function(){
  list(sigma=rlnorm(1))
}
#Model 2: Amphib_MeHg ~ Spp*Water_MeHg + Life_Stage
m2 <- function () {
  for(i in 1:n){
    hg[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha[spp[i]] + beta.waterhg*waterhg[i] + beta.stage[stage[i]]
    #esthg[i] <- exp(mu[i])
  }
  for(i in 1:n.groups){
    alpha[i] ~ dnorm(0, 0.001) 
  }
  for(i in 1:n.stages){
    beta.stage[i] ~ dnorm(0, 0.001)
  }
  for(i in 1:n){
    residual[i] <- hg[i]-mu[i]
    predicted[i] <- mu[i]
    sq[i] <- pow(residual[i],2)
    y.new[i] ~ dnorm(mu[i], tau)
    sq.new[i] <- pow(y.new[i] - predicted[i],2)
  }
  fit<- sum(sq[])
  fit.new <- sum(sq.new[])
  test <- step(fit.new - fit)
  bpvalue <- mean(test)
  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)
  beta.waterhg ~ dnorm (0,0.001)
}

jagsfitm2 <- jags(data = c("spp","hg","waterhg","n","n.groups","stage", "n.stages"), inits = jags.inits, jags.params,
                     n.iter = 50000, model.file = m2)
print(jagsfitm2,digits = 5)
traceplot(jagsfitm2)
dev.off()

##Bayesian P-value suggests a good fit (0.525):
jagsfitm2$BUGSoutput$mean$bpvalue

##Compare expected to observed values:
plot(jagsfitm2$BUGSoutput$median$mu,log(data_subset2$Amphib_MeHg))
abline(0,1)

##Confirming that the JAGS and ML models are the same:
(m2_ML <- lm(log(Amphib_MeHg) ~ Spp + log(Water_MeHg) + Life_Stage, data = data_subset2))
plot(jagsfitm2$BUGSoutput$median$mu,predict(m2_ML, newdata = data_subset2, type = "response"),
     xlab = "JAGS predictions", 
     ylab = "ML predictions")
abline(0,1)

#The most complex model:  Methylmercury in amphibians ~ Species + Water MeHg + Life_Stage + Habitat
spp <- data_subset2$Spp
hg <- log(data_subset2$Amphib_MeHg)
stage <- data_subset2$Life_Stage
stage <- droplevels(stage)
waterhg <- log(data_subset2$Water_MeHg) #could consider scaling this to get smaller Betas
habitat <- data_subset2$Habitat
n <- nrow(data_subset2)
n.groups <- length(levels(data_subset2$Spp))
n.habitat <- 2
n.stages <- 3
jags.params <- c("alpha","beta.waterhg", "beta.stage", "beta.habitat","sigma","esthg",
                 "fit","fit.new","bpvalue","residual","predicted", "mu") #add mu to monitor predictions!
jags.inits <- function(){
  list(sigma=rlnorm(1))
}
#Model 3: Amphib_MeHg ~ Spp*Water_MeHg + Life_Stage + Habitat
m3 <- function () {
  for(i in 1:n){
    hg[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha[spp[i]] + beta.waterhg*waterhg[i] + beta.stage[stage[i]] + beta.habitat[habitat[i]]
    esthg[i] <- exp(mu[i])
  }
  for(i in 1:n.groups){
    alpha[i] ~ dnorm(0, 0.001) 
  }
  for(i in 1:n.stages){
    beta.stage[i] ~ dnorm(0, 0.001)
  }
  for(i in 1:n.habitat){
    beta.habitat[i] ~ dnorm(0, 0.001)
  }
  for(i in 1:n){
    residual[i] <- hg[i]-mu[i]
    predicted[i] <- mu[i]
    sq[i] <- pow(residual[i],2)
    y.new[i] ~ dnorm(mu[i], tau)
    sq.new[i] <- pow(y.new[i] - predicted[i],2)
  }
  fit<- sum(sq[])
  fit.new <- sum(sq.new[])
  test <- step(fit.new - fit)
  bpvalue <- mean(test)
  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)
  beta.waterhg ~ dnorm(0,0.001)
}

jagsfitm3 <- jags(data = c("spp","hg","waterhg","habitat","n","n.groups","n.habitat","stage", "n.stages"), inits = jags.inits, jags.params,
                     n.iter = 50000, model.file = m3)

print(jagsfitm3,digits = 5) 

traceplot(jagsfitm3)
dev.off()

##Bayesian P-value suggests a good fit (0.509):
jagsfitm3$BUGSoutput$mean$bpvalue

##Compare expected to observed values:
plot(jagsfitm3$BUGSoutput$median$mu,log(data_subset2$Amphib_MeHg))
abline(0,1)

##Confirming that the JAGS and ML models are the same:
(m3_ML <- lm(log(Amphib_MeHg) ~ Spp + log(Water_MeHg) + Life_Stage + habitat, data = data_subset2))
plot(jagsfitm3$BUGSoutput$median$mu,predict(m3_ML, newdata = data_subset2, type = "response"),
     xlab = "JAGS predictions", 
     ylab = "ML predictions")
abline(0,1)

#Best model so far is Spp*WaterHg + Stage:
cat(paste(c("DIC for Spp + WaterHg:"), jagsfitm1$BUGSoutput$DIC))
cat(paste(c("DIC for Spp + WaterHg + Stage:"), jagsfitm2$BUGSoutput$DIC))
cat(paste(c("DIC for Spp + WaterHg + Stage + Habitat:"),jagsfitm3$BUGSoutput$DIC))

#Overall fit is good. Can visualize this graphically: not systematically under- or over-estimating error.
plot(jagsfitm2$BUGSoutput$sims.list$fit,jagsfitm2$BUGSoutput$sims.list$fit.new, main = 
       "Graphical posterior predictive check", las = 1, xlab = "SSQ for actual data set",
     ylab = "SSQ for ideal data sets")
abline(0,1)

##Generating predictions for plotting and inference:
spp <- data_subset2$Spp
hg <- log(data_subset2$Amphib_MeHg)
stage <- data_subset2$Life_Stage
stage <- droplevels(stage)
waterhg <- log(data_subset2$Water_MeHg)
waterhgseq <- seq(min(waterhg), max(waterhg),length.out = 111)
n <- nrow(data_subset2)
n.groups <- length(levels(data_subset2$Spp))
n.stages <- 3
jags.params <- c("alpha", "beta.stage","beta.waterhg", "waterhg","waterhgseq","sigma","predWOFRegg",
                 "predWOFRearly","predWOFRlate","predSPSAegg","predSPSAearly","predSPSAlate","fit",
                 "fit.new","bpvalue","residual","predicted") #added mu to monitor predictions!
jags.inits <- function(){
  list(sigma=rlnorm(1))
}
#Model Amphib_MeHg ~ Spp + Water_MeHg + Life_Stage
predmodel <- function () {
  for(i in 1:n){
    hg[i] ~ dnorm(mu[i], tau)
    #esthg[i] <- exp(mu[i])
    mu[i] <- alpha[spp[i]] + beta.waterhg*waterhg[i] + beta.stage[stage[i]]
    predWOFRegg[i] <- alpha[spp[5]] + beta.waterhg*waterhgseq[i] + beta.stage[stage[1]]
    predWOFRearly[i] <-alpha[spp[5]] + beta.waterhg*waterhgseq[i] + beta.stage[stage[5]]
    predWOFRlate[i] <- alpha[spp[5]] + beta.waterhg*waterhgseq[i] + beta.stage[stage[13]]
    predSPSAegg[i] <- alpha[spp[1]] + beta.waterhg*waterhgseq[i] + beta.stage[stage[1]]
    predSPSAearly[i] <- alpha[spp[1]] + beta.waterhg*waterhgseq[i] + beta.stage[stage[5]]
    predSPSAlate[i] <- alpha[spp[1]] + beta.waterhg*waterhgseq[i] + beta.stage[stage[13]]
  }
  for(i in 1:n.groups){
    alpha[i] ~ dnorm(0, 0.001) 
  }
  for(i in 1:n.stages){
    beta.stage[i] ~ dnorm(0, 0.001)
  }
  for(i in 1:n){
    residual[i] <- hg[i]-mu[i]
    predicted[i] <- mu[i]
    sq[i] <- pow(residual[i],2)
    y.new[i] ~ dnorm(mu[i], tau)
    sq.new[i] <- pow(y.new[i] - predicted[i],2)
  }
  fit<- sum(sq[])
  fit.new <- sum(sq.new[])
  test <- step(fit.new - fit)
  bpvalue <- mean(test)
  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)
  beta.waterhg ~ dnorm(0,0.001)
}

jagsfitpredmodel <- jags(data = c("spp","hg","waterhg","n","n.groups","stage", "n.stages","waterhgseq"), inits = jags.inits, jags.params,
                         n.iter = 50000, model.file = predmodel)
#These should be the same:
print(jagsfitpredmodel,digits = 5)
print(jagsfitm2, digits = 5)

pred.matrix <- matrix(nrow = 111,ncol=19,dimnames = list(c(),c("WOFRegg","WOFRegglowci","WOFRegghighci",
                                                               "SPSAegg","SPSAegglowci","SPSAegghighci",
                                                               "WOFRearly","WOFRearlylowci","WOFRearlyhighci",
                                                               "SPSAearly","SPSAearlylowci","SPSAearlyhighci",
                                                               "WOFRlate","WOFRlatelowci","WOFRlatehighci",
                                                               "SPSAlate","SPSAlatelowci","SPSAlatehighci",
                                                               "WaterHg")))

for (i in 1:111){
  pred.matrix[i,2] <- quantile(jagsfitpredmodel$BUGSoutput$sims.list$predWOFRegg[,i],probs = 0.025)
  pred.matrix[i,3] <- quantile(jagsfitpredmodel$BUGSoutput$sims.list$predWOFRegg[,i],probs = 0.975)
  pred.matrix[i,5] <- quantile(jagsfitpredmodel$BUGSoutput$sims.list$predSPSAegg[,i],probs = 0.025)
  pred.matrix[i,6] <- quantile(jagsfitpredmodel$BUGSoutput$sims.list$predSPSAegg[,i],probs = 0.975)
  pred.matrix[i,8] <- quantile(jagsfitpredmodel$BUGSoutput$sims.list$predWOFRearly[,i],probs = 0.025)
  pred.matrix[i,9] <- quantile(jagsfitpredmodel$BUGSoutput$sims.list$predWOFRearly[,i],probs = 0.975)
  pred.matrix[i,11] <- quantile(jagsfitpredmodel$BUGSoutput$sims.list$predSPSAearly[,i],probs = 0.025)
  pred.matrix[i,12] <- quantile(jagsfitpredmodel$BUGSoutput$sims.list$predSPSAearly[,i],probs = 0.975)
  pred.matrix[i,14] <- quantile(jagsfitpredmodel$BUGSoutput$sims.list$predWOFRlate[,i],probs = 0.025)
  pred.matrix[i,15] <- quantile(jagsfitpredmodel$BUGSoutput$sims.list$predWOFRlate[,i],probs = 0.975)
  pred.matrix[i,17] <- quantile(jagsfitpredmodel$BUGSoutput$sims.list$predSPSAlate[,i],probs = 0.025)
  pred.matrix[i,18] <- quantile(jagsfitpredmodel$BUGSoutput$sims.list$predSPSAlate[,i],probs = 0.975)
}
pred.matrix[,1] <- jagsfitpredmodel$BUGSoutput$median$predWOFRegg
pred.matrix[,4] <- jagsfitpredmodel$BUGSoutput$median$predSPSAegg
pred.matrix[,7] <- jagsfitpredmodel$BUGSoutput$median$predWOFRearly
pred.matrix[,10] <- jagsfitpredmodel$BUGSoutput$median$predSPSAearly
pred.matrix[,13] <- jagsfitpredmodel$BUGSoutput$median$predWOFRlate
pred.matrix[,16] <- jagsfitpredmodel$BUGSoutput$median$predSPSAlate
pred.matrix[,19] <- waterhgseq

##All stages
par(mfrow = c(3,2))
par(mar = c(4,5,1,1))
##Eggs
plot(exp(pred.matrix[,19]),exp(pred.matrix[,"WOFRegg"]), type = "l", xlab = "Mercury (ng/ml), water",
     ylab = "Mercury (ng/g),\n Wood Frog eggs", ylim = c(0,12))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRegglowci"]), type = "l",lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRegghighci"]), type = "l", lty = 2)
#Add observed values
par(new=T)
plot(data_subset2$Water_MeHg[data_subset2$Spp=="WOFR"&data_subset2$Life_Stage=="Eggs"],data_subset2$Amphib_MeHg[data_subset2$Spp=="WOFR"&data_subset2$Life_Stage=="Eggs"],
     xaxt="n",yaxt="n",xlab="",ylab="")
par(new=F)

plot(exp(pred.matrix[,19]),exp(pred.matrix[,"SPSAegg"]), type = "l", xlab = "Mercury (ng/ml), water",
     ylab = "Mercury (ng/g),\n Spotted Salamander eggs", ylim = c(0,12))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAegglowci"]), type = "l", lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAegghighci"]), type = "l", lty = 2)
par(new=T)
plot(data_subset2$Water_MeHg[data_subset2$Spp=="SPSA"&data_subset2$Life_Stage=="Eggs"],data_subset2$Amphib_MeHg[data_subset2$Spp=="SPSA"&data_subset2$Life_Stage=="Eggs"],
     xaxt="n",yaxt="n",xlab="",ylab="")
par(new=F)
##Early larvae
plot(exp(pred.matrix[,19]),exp(pred.matrix[,"WOFRearly"]), type = "l", xlab = "Mercury (ng/ml), water",
     ylab = "Mercury (ng/g),\n Wood Frog early larvae", ylim = c(10,375))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRearlylowci"]), type = "l",lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRearlyhighci"]), type = "l", lty = 2)
par(new=T)
plot(data_subset2$Water_MeHg[data_subset2$Spp=="WOFR"&data_subset2$Life_Stage=="Early Larvae"],data_subset2$Amphib_MeHg[data_subset2$Spp=="WOFR"&data_subset2$Life_Stage=="Early Larvae"],
     xaxt="n",yaxt="n",xlab="",ylab="")
par(new=F)

plot(exp(pred.matrix[,19]),exp(pred.matrix[,"SPSAearly"]), type = "l", xlab = "Mercury (ng/ml), water",
     ylab = "Mercury (ng/g),\n Spotted Salamander early larvae", ylim = c(10,375))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAearlylowci"]), type = "l", lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAearlyhighci"]), type = "l", lty = 2)
par(new=T)
plot(data_subset2$Water_MeHg[data_subset2$Spp=="SPSA"&data_subset2$Life_Stage=="Early Larvae"],data_subset2$Amphib_MeHg[data_subset2$Spp=="SPSA"&data_subset2$Life_Stage=="Early Larvae"],
     xaxt="n",yaxt="n",xlab="",ylab="")
par(new=F)
##Late larvae
plot(exp(pred.matrix[,19]),exp(pred.matrix[,"WOFRlate"]), type = "l", xlab = "Mercury (ng/ml), water",
     ylab = "Mercury (ng/g),\n Wood Frog late larvae", ylim = c(10,375))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRlatelowci"]), type = "l",lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRlatehighci"]), type = "l", lty = 2)
par(new=T)
plot(data_subset2$Water_MeHg[data_subset2$Spp=="WOFR"&data_subset2$Life_Stage=="Late Larvae"],data_subset2$Amphib_MeHg[data_subset2$Spp=="WOFR"&data_subset2$Life_Stage=="Late Larvae"],
     xaxt="n",yaxt="n",xlab="",ylab="")
par(new=F)

plot(exp(pred.matrix[,19]),exp(pred.matrix[,"SPSAlate"]), type = "l", xlab = "Mercury (ng/ml), water",
     ylab = "Mercury (ng/g),\n Spotted Salamander late larvae", ylim = c(10,375))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAlatelowci"]), type = "l", lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAlatehighci"]), type = "l", lty = 2)
par(new=T)
plot(data_subset2$Water_MeHg[data_subset2$Spp=="SPSA"&data_subset2$Life_Stage=="Late Larvae"],data_subset2$Amphib_MeHg[data_subset2$Spp=="SPSA"&data_subset2$Life_Stage=="Late Larvae"],
     xaxt="n",yaxt="n",xlab="",ylab="")
par(new=F)

library(ggfortify)
ggdistribution(dnorm, seq(-100, 100, 0.1), mean = 0, sd = 100)


REmodGlobal <- lmer(log(Amphib_MeHg) ~ (1|Pool) + Spp-1 + log(Water_MeHg) + Life_Stage + Habitat, 
data = data_subset2)
#install.packages("MuMIn")
library(MuMIn)
modsele <- dredge(REmodGlobal, trace = TRUE, rank = "AICc", REML = F)
fmlist <- get.models(modsele,1:15)
summary(model.avg(fmlist))
print(fmlist)
summary(modsele)
library(AICcmodavg)
rm(list = ls())
