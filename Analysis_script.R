data.df <- read.csv(url("https://raw.githubusercontent.com/5355693/VP_Mercury/master/hg_data.csv"))
str(data.df)
library(ggplot2)
library(plyr)
library(lubridate)
data.df$Sample_Date2 <- as.POSIXct(data.df$Sample_Date, format = "%m/%d/%y") #create a true date variable.
data.df$Year <- year(data.df$Sample_Date2) #Create a year variable to allow us to separate the 2016 observation.
data.df$individual <- seq(from = 1,to = nrow(data.df), by = 1) #create a unique identifier for each individual.

##MeHg
###No clear temporal pattern.
p1<-ggplot(data = subset(data.df, Year == 2015))
p1 + geom_smooth(aes(x = Sample_Date2, y = Water_MeHg), method = "lm") + 
  geom_point(aes(x = Sample_Date2, y = Water_MeHg, color = Pool)) + facet_wrap(~Habitat)

##Water_TotalHg
###Decreases over time. Higher in coniferous forests.
p1a <-ggplot(data = subset(data.df, Year == 2015))
p1a + geom_smooth(aes(x = Sample_Date2, y = Water_THg), method = "lm") + 
  geom_point(aes(x = Sample_Date2, y = Water_THg, color = Pool)) + facet_wrap(~Habitat)

#Water Me_mercury levels tend to be higher in Coniferous pools, but substantial variation exists among pools.
p1 + geom_line(aes(x = Sample_Date2, y = Water_MeHg, color = Pool)) + 
  geom_point(aes(x = Sample_Date2, y = Water_MeHg, color = Pool)) + facet_wrap(~Habitat) + 
  geom_smooth(aes(x = Sample_Date2, y = Water_MeHg),method = "lm") + ylab("Water methylmercury levels") + 
  xlab("Sample date") + theme(axis.text.x = element_text(angle = 340, vjust = 0.5))

#Similar for Total mercury in water.
p1a + geom_line(aes(x = Sample_Date2, y = Water_THg, color = Pool)) + 
  geom_point(aes(x = Sample_Date2, y = Water_THg, color = Pool)) + facet_wrap(~Habitat) + 
  geom_smooth(aes(x = Sample_Date2, y = Water_THg),method = "lm") + ylab("Water total mercury levels") + 
  xlab("Sample date") + theme(axis.text.x = element_text(angle = 340, vjust = 0.5))

summary(data.df)
library(R2jags)
#Data
data_subset <- subset(data.df, !duplicated(Water_MeHg))
data_subset <- subset(data_subset, (!is.na(data_subset$Water_MeHg)))
x <- data_subset$Habitat
y <- data_subset$Water_MeHg
n <- length(data_subset$Habitat)
jags.params <- c("mu1","mu2","delta","sigma","residual")
jags.inits <- function(){
  list(mu1=rnorm(1), delta=0,sigma=rlnorm(1))
}
#Model
mercuryttest <- function () {
  for(i in 1:n){
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- mu1 + delta*x[i]
  residual[i] <- y[i] - mu[1]
  }
  mu1 ~ dnorm(0,0.001)
  delta ~ dnorm (0,0.001)
  tau <- 1/(sigma*sigma)
  sigma ~ dunif(0,10)
  mu2 <- mu1 + delta #difference in mean water mercury
}

jagsfit <- jags(data = c("x","y","n"), inits = jags.inits, jags.params,
                n.iter = 5000, model.file = mercuryttest)
print(jagsfit,digits = 5)
traceplot(jagsfit)

jagsfit.mcmc <- as.mcmc(jagsfit)
plot(1:196, jagsfit$BUGSoutput$mean$residual)
boxplot(jagsfit$BUGSoutput$mean$residual ~ x)
require(lattice)

#Model without assuming equal variances:
#Data
x1 <- data_subset$Habitat[data_subset$Habitat == "Deciduous"]
x2 <- data_subset$Habitat[data_subset$Habitat == "Coniferous"]
y1 <- data_subset$Water_MeHg[data_subset$Habitat == "Deciduous"]
y2 <- data_subset$Water_MeHg[data_subset$Habitat == "Coniferous"]
n1 <- length(data_subset$Habitat[data_subset$Habitat == "Deciduous"])
n2 <- length(data_subset$Habitat[data_subset$Habitat == "Coniferous"])
jags.params <- c("mu1","mu2","delta","sigma1","sigma2")
jags.inits <- function(){
  list(mu1 = rnorm(0.001), mu2 = rnorm(0.001), sigma1 = rlnorm(1), sigma2 = rlnorm(1))
}
#Model
mercuryttest <- function () {
  for(i in 1:n1){
    y1[i] ~ dnorm(mu1, tau1)
  }
  for(i in 1:n2){
    y2[i] ~ dnorm(mu2, tau2)
  }
  mu1 ~ dnorm(0,0.001)
  mu2 ~ dnorm(0, 0.001)
  tau1 <- 1/(sigma1*sigma1)
  sigma1 ~ dunif(0,1000)
  tau2 <- 1/(sigma2*sigma2)
  sigma2 ~ dunif(0,1000)
  delta <- mu2 - mu1 #difference in mean water mercury
}

jagsfit <- jags(data = c("y1","y2","n1","n2"), inits = jags.inits, jags.params,
                n.iter = 5000, model.file = mercuryttest)
print(jagsfit,digits = 5)
traceplot(jagsfit)
library(mcmcplots)
plot(jagsfit.mcmc, trace = F)
denplot(jagsfit.mcmc,"delta", ci = 0.95)
jagsfit.mcmc <- as.mcmc(jagsfit)
densityplot(jagsfit.mcmc)
hist(jagsfit$BUGSoutput$sims.list$delta,
     xlab = "Posterior distribution of the difference in mean methylmercury\nbetween coniferous and deciudous vernal pools",
     main = NULL)
abline(v = quantile(jagsfit$BUGSoutput$sims.list$delta,probs = c(0.025,0.975)), col = "red")

##Repeat above analysis for Total mercury:
#Model without assuming equal variances:
#Data
data_subsetThg <- subset(data.df, !duplicated(Water_THg))
data_subsetThg <- subset(data_subsetThg, (!is.na(data_subsetThg$Water_THg)))
x1 <- data_subsetThg$Habitat[data_subsetThg$Habitat == "Deciduous"]
x2 <- data_subsetThg$Habitat[data_subsetThg$Habitat == "Coniferous"]
y1 <- data_subsetThg$Water_MeHg[data_subsetThg$Habitat == "Deciduous"]
y2 <- data_subsetThg$Water_MeHg[data_subsetThg$Habitat == "Coniferous"]
n1 <- length(data_subsetThg$Habitat[data_subsetThg$Habitat == "Deciduous"])
n2 <- length(data_subsetThg$Habitat[data_subsetThg$Habitat == "Coniferous"])
jags.params <- c("mu1","mu2","delta","sigma1","sigma2")
jags.inits <- function(){
  list(mu1 = rnorm(0.001), mu2 = rnorm(0.001), sigma1 = rlnorm(1), sigma2 = rlnorm(1))
}
#Model
totalmercuryttest <- function () {
  for(i in 1:n1){
    y1[i] ~ dnorm(mu1, tau1)
  }
  for(i in 1:n2){
    y2[i] ~ dnorm(mu2, tau2)
  }
  mu1 ~ dnorm(0,0.001)
  mu2 ~ dnorm(0, 0.001)
  tau1 <- 1/(sigma1*sigma1)
  sigma1 ~ dunif(0,1000)
  tau2 <- 1/(sigma2*sigma2)
  sigma2 ~ dunif(0,1000)
  delta <- mu2 - mu1 #difference in mean water mercury
}

jagsfittotalmercuryttest <- jags(data = c("y1","y2","n1","n2"), inits = jags.inits, jags.params,
                                 n.iter = 50000, model.file = totalmercuryttest)
print(jagsfittotalmercuryttest,digits = 5)
traceplot(jagsfittotalmercuryttest)
library(mcmcplots)
jagsfittotalmercuryttest.mcmc <- as.mcmc(jagsfittotalmercuryttest)
plot(jagsfittotalmercuryttest.mcmc, trace = F)
denplot(jagsfittotalmercuryttest.mcmc,"delta", ci = 0.95)
densityplot(jagsfittotalmercuryttest.mcmc)

hist(jagsfittotalmercuryttest$BUGSoutput$sims.list$delta,
     xlab = "Posterior distribution of the difference in mean total mercury\nbetween coniferous and deciudous vernal pools",
     main = NULL)
abline(v = quantile(jagsfittotalmercuryttest$BUGSoutput$sims.list$delta,probs = c(0.025,0.975)), col = "red")


##Compare mercury in amphibians among pools. Differences within-species and between habitats are nominal. 
##SPSA tend to have higher levels of MeHg, with the difference most pronounced in the deciduous habitat.

###Adults
#For adults, we need to choose which measure of mercury to use. 
#More samples are missing blood mercury (n = 17) than are missing tissue mercury (n = 3). 
#Therefore, for consistency, we should use tissue mercury as our measure.

summary(data.df$Blood_MeHg[data.df$Life_Stage=="Adult"])
summary(data.df$Tissue_MeHg[data.df$Life_Stage=="Adult"])

library(tidyr)

#In addition, mercury levels measured in blood and tissue taken from the same individual show little
#correspondence and are not interchangeable:
data.df[data.df$Life_Stage == "Adult",] %>%
gather(., key = Meas_type,value = Mercury, 9:10) %>%
  ggplot(., aes(x = Meas_type, y = Mercury, group = individual, color = Spp)) + geom_point() + geom_path() + 
  xlab("Measurement type") + ylab("Mercury level")

#Boxplots for Wood Frog by pool and habitat.
p2<-ggplot(data = subset(data.df, Year == 2015 & Spp == "WOFR"))
p2 + geom_boxplot(aes(x = Pool, y = Tissue_MeHg)) + facet_wrap(~Habitat, scales = "free_x") + 
  ylab("Methylmercury in Wood Frog tissue")

#Boxplots for Spotted Salamander by pool and habitat.
p3<-ggplot(data = subset(data.df, Year == 2015 & Spp == "SPSA"))
p3 + geom_boxplot(aes(x = Pool, y = Tissue_MeHg)) + facet_wrap(~Habitat, scales = "free_x") + 
  ylab("Methylmercury in Spotted Salamander tissue")

#Both species together, by habitat.
p4 <- ggplot(data = subset(data.df, Year == 2015))
p4 + geom_boxplot(aes(x = Spp, y = Tissue_MeHg)) + facet_wrap(~Habitat) + ylab ("Tissue methylmercury") + 
  xlab("Species")

#Not much relationship between mercury in the water and mercury in the tissue, for either species:
p5 <- ggplot(data = subset(data.df, Year == 2015 & Life_Stage == "Adult"))
p5 + geom_smooth(aes(x = Water_MeHg, y = Tissue_MeHg), method = "lm") + 
  geom_point(aes(x = Water_MeHg, y = Tissue_MeHg, color = Pool)) + facet_wrap(~Spp)


#Unclear pattern of accumulation of mercury for juvenile stages. Eggs tend to have low values for both species.
#For Wood Frog, mercury levels accumulate as might be expected if tissue mercury depends on cumulative
#exposure to water mercury. In Spotted Salamanders, however, early larvae in deciduous pools tend to have
#higher levels than late larvae.
p6 <- ggplot(data = subset(data.df, Year == 2015 & Life_Stage != "Adult"))
p6 + geom_boxplot(aes(x = Life_Stage, y = Amphib_MeHg, color = Spp)) + 
  scale_x_discrete(limits = c("Eggs","Early Larvae","Late Larvae")) + facet_wrap(~Habitat) +
  theme(axis.text.x = element_text(angle = 340, vjust = 0.5)) + 
  xlab("Life stage") + ylab("Methylmercury levels")

p7 <- ggplot(data = subset(data.df, Year == 2015 & Life_Stage != "Adult"))
p7 + geom_boxplot(aes(x = Life_Stage, y = Amphib_MeHg, color = Spp)) + 
  scale_x_discrete(limits = c("Eggs","Early Larvae","Late Larvae")) +
  theme(axis.text.x = element_text(angle = 340, vjust = 0.5)) + 
  xlab("Life stage") + ylab("Methylmercury levels")

#The odd pattern of higher mercury in early-stage larvae for SPSA doesn't appear
#to be strictly a pool effect. For example, the only two pools with early larvae
#sampled - Downer and Mauran - also show a decreased level of mercury in late-stage larvae.
p8 <- ggplot(data = subset(data.df, Year == 2015 & Life_Stage != "Adult"))
p8 + geom_point(aes(x = Life_Stage, y = Amphib_MeHg, color = Pool)) + 
  scale_x_discrete(limits = c("Eggs","Early Larvae","Late Larvae")) +
  theme(axis.text.x = element_text(angle = 340, vjust = 0.5)) + 
  xlab("Life stage") + ylab("Methylmercury levels") + facet_wrap(~Spp)

#There is no consistent relationship between water mercury levels and levels of methylmercury in the 
#tissue or blood of larval amphibians. There is a weak tendendcy for Spotted Salamander larvae to show
#a negative relationship with water mercury, but Wood Frogs show a contrasting pattern: lower levels of
#mercury in the body when water mercury is higher.
p9 <- ggplot(data = subset(data.df, Year == 2015 & Life_Stage != "Adult"))
p9 + geom_smooth(aes(x = Water_MeHg, y = Amphib_MeHg), method = "lm") + 
  geom_point(aes(x = Water_MeHg, y = Amphib_MeHg, color = Pool)) + facet_wrap(~Spp+Life_Stage, scales = "free")

p10 <- ggplot(data = subset(data.df, Year ==2015 & Life_Stage != "Adult"))
p10 + geom_smooth(aes(x = Water_MeHg, y = Amphib_MeHg), method = "lm") + 
  geom_point(aes(x = Water_MeHg, y = Amphib_MeHg, color = Life_Stage)) + facet_wrap(~Spp)

#Playing with Bayesian analysis of mercury levels in larvae.
#Dumb model: Amphib_Hg ~ Spp*Water_Hg
#Data
data_subset2 <- subset(data.df, Life_Stage != "Adult")
data_subset2 <- subset(data_subset2, (!is.na(data_subset2$Amphib_MeHg)))

spp <- data_subset2$Spp
hg <- log(data_subset2$Amphib_MeHg)
waterhg <- log(data_subset2$Water_MeHg) 
n <- nrow(data_subset2)
n.groups <- length(levels(data_subset2$Spp))
jags.params <- c("alpha","beta", "sigma")
jags.inits <- function(){
  list(sigma=rlnorm(1))
}
#Model
mercuryanova <- function () {
  for(i in 1:n){
    hg[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha[spp[i]] + beta[spp[i]]*waterhg[i]
  }
  for(i in 1:n.groups){
  alpha[i] ~ dnorm(0, 1e-11) #choice of a prior is key - if the SD of the distribution >1e-10, results bad. Assume this is
  beta [i] ~ dnorm(0, 1e-11)
  }
  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)
}

jagsfitm1 <- jags(data = c("spp","hg","waterhg","n","n.groups"), inits = jags.inits, jags.params,
                n.iter = 50000, model.file = mercuryanova)
print(jagsfitm1,digits = 5)
#traceplot(jagsfit)

##Compare with ML model. Close match, as long as careful with priors.
(mercurytest_ML <- lm(Amphib_MeHg ~ Spp*Water_MeHg, data = data_subset2))
summary(mercurytest_ML)

rm(spp,hg,waterhg,n,n.groups,jags.params,jags.inits, jagsfit)
rm(mercurytest_ML)

##Slightly smarter model: Amphib_Hg ~ Spp*Water_Hg + Life_Stage

#data_subset2 <- subset(data.df, Life_Stage != "Adult")
#data_subset2 <- subset(data_subset2, (!is.na(data_subset2$Amphib_MeHg)))

data_subset2$Life_Stage <- factor(data_subset2$Life_Stage, levels = c("Eggs", "Early Larvae", "Late Larvae", "Adult"))
spp <- data_subset2$Spp
hg <- data_subset2$Amphib_MeHg
stage <- data_subset2$Life_Stage
stage <- droplevels(stage)
waterhg <- data_subset2$Water_MeHg #could consider scaling this to get smaller Betas
n <- nrow(data_subset2)
n.groups <- length(levels(data_subset2$Spp))
n.stages <- 3
jags.params <- c("alpha","beta.spp", "beta.stage", "sigma", "mu") #added mu to monitor predictions!
jags.inits <- function(){
  list(sigma=rlnorm(1))
}
#Model
mercuryanova2 <- function () {
  for(i in 1:n){
    hg[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha[spp[i]] + beta.spp[spp[i]]*waterhg[i] + beta.stage[stage[i]]
  }
  for(i in 1:n.groups){
    alpha[i] ~ dnorm(0, 1e-11) #choice of a prior is key - if the SD of the distribution >1e-10, results bad. Assume this is
    beta.spp[i] ~ dnorm(0, 1e-11)
  }
  for(i in 1:n.stages){
    beta.stage[i] ~ dnorm(0,0.001)
  }
  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)
}

jagsfitm2 <- jags(data = c("spp","hg","waterhg","n","n.groups","stage", "n.stages"), inits = jags.inits, jags.params,
                n.iter = 20000, model.file = mercuryanova2)
print(jagsfitm2,digits = 5)
traceplot(jagsfitm2)



##Compare ML model
mercuryanova2_ML <- lm(Amphib_MeHg ~ Spp*Water_MeHg-1-Spp-Water_MeHg + Life_Stage, data = data_subset2)
summary(mercuryanova2_ML)

##Problems with prediction...getting negative values of mercury. Looking more closely at predictions, 
##these are eggs. There seems to be a boundary-space issue, where the model is struggling to 
##estimate values of mercury in eggs.
plot(predict(mercuryanova2_ML, newdata = data_subset2, type = "response"), data_subset2$Amphib_MeHg)

###Unfortunately,the BUGS output is suffering the same problem. 
plot(jagsfit$BUGSoutput$median$mu,data_subset2$Amphib_MeHg)

rm(spp, hg, stage, waterhg, n, n.groups, n.stages, jags.params, jags.inits)
rm(mercurytest2_ML)

###Try using the log of Hg levels:
(logmercuryanova2_ML <- lm(log(Amphib_MeHg) ~ Spp*log(Water_MeHg)+Life_Stage, data = data_subset2))
plot(predict(logmercuryanova2_ML, newdata = data_subset2, type = "response"),log(data_subset2$Amphib_MeHg))
summary(logmercuryanova2_ML)
#Compare adjusted r2 - this model fits much better.
abline(0,1)
plot(mercurytest2_ML)

##Try similar transform in the Bayesian model:
spp <- data_subset2$Spp
hg <- log(data_subset2$Amphib_MeHg)
stage <- data_subset2$Life_Stage
stage <- droplevels(stage)
waterhg <- log(data_subset2$Water_MeHg) #could consider scaling this to get smaller Betas
n <- nrow(data_subset2)
n.groups <- length(levels(data_subset2$Spp))
n.stages <- 3
jags.params <- c("alpha","beta.spp", "beta.stage", "sigma", "mu") #added mu to monitor predictions!
jags.inits <- function(){
  list(sigma=rlnorm(1))
}
#Model
logmercuryanova2 <- function () {
  for(i in 1:n){
    hg[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha[spp[i]] + beta.spp[spp[i]]*waterhg[i] + beta.stage[stage[i]]
  }
  for(i in 1:n.groups){
    alpha[i] ~ dnorm(0, 0.001) 
    beta.spp[i] ~ dnorm(0, 0.001)
  }
  for(i in 1:n.stages){
    beta.stage[i] ~ dnorm(0,0.001)
  }
  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)
}

jagsfitlogm2 <- jags(data = c("spp","hg","waterhg","n","n.groups","stage", "n.stages"), inits = jags.inits, jags.params,
                n.iter = 20000, model.file = logmercuryanova2)
print(jagsfitlogm2,digits = 5)
traceplot(jagsfitlogm2)
plot(jagsfit$BUGSoutput$median$mu,log(data_subset2$Amphib_MeHg))
abline(0,1)

##It is hard to compare the coefficients (and I can't figure out how to code them to be equivalent), 
#but this model is the same as the LM:
plot(jagsfit$BUGSoutput$median$mu,predict(logmercuryanova2_ML, newdata = data_subset2, type = "response"))
abline(0,1)

##Re-run with a variable to monitor predicted values and to generate stats for GOF
spp <- data_subset2$Spp
hg <- log(data_subset2$Amphib_MeHg)
stage <- data_subset2$Life_Stage
stage <- droplevels(stage)
waterhg <- log(data_subset2$Water_MeHg) 
n <- nrow(data_subset2)
n.groups <- length(levels(data_subset2$Spp))
n.stages <- 3
jags.params <- c("alpha","beta.spp", "beta.stage", "waterhg","sigma", "mu","esthg",
                 "fit","fit.new","bpvalue","residual","predicted") #added mu to monitor predictions!
jags.inits <- function(){
  list(sigma=rlnorm(1))
}
#Model Amphib_MeHg ~ Spp*Water_MeHg + Life_Stage
logmercuryanova2 <- function () {
  for(i in 1:n){
    hg[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha[spp[i]] + beta.spp[spp[i]]*waterhg[i] + beta.stage[stage[i]]
    esthg[i] <- exp(mu[i])
  }
  for(i in 1:n.groups){
    alpha[i] ~ dnorm(0, 0.001) 
    beta.spp[i] ~ dnorm(0, 0.001)
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
}

jagsfitlogm2 <- jags(data = c("spp","hg","waterhg","n","n.groups","stage", "n.stages"), inits = jags.inits, jags.params,
                n.iter = 20000, model.file = logmercuryanova2)
print(jagsfitlogm2,digits = 5)
traceplot(jagsfitlogm2)

rm(hg, n, n.groups, n.stages, jags.params, spp, stage, waterhg, jags.inits)


###Amphib_MeHg ~ Spp*Water_MeHg + Life_Stage + Habitat
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
jags.params <- c("alpha","beta.waterhg","beta.spp", "beta.stage", "beta.habitat","sigma", "mu","esthg",
                 "fit","fit.new","bpvalue","residual","predicted") #added mu to monitor predictions!
jags.inits <- function(){
  list(sigma=rlnorm(1))
}
#Model Amphib_MeHg ~ Spp*Water_MeHg + Life_Stage + Habitat
logmercuryanova3 <- function () {
  for(i in 1:n){
    hg[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha[spp[i]] + beta.spp[spp[i]]*waterhg[i] + beta.waterhg*waterhg[i] + beta.stage[stage[i]] + beta.habitat[habitat[i]]
    esthg[i] <- exp(mu[i])
  }
  for(i in 1:n.groups){
    alpha[i] ~ dnorm(0, 0.001) 
    beta.spp[i] ~ dnorm(0, 0.001)
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

jagsfitlogm3 <- jags(data = c("spp","hg","waterhg","habitat","n","n.groups","n.habitat","stage", "n.stages"), inits = jags.inits, jags.params,
                n.iter = 20000, model.file = logmercuryanova3)

print(jagsfitlogm3,digits = 5) 
print(jagsfitlogm2, digits = 5)

##The model that includes habitat has a higher DIC.
cat(paste(c("DIC for habitat model:"),jagsfitlogm3$BUGSoutput$DIC))
cat(paste(c("DIC for model w/out habitat:"), jagsfitlogm2$BUGSoutput$DIC))

traceplot(jagsfitlogm3)
plot(jagsfitlogm3$BUGSoutput$median$esthg,data_subset2$Amphib_MeHg)
abline(0,1)

##Bayesian P-value suggests a slightly worse fit (0.53):
jagsfitlogm3$BUGSoutput$mean$bpvalue

#Can visualize this graphically: not systematically under- or over-estimating error.
plot(jagsfitlogm3$BUGSoutput$sims.list$fit,jagsfitlogm3$BUGSoutput$sims.list$fit.new, main = 
       "Graphical posterior predictive check", las = 1, xlab = "SSQ for actual data set",
     ylab = "SSQ for ideal data sets")
abline(0,1)

#The ML solution also fits slightly worse than the model w/out habitat (R2 = 0.8552 v.8549)
summary(lm(log(Amphib_MeHg) ~ Spp*log(Water_MeHg) + Life_Stage + Habitat, data = data_subset2))
logmercuryanova3_ML <- lm(log(Amphib_MeHg) ~ Spp*log(Water_MeHg) + Life_Stage + Habitat, data = data_subset2)

#Confirming that the JAGS and ML models are the same:
plot(jagsfitlogm3$BUGSoutput$median$mu,predict(logmercuryanova3_ML, newdata = data_subset2, type = "response"),
     xlab = "JAGS predictions", 
     ylab = "ML predictions")
abline(0,1)

#Best model so far is Spp*WaterHg + Stage:
cat(paste(c("DIC for Spp*WaterHg:"), jagsfitm1$BUGSoutput$DIC))
cat(paste(c("DIC for Spp*WaterHg + Stage"), jagsfitlogm2$BUGSoutput$DIC))
cat(paste(c("DIC for Spp*WaterHg + Stage + Habitat:"),jagsfitlogm3$BUGSoutput$DIC))

#It would make sense to consider adding a random effect for pool, but the ML
#results suggest doing so adds nothing in terms of inference. Pool as a fixed
#effect also seems to add little information, adding some confidence to the idea
#that the effect of Pool can be subsumed by the residual variance.
library(lme4)
REmod <- lmer(log(Amphib_MeHg) ~ Spp*log(Water_MeHg) + Life_Stage + (1|Pool), 
                 data = data_subset2)
FEmod <- lm(log(Amphib_MeHg) ~ Spp*log(Water_MeHg) + Life_Stage + Pool, 
              data = data_subset2)
summary(FEmod)

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
jags.params <- c("alpha","beta.spp", "beta.stage", "waterhg","waterhgseq","sigma", "mu","predWOFRegg",
                 "predWOFRearly","predWOFRlate","predSPSAegg","predSPSAearly","predSPSAlate","fit",
                 "fit.new","bpvalue","residual","predicted") #added mu to monitor predictions!
jags.inits <- function(){
  list(sigma=rlnorm(1))
}
#Model Amphib_MeHg ~ Spp*Water_MeHg + Life_Stage
predmodel <- function () {
  for(i in 1:n){
    hg[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha[spp[i]] + beta.spp[spp[i]]*waterhg[i] + beta.stage[stage[i]]
    #esthg[i] <- exp(mu[i])
    predWOFRegg[i] <- alpha[spp[5]] + beta.spp[spp[5]]*waterhgseq[i] + beta.stage[stage[1]]
    predWOFRearly[i] <- alpha[spp[5]] + beta.spp[spp[5]]*waterhgseq[i] + beta.stage[stage[5]]
    predWOFRlate[i] <- alpha[spp[5]] + beta.spp[spp[5]]*waterhgseq[i] + beta.stage[stage[13]]
    predSPSAegg[i] <- alpha[spp[1]] + beta.spp[spp[1]]*waterhgseq[i] + beta.stage[stage[1]]
    predSPSAearly[i] <- alpha[spp[1]] + beta.spp[spp[1]]*waterhgseq[i] + beta.stage[stage[5]]
    predSPSAlate[i] <- alpha[spp[1]] + beta.spp[spp[1]]*waterhgseq[i] + beta.stage[stage[13]]
  }
  for(i in 1:n.groups){
    alpha[i] ~ dnorm(0, 0.001) 
    beta.spp[i] ~ dnorm(0, 0.001)
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
}

jagsfitpredmodel <- jags(data = c("spp","hg","waterhg","n","n.groups","stage", "n.stages","waterhgseq"), inits = jags.inits, jags.params,
                     n.iter = 20000, model.file = predmodel)
print(jagsfitpredmodel,digits = 5)

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

##Plot on same figure
par(mfrow = c(1,1))
par(mar = c(4,5,1,1))
plot(exp(pred.matrix[,19]),exp(pred.matrix[,"WOFRegg"]), type = "l", xlab = "Mercury (ng/ml), water",
     ylab = "Mercury (ng/g),\n Wood Frog eggs", ylim = c(0,12))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRegglowci"]), type = "l",lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRegghighci"]), type = "l", lty = 2)
par(new = T)
plot(exp(pred.matrix[,19]),exp(pred.matrix[,"SPSAegg"]), type = "l", xaxt="n", yaxt="n", xlab="",
     ylab="",ylim = c(0,12), col = "red")
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAegglowci"]), type = "l", lty = 2, col = "red")
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAegghighci"]), type = "l", lty = 2, col = "red")
par(new = F)

##Plot separately
par(mfrow = c(1,2))
par(mar = c(4,5,1,1))
plot(exp(pred.matrix[,19]),exp(pred.matrix[,"WOFRegg"]), type = "l", xlab = "Mercury (ng/ml), water",
     ylab = "Mercury (ng/g),\n Wood Frog eggs", ylim = c(0,12))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRegglowci"]), type = "l",lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRegghighci"]), type = "l", lty = 2)

plot(exp(pred.matrix[,19]),exp(pred.matrix[,"SPSAegg"]), type = "l", xlab = "Mercury (ng/ml), water",
     ylab = "Mercury (ng/g),\n Spotted Salamander eggs", ylim = c(0,12))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAegglowci"]), type = "l", lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAegghighci"]), type = "l", lty = 2)

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
     ylab = "Mercury (ng/g),\n Wood Frog early larvae", ylim = c(25,250))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRearlylowci"]), type = "l",lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRearlyhighci"]), type = "l", lty = 2)
par(new=T)
plot(data_subset2$Water_MeHg[data_subset2$Spp=="WOFR"&data_subset2$Life_Stage=="Early Larvae"],data_subset2$Amphib_MeHg[data_subset2$Spp=="WOFR"&data_subset2$Life_Stage=="Early Larvae"],
     xaxt="n",yaxt="n",xlab="",ylab="")
par(new=F)

plot(exp(pred.matrix[,19]),exp(pred.matrix[,"SPSAearly"]), type = "l", xlab = "Mercury (ng/ml), water",
     ylab = "Mercury (ng/g),\n Spotted Salamander early larvae", ylim = c(25,250))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAearlylowci"]), type = "l", lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAearlyhighci"]), type = "l", lty = 2)
par(new=T)
plot(data_subset2$Water_MeHg[data_subset2$Spp=="SPSA"&data_subset2$Life_Stage=="Early Larvae"],data_subset2$Amphib_MeHg[data_subset2$Spp=="SPSA"&data_subset2$Life_Stage=="Early Larvae"],
     xaxt="n",yaxt="n",xlab="",ylab="")
par(new=F)
##Late larvae
plot(exp(pred.matrix[,19]),exp(pred.matrix[,"WOFRlate"]), type = "l", xlab = "Mercury (ng/ml), water",
     ylab = "Mercury (ng/g),\n Wood Frog late larvae", ylim = c(50,300))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRlatelowci"]), type = "l",lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"WOFRlatehighci"]), type = "l", lty = 2)
par(new=T)
plot(data_subset2$Water_MeHg[data_subset2$Spp=="WOFR"&data_subset2$Life_Stage=="Late Larvae"],data_subset2$Amphib_MeHg[data_subset2$Spp=="WOFR"&data_subset2$Life_Stage=="Late Larvae"],
     xaxt="n",yaxt="n",xlab="",ylab="")
par(new=F)

plot(exp(pred.matrix[,19]),exp(pred.matrix[,"SPSAlate"]), type = "l", xlab = "Mercury (ng/ml), water",
     ylab = "Mercury (ng/g),\n Spotted Salamander late larvae", ylim = c(50,300))
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAlatelowci"]), type = "l", lty = 2)
lines(exp(pred.matrix[,"WaterHg"]),exp(pred.matrix[,"SPSAlatehighci"]), type = "l", lty = 2)
par(new=T)
plot(data_subset2$Water_MeHg[data_subset2$Spp=="SPSA"&data_subset2$Life_Stage=="Late Larvae"],data_subset2$Amphib_MeHg[data_subset2$Spp=="SPSA"&data_subset2$Life_Stage=="Late Larvae"],
     xaxt="n",yaxt="n",xlab="",ylab="")
par(new=F)
