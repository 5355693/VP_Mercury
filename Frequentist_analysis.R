data.df <- read.csv(url("https://raw.githubusercontent.com/5355693/VP_Mercury/master/hg_data.csv"))
str(data.df)
#library(ggplot2)
#library(plyr)
#library(lubridate)
data.df$Sample_Date2 <- as.POSIXct(data.df$Sample_Date, format = "%m/%d/%y") #create a true date variable.
data.df$Year <- year(data.df$Sample_Date2) #Create a year variable to allow us to separate the 2016 observation.
data.df$individual <- seq(from = 1,to = nrow(data.df), by = 1) #create a unique identifier for each individual.
data.df$Life_Stage <- factor(data.df$Life_Stage, levels = c("Eggs","Early Larvae", "Late Larvae", "Adult"))
summary(data.df)

data_juvenile <- subset(data.df, Life_Stage != "Adult")
data_juvenile <- subset(data_juvenile, (!is.na(data_juvenile$Amphib_MeHg)))
data_juvenile <- droplevels(data_juvenile)
data_juvenile$Life_Stage <- as.factor(ifelse(data_juvenile$Life_Stage == "Late Larvae", "LateLarvae",
                                   ifelse(data_juvenile$Life_Stage == "Early Larvae", "EarlyLarvae","Eggs")))

data_juvenile$lnAmphib_MeHg <- log(data_juvenile$Amphib_MeHg)
data_juvenile$lnWater_MeHg <-  log(data_juvenile$Water_MeHg)
##Set "Egg" as the reference category for Life_Stage
data_juvenile <- within(data_juvenile, Life_Stage <- relevel(Life_Stage, ref = "Eggs"))
levels(data_juvenile$Life_Stage)

library(AICcmodavg)

cand.Mod <- list()
##global model
cand.Mod[[1]] <- lmer(lnAmphib_MeHg ~ Spp + Life_Stage + Spp:Life_Stage + lnWater_MeHg + Habitat + (1|Pool),
                      data = data_juvenile, REML = F)
##Spp + Life_Stage + WaterMeHg
cand.Mod[[2]] <- lmer(lnAmphib_MeHg ~ Spp + Life_Stage + lnWater_MeHg + (1|Pool),
                      data = data_juvenile,REML = F)
##Spp + Life_Stage + Habitat
cand.Mod[[3]] <- lmer(lnAmphib_MeHg ~ Spp + Life_Stage + Habitat + (1|Pool),
                      data = data_juvenile, REML = F)
##Spp + Habitat + WaterMeHg
cand.Mod[[4]] <- lmer(lnAmphib_MeHg ~ Spp + Habitat + lnWater_MeHg + (1|Pool),
                      data = data_juvenile, REML = F)
##Life_Stage + Habitat + WaterMeHg
cand.Mod[[5]] <- lmer(lnAmphib_MeHg ~ Life_Stage + Habitat + lnWater_MeHg + (1|Pool),
                      data = data_juvenile, REML = F)
##Spp + Life_Stage
cand.Mod[[6]] <- lmer(lnAmphib_MeHg ~ Spp + Life_Stage + (1|Pool),
                      data = data_juvenile, REML = F)
##Spp + WaterMeHg
cand.Mod[[7]] <- lmer(lnAmphib_MeHg ~ Spp + lnWater_MeHg + (1|Pool),
                      data = data_juvenile, REML = F)
##Spp + Habitat
cand.Mod[[8]] <- lmer(lnAmphib_MeHg ~ Spp + Habitat + (1|Pool),
                      data = data_juvenile, REML = F)
##Life_stage + Water_MeHg
cand.Mod[[9]] <- lmer(lnAmphib_MeHg ~ Life_Stage + lnWater_MeHg + (1|Pool),
                      data = data_juvenile, REML = F)
##Life_Stage + Habitat
cand.Mod[[10]] <- lmer(lnAmphib_MeHg ~ Life_Stage + Habitat + (1|Pool),
                      data = data_juvenile, REML = F)
##Water_MeHg + Habitat
cand.Mod[[11]] <- lmer(lnAmphib_MeHg ~ lnWater_MeHg + Habitat + (1|Pool),
                       data = data_juvenile, REML = F)
##Spp
cand.Mod[[12]] <- lmer(lnAmphib_MeHg ~ Spp + (1|Pool),
                       data = data_juvenile, REML = F)
##Life_Stage
cand.Mod[[13]] <- lmer(lnAmphib_MeHg ~ Life_Stage + (1|Pool),
                       data = data_juvenile, REML = F)
##Water_MeHg
cand.Mod[[14]] <- lmer(lnAmphib_MeHg ~ lnWater_MeHg + (1|Pool),
                       data = data_juvenile, REML = F)
##Habitat
cand.Mod[[15]] <- lmer(lnAmphib_MeHg ~ Habitat + (1|Pool),
                       data = data_juvenile, REML = F)
##Spp * Life_Stage + WaterMeHg
cand.Mod[[16]] <- lmer(lnAmphib_MeHg ~ Spp + Life_Stage + Spp:Life_Stage + lnWater_MeHg + (1|Pool),
                      data = data_juvenile,REML = F)
##Spp * Life_Stage + Habitat
cand.Mod[[17]] <- lmer(lnAmphib_MeHg ~ Spp + Life_Stage + Spp:Life_Stage + Habitat + (1|Pool),
                      data = data_juvenile, REML = F)
##Spp * Life_Stage
cand.Mod[[18]] <- lmer(lnAmphib_MeHg ~ Spp + Life_Stage + Spp:Life_Stage + (1|Pool),
                      data = data_juvenile, REML = F)

modNames <- c("Global","Spp + Life_Stage + WaterMeHg","Spp + Life_Stage + Habitat",
              "Spp + Habitat + WaterMeHg","Life_Stage + Habitat + WaterMeHg",
              "Spp + Life_Stage","Spp + WaterMeHg","Spp + Habitat","Life_Stage + Water_MeHg",
              "Life_Stage + Habitat","Water_MeHg + Habitat","Spp","Life_Stage","Water_MeHg",
              "Habitat", "Species*Life_Stage + Water_MeHg", "Species*Life_Stage + Habitat",
              "Species*Life_stage")
aictab(cand.set = cand.Mod, modnames = modNames, nobs = NULL, second.ord = F)
print(modavg(parm = "lnWater_MeHg",cand.set = cand.Mod, modnames = modNames))
print(modavg(parm = "SppWOFR",cand.set = cand.Mod, modnames = modNames,
             exclude = list("Spp:Life_Stage")))
print(modavg(parm = "HabitatDeciduous", cand.set = cand.Mod, modnames = modNames))
print(modavg(parm = "Life_StageLateLarvae", cand.set = cand.Mod, modnames = modNames,
             exclude = list("Spp:Life_Stage")))
print(modavg(parm = "Life_StageEarlyLarvae", cand.set = cand.Mod, modnames = modNames,
             exclude = list("Spp:Life_Stage")))

##Examining the global model. Note that the effect of pool is approximately zero:
## From a question about this on StackExchange: 
## "the extent of this subject variation can be fully or virtually-fully explained by 
## just the residual variance term alone. There is not enough additional subject-level variation 
## to warrant adding an additional subject-level random effect to explain all the observed variation."
## Ben Bolker says this is common, especially with small numbers of random factors. Leaving it in 
## or removing it makes no difference to tests of the effects of fixed effects, although it will effect
## AIC values b/c the variance parameter is counted...so will make probably tend to favor simpler models.
## "If one chooses for philosophical grounds to retain these parameters, it won't change any of the answers".
global.mod <- lmer(lnAmphib_MeHg ~ Spp + Life_Stage + Spp:Life_Stage + lnWater_MeHg + Habitat + (1|Pool),
                   data = data_juvenile, REML = F)
summary(global.mod)

##Residuals are centered on zero and error is largely similar across levels of the random effect:
plot(global.mod, Pool ~ resid(.), abline = 0)
#Residuals from the final weighted-least-squares regression of the IWLS procedure used to fit the model
##useful, for example, for detecting nonlinearity. 
plot(global.mod, resid(., type = "working") ~ fitted(.)|Pool, abline = 0)
plot(global.mod, resid(., type = "pearson") ~ fitted(.)|Pool, abline = 0)
plot(global.mod, resid(., type = "deviance") ~ fitted(.)|Pool, abline = 0)

##Normality of errors: symmetrical around zero, with heavier tails.
##Heavier tails tend to inflate estimate of within-group error,leading to 
##more conservative tests of fixed effects.
qqnorm(global.mod, ~resid(.))
abline(0,1)

###Errors distributed around zero with approximately constant variance
plot(global.mod)
##Look at independence:
  # i. plot residuals vs each covariate in the model
plot(global.mod,resid(.) ~lnWater_MeHg|Pool)
plot(global.mod, Spp ~ resid(.), abline = 0)
plot(global.mod, Habitat ~ resid(.), abline = 0)
plot(global.mod, Life_Stage ~ resid(.), abline = 0)

##Look at predicted vs. observed
plot(predict(global.mod),data_juvenile$lnAmphib_MeHg)
abline(0,1)
summary(global.mod)
r.squaredGLMM(global.mod)

##Generate predictions
##Full set:
#newdata <- data.frame(
#  lnWater_MeHg = rep(c(-7.9,-7.59,-7.38,-7.17,-6.96,-6.75,-6.54,-6.33,-6.12,-5.91),12),
#  Spp = rep(c(rep("WOFR",10),rep("SPSA",10)),6),
#  Life_Stage = rep(c(rep("Eggs",20),rep("EarlyLarvae",20), rep("LateLarvae",20)),2),
#  Habitat = rep(c(rep("Deciduous",60), rep("Coniferous",60)))
#)

#Holding unimportant variables constant:
newdata <- data.frame(
  lnWater_MeHg = rep(mean(data_juvenile$lnWater_MeHg),120),
  Spp = rep(c(rep("WOFR",10),rep("SPSA",10)),6),
  Life_Stage = rep(c(rep("Eggs",20),rep("EarlyLarvae",20), rep("LateLarvae",20)),2),
  Habitat = rep("Deciduous",120)
)


modAvPreds <- as.data.frame(modavgPred(cand.set = cand.Mod, modnames = modNames, newdata = newdata, type = "response",
           uncond.se = "revised"), digits = 4)
modAvPreds <- cbind(newdata,modAvPreds)
modAvPreds$bt.mod.avg.pred <- exp(modAvPreds$mod.avg.pred)
modAvPreds$bt.lower.CL <- exp(modAvPreds$lower.CL)
modAvPreds$bt.upper.CL <- exp(modAvPreds$upper.CL)

##Plot them
modAvPreds$Life_Stage <- as.factor(ifelse(modAvPreds$Life_Stage == "LateLarvae", "Late Larvae",
                                             ifelse(modAvPreds$Life_Stage == "EarlyLarvae", "Early Larvae","Eggs")))
modAvPreds$Spp <- as.factor(ifelse(modAvPreds$Spp == "SPSA", "Spotted Salamander","Wood Frog"))
modAvPreds <- within(modAvPreds, Life_Stage <- relevel(Life_Stage, ref = "Eggs"))
colnames(modAvPreds)[2] <- "Species"
p <- ggplot(aes(x = Life_Stage, y = mod.avg.pred, fill = Species), data = modAvPreds) + 
  geom_bar(stat = "identity",color = "black", position = position_dodge()) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, 
                position = position_dodge(0.9)) + 
  labs(x = "Life stage", y = "ln methlymercury (ng/g)")
p + theme_classic() + scale_fill_manual(values = c('#999999','808080'))


#Back-transformed MeHg
p <- ggplot(aes(x = Life_Stage, y = bt.mod.avg.pred, fill = Species), data = modAvPreds) + 
  geom_bar(stat = "identity",color = "black", position = position_dodge()) + 
  geom_errorbar(aes(ymin = bt.lower.CL, ymax = bt.upper.CL), width = 0.2, 
                position = position_dodge(0.9)) + 
  labs(x = "Life stage", y = "Methlymercury (ng/g)")
p + theme_classic() + scale_fill_manual(values = c('#999999','808080'))

##Formal analysis for adults. There are really only 2 sensible models:
##1) Tissue_MeHg ~ Species
##1) Tissue_MeHg ~ Habitat
##3) Tissue_MeHg ~ Species + Habitat

#For adults, we need to choose which measure of mercury to use. 
#More samples are missing blood mercury (n = 17) than are missing tissue mercury (n = 3). 
#Therefore, for consistency, we should use tissue mercury as our measure.
data_adult <- subset(data.df, Life_Stage == "Adult")
data_adult <- subset(data_adult, (!is.na(data_adult$Tissue_MeHg)))
data_adult <- droplevels(data_adult)
summary(data_adult)

cand.ModAd <- list()
##global model
cand.ModAd[[1]] <- lmer(log(Tissue_MeHg) ~ Spp + Habitat + (1|Pool),
                      data = data_adult, REML = F)
##Spp 
cand.ModAd[[2]] <- lmer(log(Tissue_MeHg) ~ Spp + (1|Pool),
                      data = data_adult,REML = F)
##Habitat
cand.ModAd[[3]] <- lmer(log(Tissue_MeHg) ~ Habitat + (1|Pool),
                      data = data_adult, REML = F)

modNamesAd <- c("Global","Spp","Habitat")
aictab(cand.set = cand.ModAd, modnames = modNamesAd, nobs = NULL, second.ord = F)
print(modavgShrink(parm = "HabitatDeciduous",cand.set = cand.ModAd, modnames = modNamesAd))
print(modavgShrink(parm = "SppWOFR",cand.set = cand.ModAd, modnames = modNamesAd))

global.modAd <- lmer(log(Tissue_MeHg) ~ Spp + Habitat + (1|Pool),
                     data = data_adult, REML = F)
summary(global.modAd)
##Residuals are centered on zero and error is largely similar across levels of the random effect.
## The Podunk pools look to have residuals shifted low, and have outliers, but it's hard to tell
## if this is real or the result of very small sample size. 
plot(cand.ModAd[[1]], Pool ~ resid(.), abline = 0)
plot(cand.ModAd[[1]], fitted(.) ~ resid(.)|Habitat)
#Residuals from the final weighted-least-squares regression of the IWLS procedure used to fit the model
##useful, for example, for detecting nonlinearity. 
plot(cand.ModAd[[1]], resid(., type = "working") ~ fitted(.)|Pool, abline = 0)
plot(cand.ModAd[[1]], resid(., type = "pearson") ~ fitted(.)|Pool, abline = 0)
plot(cand.ModAd[[1]], resid(., type = "deviance") ~ fitted(.)|Pool, abline = 0)

##Normality of errors: symmetrical around zero, with heavier tails.
##Heavier tails tend to inflate estimate of within-group error,leading to 
##more conservative tests of fixed effects. Again, the symmetry of the
##heavy tails indicates little effect on fixed effects is likely.
qqnorm(resid(cand.ModAd[[1]]), abline(0,1))
abline(0,1)

###Errors distributed around zero with approximately constant variance
plot(cand.ModAd[[1]])
##Look at independence:
# i. plot residuals vs each covariate in the model
plot(cand.ModAd[[1]], Spp ~ resid(.), abline = 0)
plot(cand.ModAd[[1]], Habitat ~ resid(.), abline = 0)

##Look at predicted vs. observed
plot(predict(cand.ModAd[[1]]),log(data_adult$Tissue_MeHg))
abline(0,1)

## This model predicts fairly poorly:
r.squaredGLMM(cand.ModAd[[1]]) #R2c = 0.274

##Generate predictions
##Full set:
newdataAd <- data.frame(
  Spp = c(rep("WOFR",2),rep("SPSA",2)),
  Habitat = c(rep("Deciduous",2),rep("Coniferous",2))
)

modAvPredsAd <- as.data.frame(modavgPred(cand.set = cand.ModAd, modnames = modNamesAd, newdata = newdataAd,
                                         type = "response",
                                       uncond.se = "revised"), digits = 4)
modAvPredsAd <- cbind(newdataAd,modAvPredsAd)
modAvPredsAd$Spp <- as.factor(ifelse(modAvPredsAd$Spp == "SPSA", "Spotted Salamander","Wood Frog"))
colnames(modAvPredsAd)[1] <- "Species"
modAvPredsAd <- within(modAvPredsAd, Species <- relevel(Species, ref = "Spotted Salamander"))


ggplot(aes(x = Species, y = mod.avg.pred, fill = Species), data = modAvPredsAd) + 
  geom_bar(stat = "identity",color = "black", position = position_dodge()) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, 
                position = position_dodge(0.9)) + labs(x = "Species", y = "ln methlymercury (ng/g)") + 
  theme_classic() + scale_fill_manual(values = c('#999999','808080'))
library(ggplot2)
