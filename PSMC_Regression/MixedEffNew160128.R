# Goal: try a mixed effects model with lots of interactions
library("ggplot2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library(nlme)
library("MASS", lib.loc="/usr/local/lib/R/library")

data <- read.csv("~/Desktop/Software/AlexBioProgs/160121_ErrorRegression/ShaunRegressionInput.txt", sep="")

# CLEANING
data = data[data$Pop_Dynamics_Type != "psmcSim2" & data$Pop_Dynamics_Type != "constantPop",] # these dynamics are bad
data = data[data$Bp > 1000000,] # remove ones where there's not enough data
data = data[,-c(3,5,6)] # remove split, recomb and mut rate
length(data$Error)
data = data[!duplicated(data$Error),] # remove duplicated error, indicative of something stuffed

# CHANGING STUFF
data$lBpPC = log(data$Bp_per_contig) # add log of Bp_per_contig
data = data[,-3] # add Bp_per_contig cos we've got log now

data$Bp = data$Bp/1000000 # scale population by 1 million

full.lm = lm(Error ~ .^2, data=data)
step(full.lm, scope = list(upper=full.lm), data=data, direction="both")

model = lm(formula = Error ~ Bp + Int + Pop_Dynamics_Type + lBpPC + Bp:Int + Bp:Pop_Dynamics_Type + Bp:lBpPC + Int:Pop_Dynamics_Type,  data = data)
data$res = residuals(model)
data$fit = fitted(model)

ggplot(data, aes(x = fit, y = res, colour = Pop_Dynamics_Type) ) + geom_point(size=2.5) # + geom_smooth() 

ggplot(data, aes(x = Bp, y = Error, colour = Pop_Dynamics_Type) ) + geom_point(size=2.5) + geom_smooth() + facet_grid(Pop_Dynamics_Type ~ Int, scale = "free_y")

ggplot(data, aes(x = 1:length(Error), y = Error, colour = Pop_Dynamics_Type) ) + geom_boxplot()

boxcox(model)
# use a square root?

data$RtError = sqrt(data$Error)

model2 = lm(formula = RtError ~ Bp + Int + Pop_Dynamics_Type + lBpPC + Bp:Int + Bp:Pop_Dynamics_Type + Bp:lBpPC + Int:Pop_Dynamics_Type,  data = data)
summary(model2)

# CHECK RESIDUALS
data$res2 = residuals(model2)
data$fit2 = fitted(model2)
ggplot(data, aes(x = fit2, y = res2, colour = Pop_Dynamics_Type) ) + geom_point(size=2.5) # DECENT
ggplot(data = data, aes(sample = res2)) + stat_qq() # fairy linear... i guess
qplot(Bp, res2, colour = Pop_Dynamics_Type, data = data) + geom_smooth() # looks fine
qplot(Int, res2, colour = Pop_Dynamics_Type, data = data) + geom_smooth() # looks fine
qplot(Pop_Dynamics_Type, res2, colour = Pop_Dynamics_Type, data = data) + geom_smooth() # homoscedasticity fails for trench
qplot(lBpPC, res2, colour = Pop_Dynamics_Type, data = data) + geom_smooth() # looks fine
qplot(Bp*Int, res2, colour = Pop_Dynamics_Type, data = data) + geom_smooth() # more variance to the left
qplot(Bp, res2, data = data[data$Pop_Dynamics_Type == "fauxHuman",]) + geom_smooth() # looks fine
qplot(Bp, res2, data = data[data$Pop_Dynamics_Type == "psmcSim1",]) + geom_smooth() # var decr a bit left to right
qplot(Bp, res2, data = data[data$Pop_Dynamics_Type == "trench",]) + geom_smooth() # var incr a bit left to right
qplot(Bp*lBpPC, res2, colour = Pop_Dynamics_Type, data = data) + geom_smooth() # var decr a bit left to right
qplot(Int, res2, data = data[data$Pop_Dynamics_Type == "fauxHuman",]) + geom_smooth() # var decr a bit left to right
qplot(Int, res2, data = data[data$Pop_Dynamics_Type == "psmcSim1",]) + geom_smooth() # a tad wavy
qplot(Int, res2, data = data[data$Pop_Dynamics_Type == "trench",]) + geom_smooth() # looks fine
# a few things a bit off but oh well mostly okay!


# also try stepwise
full.lm2 = lm(RtError ~ .^2, data=data)
step(full.lm2, scope = list(upper=full.lm2), data=data, direction="both")
model3 = lm(formula = RtError ~ Bp + Int + Pop_Dynamics_Type + Error +              lBpPC + Bp:Pop_Dynamics_Type + Bp:Error + Int:Pop_Dynamics_Type +               Int:Error + Pop_Dynamics_Type:Error + Error:lBpPC, data = data)
summary(model3)

data$res3 = residuals(model3)
data$fit3 = fitted(model3)
ggplot(data, aes(x = fit3, y = res3, colour = Pop_Dynamics_Type) ) + geom_point(size=2.5) # CRAAAAAAAAAAAAAAAP

# WE WILL USE MODEL 2 AS A BASIS FOR MIXED EFFECTS
# gonna try the lme4 package

# random intercept model
# remove all terms with Pop_Dynamics_Type
library("lme4", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")

lme1 <- lmer(RtError ~ Bp + Int + lBpPC + Bp:Int + Bp:lBpPC + (1 | Pop_Dynamics_Type), data = data)
summary(lme1)
lme2

M1 <- gls(RtError ~ Bp + Int + Pop_Dynamics_Type + lBpPC + Bp:Int + Bp:Pop_Dynamics_Type + Bp:lBpPC + Int:Pop_Dynamics_Type,  data = data, method = "REML")
M2 <- lmer(RtError ~ Bp * Int * lBpPC + (1 | Pop_Dynamics_Type), data = data)
M3 <- lmer(RtError ~ Bp * Int * lBpPC + (1 + Bp | Pop_Dynamics_Type), data = data)
M4 <- lmer(RtError ~ Bp * Int * lBpPC + (1 + Bp + Int | Pop_Dynamics_Type), data = data)
M5 <- lmer(RtError ~ Bp * Int * lBpPC + (1 + Bp + Int + lBpPC| Pop_Dynamics_Type), data = data)
AIC(M1, M2, M3, M4, M5)

# df      AIC
# M1 13 213.2850
# M2 10 272.1443
# M3 12 253.7275
# M4 15 220.7740
# M5 19 226.2794
# Warning message:
#   In AIC.default(M1, M2, M3, M4, M5) :
#   models are not all fitted to the same number of observations

# IS THIS A PROBLEM?
# Anyway, gonna choose model M1 cos it's got the lowest AIC
# Which isn't really mixed effects?

full.gls = gls(RtError ~ (Bp + Int + Pop_Dynamics_Type + lBpPC)^2, data=data)
stepAIC(full.gls, scope = list(upper=full.gls), data=data, direction="both")

full.glm = glm(RtError ~ (Bp + Int + Pop_Dynamics_Type + lBpPC)^2, data=data)
stepAIC(full.glm, scope = list(upper=full.glm), data=data, direction="both")