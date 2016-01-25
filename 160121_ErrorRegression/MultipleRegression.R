ShaunRegressionInput <- read.csv("~/Desktop/Software/Alex/160121_ErrorRegression/ShaunRegressionInput.txt", sep="")

ShaunRegressionInputOrd = ShaunRegressionInput[order(ShaunRegressionInput$Error, decreasing=TRUE),]
head(ShaunRegressionInputOrd, n=15)

# gonna remove the first 11 lines, which contain huge errors

ShaunRegressionInputOrdCleaned = ShaunRegressionInputOrd[-c(1:11),]
head(ShaunRegressionInputOrdCleaned)

# need to remove one of recombination or mutation rate because they are perfectly collinear... actually remove both, they both can be found straight from the Bp
# also remove split, just keep bp per contig (split and bp give bp per contig)
ShaunRegressionInputOrdCleanedRemov_Mut = ShaunRegressionInputOrdCleaned[,-c(3,5,6)]
head(ShaunRegressionInputOrdCleanedRemov_Mut)

# Removing duplicated errors, probably indicitave of a dodgy response from PSMC (e.g. just guesses population sizes of 1 all the time)
ShaunRegressionInputOrdCleanedRemov_MutUnique<-ShaunRegressionInputOrdCleanedRemov_Mut[!duplicated(ShaunRegressionInputOrdCleanedRemov_Mut$Error),]
#ShaunRegressionInputOrdCleanedRemov_MutUnique = ShaunRegressionInputOrdCleanedRemov_Mut

# add a log error column
ShaunRegressionInputOrdCleanedRemov_MutUnique$LogError = log(ShaunRegressionInputOrdCleanedRemov_MutUnique$Error)
head(ShaunRegressionInputOrdCleanedRemov_MutUnique)

##############################################
######### MODELLING ##########################
##############################################

# Ben's suggestions:
# look at correlation between vars
# pairwise plots, with colors for diff pop dynamics
# note errors can be forced by params

# compares individual predictors with logError. Produces 5 sets of plots, one for each type of population dynamic
par(mfrow=c(1,3))
for (lvl in levels(ShaunRegressionInputOrdCleanedRemov_MutUnique$Pop_Dynamics_Type)) {
  x = lvl
  plot(ShaunRegressionInputOrdCleanedRemov_MutUnique$Bp[ShaunRegressionInputOrdCleanedRemov_MutUnique$Pop_Dynamics_Type == x],ShaunRegressionInputOrdCleanedRemov_MutUnique$LogError[ShaunRegressionInputOrdCleanedRemov_MutUnique$Pop_Dynamics_Type == x], xlab= "Bp", ylab = "logError")
  plot(ShaunRegressionInputOrdCleanedRemov_MutUnique$Int[ShaunRegressionInputOrdCleanedRemov_MutUnique$Pop_Dynamics_Type == x],ShaunRegressionInputOrdCleanedRemov_MutUnique$LogError[ShaunRegressionInputOrdCleanedRemov_MutUnique$Pop_Dynamics_Type == x], main = x, xlab= "Int", ylab = "")
  plot(ShaunRegressionInputOrdCleanedRemov_MutUnique$Bp_per_contig[ShaunRegressionInputOrdCleanedRemov_MutUnique$Pop_Dynamics_Type == x],ShaunRegressionInputOrdCleanedRemov_MutUnique$LogError[ShaunRegressionInputOrdCleanedRemov_MutUnique$Pop_Dynamics_Type == x], xlab= "Bp per contig", ylab = "")
}

# too messy:
# pairs(ShaunRegressionInputOrdCleanedRemov_MutUnique[c(1:3,6)], pch = 21, bg = c("red", "green", "blue", "brown", "yellow")[unclass(ShaunRegressionInputOrdCleanedRemov_MutUnique$Pop_Dynamics_Type)])
# I think this is constant = red, faux = green, sim1 =  blue, sim2 = brown, trench = yellow

### IGNORE BELOW THIS


#require(ggplot2)



# create range of models
null.lm = lm(Error ~ 1, data=ShaunRegressionInputOrdCleanedRemov_Mut)
summary(null.lm)
full.lm = lm(Error ~ ., data=ShaunRegressionInputOrdCleanedRemov_Mut)
summary(full.lm)

# perform forward stepwise regression according to http://www.stat.columbia.edu/~martin/W2024/R10.pdf
step(null.lm, scope = list(upper=full.lm), data=ShaunRegressionInputOrdCleanedRemov_Mut, direction="both")

# we get this
m1 = lm(formula = Error ~ Pop_Dynamics_Type + Bp, data = ShaunRegressionInputOrdCleanedRemov_Mut)

# m1 ASSUMPTION CHECKING: LINEARITY
m1.res = residuals(m1)
m1.fit = fitted(m1)
plot(m1.fit,m1.res)

#Removing duplicated errors, probably indicitave of a dodgy response
ShaunRegressionInputOrdCleanedRemov_MutUnique<-ShaunRegressionInputOrdCleanedRemov_Mut[!duplicated(ShaunRegressionInputOrdCleanedRemov_Mut$Error),]

null.lm = lm(Error ~ 1, data=ShaunRegressionInputOrdCleanedRemov_MutUnique)
summary(null.lm)
full.lm = lm(Error ~ ., data=ShaunRegressionInputOrdCleanedRemov_MutUnique)
summary(full.lm)

# perform forward stepwise regression according to http://www.stat.columbia.edu/~martin/W2024/R10.pdf
step(null.lm, scope = list(upper=full.lm), data=ShaunRegressionInputOrdCleanedRemov_MutUnique, direction="both")
step(full.lm, scope = list(upper=null.lm), data=ShaunRegressionInputOrdCleanedRemov_MutUnique, direction="both")

# we get this
m2 = lm(formula = Error ~ Pop_Dynamics_Type + Bp + Int, data = ShaunRegressionInputOrdCleanedRemov_MutUnique)

# m1 ASSUMPTION CHECKING: LINEARITY
m2.res = residuals(m1)
m2.fit = fitted(m1)
plot(m2.fit,m2.res)
qqnorm(m2.res)

# try a log transformation
m3 = lm(formula = log(Error) ~ Pop_Dynamics_Type + Bp + Int, data = ShaunRegressionInputOrdCleanedRemov_MutUnique)

# m1 ASSUMPTION CHECKING: LINEARITY
m3.res = residuals(m1)
m3.fit = fitted(m1)
plot(m3.fit,m3.res)
qqnorm(m3.res)

#Consider only fauxhuman
fauxHumanData<-ShaunRegressionInput[which(ShaunRegressionInput$Pop_Dynamics_Type=='fauxHuman'),]
fauxHumanDataUnique<-fauxHumanData[!duplicated(fauxHumanData$Error),]

null.lm = lm(Error ~ 1, data=fauxHumanDataUnique)
summary(null.lm)
full.lm = lm(Error ~ Bp+Int+Split+Bp_per_contig, data=fauxHumanDataUnique)
summary(full.lm)
plot(fitted(full.lm),residuals(full.lm))
logFull.lm = lm(log(Error) ~ Bp+Int+Bp_per_contig, data=fauxHumanDataUnique)
summary(logFull.lm)
plot(fitted(logFull.lm),residuals(logFull.lm))
plot(fauxHumanDataUnique$Bp_per_contig,fauxHumanDataUnique$Error)


#Consider only fauxhuman
psmcSim1Data<-ShaunRegressionInput[which(ShaunRegressionInput$Pop_Dynamics_Type=='psmcSim1'),]
psmcSim1DataUnique<-psmcSim1Data[!duplicated(psmcSim1Data$Error),]

null.lm = lm(Error ~ 1, data=psmcSim1DataUnique)
summary(null.lm)
full.lm = lm(Error ~ Bp+Int+Split+Bp_per_contig, data=psmcSim1DataUnique)
summary(full.lm)
plot(fitted(full.lm),residuals(full.lm))
logFull.lm = lm(log(Error) ~ Bp+Int, data=psmcSim1DataUnique)
summary(logFull.lm)
plot(fitted(logFull.lm),residuals(logFull.lm))
plot(psmcSim1DataUnique$Bp_per_contig,psmcSim1DataUnique$Error)