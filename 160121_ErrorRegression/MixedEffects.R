ShaunRegressionInput <- read.csv("~/Desktop/Analyses/160121_ErrorRegression/ShaunRegressionInput.txt", sep="")

ShaunRegressionInputOrd = ShaunRegressionInput[order(ShaunRegressionInput$Error, decreasing=TRUE),]
head(ShaunRegressionInputOrd, n=15)

# gonna remove the first 11 lines, which contain huge errors

ShaunRegressionInputOrdCleaned = ShaunRegressionInputOrd[-c(1:11),]
head(ShaunRegressionInputOrdCleaned)

# need to remove one of recombination or mutation rate because they are perfectly collinear... actually remove both, they both depend upon Bp
# also remove split, just keep bp per contig
ShaunRegressionInputOrdCleanedRemov_Mut = ShaunRegressionInputOrdCleaned[,-c(3,5,6)]
head(ShaunRegressionInputOrdCleanedRemov_Mut)

#Removing duplicated errors, probably indicitave of a dodgy response
ShaunRegressionInputOrdCleanedRemov_MutUnique<-ShaunRegressionInputOrdCleanedRemov_Mut[!duplicated(ShaunRegressionInputOrdCleanedRemov_Mut$Error),]
head(ShaunRegressionInputOrdCleanedRemov_MutUnique)

#####################################
####### MIXED EFFECTS MODELLING #####
#####################################
library(nlme)

class(ShaunRegressionInputOrdCleanedRemov_MutUnique$Pop_Dynamics_Type) # it's a factor. good, need this for lme

# no random part
M1 = gls(log(Error) ~ 1 + Bp * Int * Bp_per_contig, method = "REML", data = ShaunRegressionInputOrdCleanedRemov_MutUnique)
# random intercept
M2 = lme(log(Error) ~ Bp * Int, random = ~ 1 | Pop_Dynamics_Type, method = "REML", data = ShaunRegressionInputOrdCleanedRemov_MutUnique)
# random slope
M3 = lme(Error ~ Bp * Int, random = ~ 1 + Bp | Pop_Dynamics_Type, method = "REML", data = ShaunRegressionInputOrdCleanedRemov_MutUnique)
M4 = lme(Error ~ Bp * Int, random = ~ 1 + Bp + Int | Pop_Dynamics_Type, method = "REML", data = ShaunRegressionInputOrdCleanedRemov_MutUnique)
M5 = lme(Error ~ Bp * Int, random = ~ 1 + Bp + Int + Bp_per_contig | Pop_Dynamics_Type, method = "REML", data = ShaunRegressionInputOrdCleanedRemov_MutUnique)

# Error in solve.default(estimates[dimE[1L] - (p:1), dimE[2L] - (p:1), drop = FALSE]) : 
# System is computationally singular: reciprocal condition number

BIC(M1, M2)
anova(M1,M2)
