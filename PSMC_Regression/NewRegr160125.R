ShaunRegressionInput <- read.csv("~/Desktop/Software/AlexBioProgs/160121_ErrorRegression/ShaunRegressionInput.txt", sep="")
library("ggplot2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library("RColorBrewer", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library(nlme)

Regr20Mbp = ShaunRegressionInput[ShaunRegressionInput$Bp == 20000000,] # only look at one length of genome
Regr20Mbp = Regr20Mbp[,-c(3,5,6)] # remove mutation and recomb rates, remove split but keep bp per contig
head(Regr20Mbp)
summary(Regr20Mbp$Pop_Dynamics_Type)
summary(Regr20Mbp$Bp)
ggplot(data = Regr20Mbp, aes(log(Bp_per_contig), Error, colour = Pop_Dynamics_Type)) + geom_point(size=2.5)

FH_Regr20Mbp = Regr20Mbp[Regr20Mbp$Pop_Dynamics_Type == "fauxHuman",] # only look at faxuHuman
head(FH_Regr20Mbp)
qplot(log(Error), data = FH_Regr20Mbp, bins=60) # doesn't look normal....
# individual predictors: well there's only 4 levels of Int: 10, 20, 30, 40
qplot(Bp_per_contig, data = FH_Regr20Mbp)
# a right-skewed distribution, but there's the same number of the diff Bp_per_contig
qplot(log(Bp_per_contig), data = FH_Regr20Mbp)
# now even

qplot(Int, log(Error), data = FH_Regr20Mbp) # nothin
qplot(Bp_per_contig, log(Error), data = FH_Regr20Mbp) # nothin, spread a lot
qplot(log(Bp_per_contig), log(Error), data = FH_Regr20Mbp)

# plot with colours
myColors <- brewer.pal(11,"Set3") # diff colors for diff Bp_per_contig
names(myColors) <- levels(factor(FH_Regr20Mbp$Bp_per_contig))
colScale <- scale_colour_manual(name = "Bp_per_contig",values = myColors)
Int_v_logE = ggplot(FH_Regr20Mbp, aes(x = Int, y = log(Error), colour = factor(Bp_per_contig))) + geom_point()
Int_v_logE
Int_v_logE = ggplot(FH_Regr20Mbp, aes(x = Int, y = log(Error), colour = Bp_per_contig)) + geom_point()
Int_v_logE

logBpPC_v_logE = ggplot(FH_Regr20Mbp, aes(x = log(Bp_per_contig), y = log(Error), colour = factor(Int))) + geom_point()
logBpPC_v_logE


#############33

# NoCon_Regr20Mbp = Regr20Mbp[Regr20Mbp$Pop_Dynamics_Type != "constantPop",] # forget about constant pop
#NoCon_Regr20Mbp = Regr20Mbp[Regr20Mbp$Pop_Dynamics_Type != "psmcSim2" & Regr20Mbp$Pop_Dynamics_Type != "constantPop",] # forget about sim2
NoCon_Regr20Mbp = Regr20Mbp[Regr20Mbp$Pop_Dynamics_Type != "psmcSim2",] # forget about sim2

droplevels(NoCon_Regr20Mbp$Pop_Dynamics_Type) # gets rid of the constantPop level
head(NoCon_Regr20Mbp)
summary(NoCon_Regr20Mbp$Pop_Dynamics_Type)

# remove repeated errors?
NoCon_Regr20Mbp_Uniq <- NoCon_Regr20Mbp[!duplicated(NoCon_Regr20Mbp$Error),] 
# I don't reckon that's done anything

#plot(log(NoCon_Regr20Mbp$Error))
#plot(NoCon_Regr20Mbp$Error)
# qplot()

# error vs log bp per cont
ggplot(NoCon_Regr20Mbp, aes(x = log(Bp_per_contig), y = Error, colour = Pop_Dynamics_Type) ) + geom_point(size=2.5) + geom_smooth() + facet_grid(Int ~ Pop_Dynamics_Type)
lBpPC_v_lE # there appears to be no relationship b/w Bp_per_contig, error val depends on Pop_Dynamics_Type

# error vs int
ggplot(NoCon_Regr20Mbp, aes(x = factor(Int), y = Error, colour = Pop_Dynamics_Type) ) + geom_boxplot() + facet_wrap(~Pop_Dynamics_Type, scale = "free_y")
Int_v_lE

# straight interaction
IntBp_v_lE = ggplot(NoCon_Regr20Mbp, aes(x = Int * log(Bp_per_contig), y = Error, colour = Pop_Dynamics_Type) ) + geom_point(size=2.5)
IntBp_v_lE

##### TAKEN STRAIGHT FROM JONO'S BULL EEG PAPER

### plots look crap but lets try fitting models anyway
# just fixed?
NoCon_Regr20Mbp$binvals = 

M1 <- gls(Error ~ Int * log(Bp_per_contig) , data = NoCon_Regr20Mbp, method = "REML")
summary(M1)

# random int
M2 <- lme(Error ~ Int * log(Bp_per_contig) , random = ~ 1 | Pop_Dynamics_Type, data = NoCon_Regr20Mbp, method = "REML")
summary(M2)



# random eff
M3 = lme(Error ~ 1, random = ~ 1 | Pop_Dynamics_Type, data = NoCon_Regr20Mbp, method = "REML")
#summary(M3) # WHY DOES THIS SAY DIFF NO OF OBS ON AIC???

AIC(M1, M2) 
BIC(M1, M2)
anova(M1,M2) # M2 is better by all accounts...
summary(M2)
anova(M2)
 # looks like we can remove the interaction term
M2 <- update(M2, .~. - Int:log(Bp_per_contig))
anova(M2) # all signif

# not sure what this does
as.data.frame(allEffects(M2)[[1]]) %>%
  ggplot(aes(TS, fit)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(x = "Time stamp", y = "AUC for F50 moving average")

# check residuals I guess...

mod_checks = data.frame(res = residuals(M2, type = "normalized"),
                        fit = fitted(M2),
                        lab = labels(residuals(M2)),
                        logBpPC = log(NoCon_Regr20Mbp$Bp_per_contig),
                        psmc_ints = NoCon_Regr20Mbp$Int)

p1 = ggplot(data=mod_checks, aes(fit, res, colour = lab)) + geom_point(size=1.5)
p1 # does not look particularly good...
p2 = ggplot(data=mod_checks, aes(logBpPC, res, colour = lab) ) + geom_point(size = 1.5)
p2 # bitof an upwards trend
p3 = ggplot(data=mod_checks, aes(psmc_ints, res, colour = lab) ) + geom_point(size = 1.5)
p3 # fine
p4 = ggplot(data = mod_checks, aes(sample = res)) + stat_qq()
p4 # good but for tail

summary(M2)


###################################################3
############# CHANGE POINT ######################
##################################################

# based off this: there is some kind of upward trend for small bp per contig
ggplot(NoCon_Regr20Mbp, aes(x = log(Bp_per_contig), y = Error, colour = Pop_Dynamics_Type) ) + geom_point(size=2.5) + geom_smooth() + facet_grid(Int ~ Pop_Dynamics_Type)

# maybe just try with psmcSim1 first
s1_data = NoCon_Regr20Mbp[NoCon_Regr20Mbp$Pop_Dynamics_Type == "psmcSim1",]

#binvals = log(s1_data$Bp_per_contig) < 13
#binvals = binvals + 0 # converts T.F to a numerical vector

cutoffs = c(40:80)/4
AICs = NULL
for (cutoff_val in cutoffs) {
  s1_data$binvals = log(s1_data$Bp_per_contig) < cutoff_val
  s1_data$binvals = s1_data$binvals + 0 # converts T/F to numeric
  
  chng_lme = lme(Error ~ log(Bp_per_contig) * binvals, random = ~ 1 | Int, data = s1_data)
  AICs = append(AICs, AIC(chng_lme))
}

#AIC_vs_cut = data.frame(cutoffs, AICs)
plot(cutoffs, AICs)

# cutoff of 11 looks okay
# exp(11)
binvals11 = log(s1_data$Bp_per_contig) < 11
binvals11 = binvals11 + 0
lm11 = lm(s1_data$Error ~ log(s1_data$Bp_per_contig) * binvals11)
summary(lm11)
ggplot(s1_data, aes(x = log(Bp_per_contig), y = Error, colour = factor(Int) )) + geom_point(size=2.5) + geom_smooth() + facet_grid(~ Int)
# maybe i jeust need to loo kat one int

# THIS DOESNT WORK BECAUSE OF SINGULARITIES
#################################################3
####################################################3
################################################

# C = c(40:80)/4
# for (cutoff in C) {
#   df_L = s1_data[log(s1_data$Bp_per_contig) < C,]
#   df_U = s1_data[log(s1_data$Bp_per_contig) >= C,]
#   
#   lm_L = lm(Error ~ Bp_per_contig, data = df_L)
#   lm_U = lm(Error ~ Bp_per_contig, data = df_U)
#   
# }

###########################################
tr_data = NoCon_Regr20Mbp[NoCon_Regr20Mbp$Pop_Dynamics_Type == "trench",]
ggplot(tr_data, aes(x = log(Bp_per_contig), y = Error, colour = factor(Int) )) + geom_point(size=2.5) + geom_smooth() + facet_grid(~ Int)

###############
### IDEA:
###############

# Look at the three plots with the most obvious change point, which are fH Int=10, s1 Int=20, and s1 Int=30. See if I can find a C value (change point) which is optimal for all. Procedure:
# for a C value in a range
  # fit the Error ~ Bp_per_contig * binvals model 
  # add the AICs for all three (don't need scaling as from same number of observations for all)
# find C which minimises the summed AIC

ggplot(NoCon_Regr20Mbp, aes(x = log(Bp_per_contig), y = Error, colour = Pop_Dynamics_Type) ) + geom_point(size=2.5) + geom_smooth() + facet_grid(Pop_Dynamics_Type ~ Int, scale = "free_y")

cutoffs = c(40:80)/4
sum_AICs = NULL

fHInt10 = NoCon_Regr20Mbp[NoCon_Regr20Mbp$Pop_Dynamics_Type == "fauxHuman" & NoCon_Regr20Mbp$Int == 10,]
s1Int20 = NoCon_Regr20Mbp[NoCon_Regr20Mbp$Pop_Dynamics_Type == "psmcSim1" & NoCon_Regr20Mbp$Int == 20,]
s1Int30 = NoCon_Regr20Mbp[NoCon_Regr20Mbp$Pop_Dynamics_Type == "psmcSim1" & NoCon_Regr20Mbp$Int == 30,]

for (C_val in cutoffs) {
  fHInt10_binvals = log(fHInt10$Bp_per_contig) < C_val
  fHInt10_binvals = fHInt10_binvals + 0
  lm_fHInt10 = lm(fHInt10$Error ~ log(fHInt10$Bp_per_contig) * fHInt10_binvals)
  
  s1Int20_binvals = log(s1Int20$Bp_per_contig) < C_val
  s1Int20_binvals = s1Int20_binvals + 0
  lm_s1Int20 = lm(s1Int20$Error ~ log(s1Int20$Bp_per_contig) * s1Int20_binvals)
  
  s1Int30_binvals = log(s1Int30$Bp_per_contig) < C_val
  s1Int30_binvals = s1Int30_binvals + 0
  lm_s1Int30 = lm(s1Int30$Error ~ log(s1Int30$Bp_per_contig) * s1Int30_binvals)
  
  sum_AICs = append(sum_AICs, AIC(lm_fHInt10) + AIC(lm_s1Int20) + AIC(lm_s1Int30) ) 
}
AICs_vs_Cval = data.frame(cutoffs, sum_AICs)
plot(AICs_vs_Cval$cutoffs, AICs_vs_Cval$sum_AICs) # LOOKS LIKE 11

C_val = 11

fHInt10_binvals = log(fHInt10$Bp_per_contig) < C_val
fHInt10_binvals = fHInt10_binvals + 0
lm_fHInt10 = lm(fHInt10$Error ~ log(fHInt10$Bp_per_contig) * fHInt10_binvals)

s1Int20_binvals = log(s1Int20$Bp_per_contig) < C_val
s1Int20_binvals = s1Int20_binvals + 0
lm_s1Int20 = lm(s1Int20$Error ~ log(s1Int20$Bp_per_contig) * s1Int20_binvals)

s1Int30_binvals = log(s1Int30$Bp_per_contig) < C_val
s1Int30_binvals = s1Int30_binvals + 0
lm_s1Int30 = lm(s1Int30$Error ~ log(s1Int30$Bp_per_contig) * s1Int30_binvals)

####

###
