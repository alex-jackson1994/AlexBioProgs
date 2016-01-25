ShaunRegressionInput <- read.csv("~/Desktop/Software/AlexBioProgs/160121_ErrorRegression/ShaunRegressionInput.txt", sep="")
library("ggplot2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library("RColorBrewer", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")

Regr20Mbp = ShaunRegressionInput[ShaunRegressionInput$Bp == 20000000,] # only look at one length of genome
Regr20Mbp = Regr20Mbp[,-c(3,5,6)] # remove mutation and recomb rates, remove split but keep bp per contig
head(Regr20Mbp)
summary(Regr20Mbp$Pop_Dynamics_Type)
summary(Regr20Mbp$Bp)

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
NoCon_Regr20Mbp = Regr20Mbp[Regr20Mbp$Pop_Dynamics_Type != "psmcSim2",] # forget about sim2
droplevels(NoCon_Regr20Mbp$Pop_Dynamics_Type) # gets rid of the constantPop level
head(NoCon_Regr20Mbp)
levels(NoCon_Regr20Mbp$Pop_Dynamics_Type)


#plot(log(NoCon_Regr20Mbp$Error))
#plot(NoCon_Regr20Mbp$Error)
# qplot()

lBpPC_v_lE = ggplot(NoCon_Regr20Mbp, aes(x = log(Bp_per_contig), y = log(Error), colour = Pop_Dynamics_Type) ) + geom_point(size=2.5)
lBpPC_v_lE # there appears to be no relationship b/w Bp_per_contig, error val depends on Pop_Dynamics_Type

Int_v_lE = ggplot(NoCon_Regr20Mbp, aes(x = Int, y = log(Error), colour = Pop_Dynamics_Type) ) + geom_point(size=2.5)
Int_v_lE
