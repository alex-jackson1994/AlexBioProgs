library(ggplot2)

# the redKangarooInt4+5*3+4.psmc run
kanga1 <- read.delim("~/Desktop/Data/redKangaroo/redKangarooInt4+5*3+4.txt")
# the redKangarooInt40.psmc run
# THIS ONE JUST GAVE CONSTANT POP
# kanga2 = read.delim("/media/alex/My Passport/bioinformaticsScholarship/Data/redKangarooInt40.txt")


# recall that time (x) must be divided by 2 to give in terms on N_0 generations having past
# taking log of time for clarity
#ggplot(data=kanga1, mapping=aes(x = log(t_k/2), y = lambda_k)) + geom_step(colour="red") + labs(xlab("log(t_k/2) (log N_0 generations)"), ylab("lambda_k")) + ggtitle("History of the red kangaroo's population dynamics") # + geom_step(data = kanga2, mapping=aes(x=log(t_k/2), y=lambda_k), colour="black")
#ggplot(data=kanga1, mapping=aes(x = t_k/2, y = lambda_k)) + geom_step(colour="red") + labs(xlab("t_k/2 (N_0 generations)"), ylab("lambda_k")) + ggtitle("History of the red kangaroo's population dynamics")

# Looking at the PSMC Paper's methods section we find N_0 = Theta/4Mu or see README
# Assume we know the per-site per-generation mutation rate as 2.5*10^(-8)
# Then (from README):
# Firstly, suppose we know the per-site per-generation mutation rate \mu, we can
# compute N_0 as:
#   
#   N_0 = \theta_0 / (4\mu) / s
# 
# where \theta_0 is given at the 2nd column of "TR" lines, and s is the bin size
# we use for generating the PSMC input.
theta_0 = 0.018540 # from the file redKangarooInt4+5*3+4.psmc. We take theta from the final iteration.
mu = 2.5*10^(-8)
s = 100
N_0 = theta_0 / (4*mu) / s # = 1854 apparently...

# Generation time: 7-10 years (Dawson 2012) http://www.environment.nsw.gov.au/resources/threatenedspecies/determinations/PDRedKangReject.pdf
gen_low = 7
gen_hi = 10

kanga1$t_unscaled_low = kanga1$t_k*2 * N_0 * gen_low
kanga1$t_unscaled_hi = kanga1$t_k*2 * N_0 * gen_hi
kanga1$pop_hist = kanga1$lambda_k * N_0

# when was the LGM... according to wikipedia lol,
LGM_max = 26500 # time for maximum pos the ice sheets reached
LGM_deg = 14500 # deglaciation
IA_start = 110000 # Ice Age, not LGM... from https://en.wikipedia.org/wiki/Last_glacial_period

cutoff7gen = c(20000/25*7, 3000000/25*7) # cutoffs of PSMC "goodness" based upon Heng Li's 20 kyr and 3 Myr estimates, but scaled to 7 year generations instead of 25 year
cutoff10gen = c(20000/25*10, 3000000/25*10)

# bit of a mess
ggplot(data=kanga1, mapping=aes(x = log(t_unscaled_low, base=10), y = pop_hist)) + geom_step(colour="red") + labs(xlab("log10(years in past)"), ylab("population size")) + ggtitle("The red kangaroo's population dynamics, and the LGM") + geom_step(data = kanga1, mapping=aes(x=log(t_unscaled_hi, base=10), y=pop_hist), colour="black") + geom_vline(xintercept = log(LGM_max, base=10), colour = "blue") + geom_vline(xintercept = log(LGM_deg, base=10), colour = "blue") + geom_vline(xintercept = log(IA_start, base = 10), colour = "cyan") + geom_vline(xintercept = log(cutoff7gen, base = 10), colour = "red") + geom_vline(xintercept = log(cutoff10gen, base = 10), colour = "black")



############## BEN'S WORK (modified) ############
gen_ave = 8.5
kanga1$t_unscaled_ave = kanga1$t_k*2 * N_0 * gen_ave
cutoff8.5gen = c(20000/25*8.5, 3000000/25*8.5)

n <- dim(kanga1)[1] # how many time points there are
lgTime <- seq(LGM_deg,IA_start,length.out=n) # time between ice age start and antartic ice sheets declining
lgmupper <- rep(max(kanga1$pop_hist)*1.1,n) # lgm upper/lower are just things to fill in
lgmlower <- rep(0,n)

g_hi <- ggplot(kanga1,aes(x=t_unscaled_hi,y=pop_hist))+theme_bw()+geom_step()+ggtitle('Upper Estimate')+coord_cartesian()+
  xlab("\nYears Before Present")+ylab(bquote(''*N[e]*'\n'))+geom_vline(xintercept=c(cutoff10gen[1]),colour='red',linetype='dashed')+
  geom_ribbon(aes(x=lgTime,ymax=lgmupper,ymin=lgmlower),colour='blue',alpha=0.4,fill='blue') # alpha is transparency
ggsave(file="/home/alex/Desktop/Data/redKangaroo/rPlots/redKangaroo_10gen.pdf",g_hi)

g_ave <- ggplot(kanga1,aes(x=t_unscaled_ave,y=pop_hist))+theme_bw()+geom_step()+ggtitle('Average Estimate')+coord_cartesian(ylim=c(0,max(kanga1$pop_hist)*1.1))+
  xlab("\nYears Before Present")+ylab(bquote(''*N[e]*'\n'))+geom_vline(xintercept=c(cutoff8.5gen[1]),colour='red',linetype='dashed')+
  geom_ribbon(aes(x=lgTime,ymax=lgmupper,ymin=lgmlower),colour='blue',alpha=0.4,fill='blue') # alpha is transparency
ggsave(file="/home/alex/Desktop/Data/redKangaroo/rPlots/redKangaroo_8.5gen.pdf",g_ave)

g_low <- ggplot(kanga1,aes(x=t_unscaled_low,y=pop_hist))+theme_bw()+geom_step()+ggtitle('Lower Estimate')+coord_cartesian(ylim=c(0,max(kanga1$pop_hist)*1.1))+
  xlab("\nYears Before Present")+ylab(bquote(''*N[e]*'\n'))+geom_vline(xintercept=c(cutoff7gen[1]),colour='red',linetype='dashed')+
  geom_ribbon(aes(x=lgTime,ymax=lgmupper,ymin=lgmlower),colour='blue',alpha=0.4,fill='blue') # alpha is transparency
ggsave(file="/home/alex/Desktop/Data/redKangaroo/rPlots/redKangaroo_7gen.pdf",g_low)

######################################

#############################3
################## CRAP BELOW...

# ### Can we work out when back in time this was? NOTE: STRONG POSSIBILITY THIS IS COMPLETELY WRONG (update: yes it is....)
# 
# # According to the Department of the Environment https://www.environment.gov.au/biodiversity/wildlife-trade/natives/wild-harvest/kangaroo-wallaby-statistics/kangaroo-population (accessed 3/2/15), the 2011 population estimates for kangaroos within the commercial harvest areas was 11514298.
# # Would it be reasonable to divide by 2 to get effective population size? I have no idea.
# N_0 = ceiling(11514298/2) # I haven't found an upper and lower estimate at the moment.
# # Generation time: 7-10 years (Dawson 2012) http://www.environment.nsw.gov.au/resources/threatenedspecies/determinations/PDRedKangReject.pdf
# gen_low = 7
# gen_hi = 10
# 
# # let's unscale time to years. At the moment, 1 unit of scaled time is equal to 2N_0 generations having passed.
# kanga1$t_unscaled_low = kanga1$t_k*2 * N_0 * gen_low
# kanga1$t_unscaled_hi = kanga1$t_k*2 * N_0 * gen_hi
# kanga1$pop_hist = kanga1$lambda_k * N_0
# 
# # when was the LGM... according to wikipedia lol,
# LGM_start = 26500
# LGM_end = 19000
# 
# # let's see if the population dynamics can be matched with the LGM? LOG VER
# ggplot(data=kanga1, mapping=aes(x = log(t_unscaled_low, base=10), y = pop_hist)) + geom_step(colour="red") + labs(xlab("log10(years in past)"), ylab("population size")) + ggtitle("The red kangaroo's population dynamics, and the LGM") + geom_step(data = kanga1, mapping=aes(x=log(t_unscaled_hi, base=10), y=pop_hist), colour="black") + geom_vline(xintercept = log(LGM_end, base=10)) + geom_vline(xintercept = log(LGM_start, base=10))
# # let's see if the population dynamics can be matched with the LGM? UNLOGGED VER
# #ggplot(data=kanga1, mapping=aes(x = t_unscaled_low, y = pop_hist)) + geom_step(colour="red") + labs(xlab("years in past"), ylab("population size")) + ggtitle("The red kangaroo's population dynamics, and the LGM") + geom_step(data = kanga1, mapping=aes(x=t_unscaled_hi, y=pop_hist), colour="black") + geom_vline(xintercept = LGM_end) + geom_vline(xintercept = LGM_start)
