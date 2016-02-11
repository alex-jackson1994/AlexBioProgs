library(ggplot2)

# the redKangarooInt4+5*3+4.psmc run
taz <- read.delim("/media/alex/My Passport/bioinformaticsScholarship/tasDevilData/tasmanianDevilInt4+5*3+4.txt")

theta_0 = 0.002679 # from the file tasmanianDevilInt4+5*3+4.psmc. We take theta from the final iteration.
mu = 2.5*10^(-8)
s = 100
N_0 = theta_0 / (4*mu) / s # 

# Generation time: 7-10 years (Dawson 2012) http://www.environment.nsw.gov.au/resources/threatenedspecies/determinations/PDRedKangReject.pdf
gen_low = 2
gen_hi = 3.5

taz$t_unscaled_low = taz$t_k*2 * N_0 * gen_low
taz$t_unscaled_hi = taz$t_k*2 * N_0 * gen_hi
taz$pop_hist = taz$lambda_k * N_0

# when was the LGM... according to wikipedia lol,
LGM_max = 26500 # time for maximum pos the ice sheets reached
LGM_deg = 14500 # deglaciation
IA_start = 110000 # Ice Age, not LGM... from https://en.wikipedia.org/wiki/Last_glacial_period

Humans = 50000

cutoff2gen = c(20000/25*2, 3000000/25*2) # cutoffs of PSMC "goodness" based upon Heng Li's 20 kyr and 3 Myr estimates, but scaled to 2 year generations instead of 25 year
cutoff3.5gen = c(20000/25*3.5, 3000000/25*3.5)

n <- dim(taz)[1] # how many time points there are
lgTime <- seq(LGM_deg,IA_start,length.out=n) # time between ice age start and antartic ice sheets declining
lgmupper <- rep(max(taz$pop_hist)*1.1,n) # lgm upper/lower are just things to fill in
lgmlower <- rep(0,n)

# high estimate
 g_hi_low <- ggplot(taz,aes(x=t_unscaled_hi,y=pop_hist))+theme_bw()+geom_step(colour="blue")+ ggtitle('Upper And Lower Estimates')+coord_cartesian(xlim=c(0,1.3*10^5))+geom_point(colour="blue")+
  xlab("\nYears Before Present")+ylab(bquote(''*N[e]*'\n'))+geom_vline(xintercept=c(cutoff10gen[1]),colour='red', linetype="dashed")+
  geom_ribbon(aes(x=lgTime,ymax=lgmupper,ymin=lgmlower),colour='blue',alpha=0.4,fill='blue') +# alpha is transparency
   # low estimate
   geom_step(data=taz, aes(x=t_unscaled_low, y=pop_hist ), colour="green")+geom_point(data=taz, aes(x=t_unscaled_low, y=pop_hist ), colour="green")
   # mark the PSMC time intervals
#      geom_vline( xintercept=psmc_int_ts_low, colour="green", linetype="dashed" ) + geom_vline( xintercept=psmc_int_ts_hi, colour="blue", linetype="dashed" )

 g_hi_low

 
 
 
#ggsave(file="/home/alex/Desktop/Data/redKangaroo/rPlots/redKangaroo_bothgens.pdf",g_hi_low)

############# PLOT ALL
setwd("/home/alex/Desktop/Data/tasmanianDevil")
filelist = list.files(pattern = "*.txt")

# We need a script to grab the theta_0 values but I'm not good enough to do that.
theta_0_s  = c(0.003810, 0.003654, 0.003574, 0.003519, 0.003484, 0.003462, 0.003445, 0.002679) # taking these manually, which is a pain in the butt...
N_0_s = theta_0_s / (4*mu) / s # = 1854 apparently...

count = 1
for (file in filelist) {
  dat <- read.delim(file)
  
  dat$t_unscaled_low = dat$t_k*2 * N_0_s[count] * gen_low
  dat$t_unscaled_hi = dat$t_k*2 * N_0_s[count] * gen_hi
  dat$pop_hist = dat$lambda_k * N_0_s[count]
  dat$type = file
  
  n <- dim(dat)[1] # how many time points there are
  lgTime <- seq(LGM_deg,IA_start,length.out=n) # time between ice age start and antartic ice sheets declining
  lgmupper <- rep(max(dat$pop_hist)*1.1,n) # lgm upper/lower are just things to fill in
  lgmlower <- rep(0,n)
  
  dat_plot <- ggplot(dat,aes(x=t_unscaled_hi,y=pop_hist))+theme_bw()+geom_step(colour="red")+ ggtitle('Upper Estimate (red) and Lower Estimate (blue)')+coord_cartesian(xlim=c(0,1*10^5),ylim=c(0,30000))+geom_point(colour="red")+
    xlab("\nYears Before Present")+ylab(bquote(''*N[e]*'\n'))+geom_vline(xintercept=c(cutoff10gen[1]),colour='red', linetype="dashed")+
    geom_ribbon(aes(x=lgTime,ymax=lgmupper,ymin=lgmlower),colour='purple',alpha=0.3,fill='purple') +
   geom_step(data=dat, aes(x=t_unscaled_low, y=pop_hist ), colour="blue")+geom_point(data=dat, aes(x=t_unscaled_low, y=pop_hist ), colour="blue")+geom_vline(xintercept=c(cutoff7gen[1]),colour='blue', linetype="dashed")+
    geom_vline(xintercept=Humans,colour='black', linetype="dashed")
  
  #setwd("/home/alex/Desktop/Data/redKangaroo/rPlots/Ints/")
  print(dat$pop_hist)
  print(file)
  ggsave(paste("plots/",file, ".pdf", sep=""), dat_plot)
#   print(file)
#   print(dat)
  #alloutput = rbind(alloutput, dat)
  
  count = count + 1
}

#####
# DEVIL + KANGA ON SAME PLOT
rk_dat <- read.delim("/home/alex/Desktop/Data/redKangaroo/Ints/redKangaroo/redKangarooBp1020215700Int40*1.txt")
td_dat <- read.delim("/home/alex/Desktop/Data/tasmanianDevil/tasmanianDevilBp3053144200Int40*1.txt")

rk_theta0 = 0.022849
rk_N0 = rk_theta0 / (4*mu) / s
rk_gen_low = 7
rk_gen_hi = 10
cutoff7gen = 20000/25*7
cutoff10gen = 20000/25*10

td_theta0 = 0.003445
td_N0 = td_theta0 / (4*mu) / s
td_gen_low = 2
td_gen_hi = 3.5
cutoff2gen = 20000/25*2
cutoff3.5gen = 20000/25*3.5

rk_dat$t_unscaled_low = rk_dat$t_k*2 * rk_N0 * rk_gen_low
rk_dat$t_unscaled_hi = rk_dat$t_k*2 * rk_N0 * rk_gen_hi
rk_dat$pop_hist = rk_dat$lambda_k * rk_N0

td_dat$t_unscaled_low = td_dat$t_k*2 * td_N0 * td_gen_low
td_dat$t_unscaled_hi = td_dat$t_k*2 * td_N0 * td_gen_hi
td_dat$pop_hist = td_dat$lambda_k * td_N0
############## ACAD PLOTS

stacktimes = c(rk_dat$t_unscaled_low, rk_dat$t_unscaled_hi, td_dat$t_unscaled_low, td_dat$t_unscaled_hi)
stackpops = c(rk_dat$pop_hist, rk_dat$pop_hist, td_dat$pop_hist, td_dat$pop_hist)
types = c(rep("RK7Gen", length(rk_dat$t_unscaled_low)), rep("RK10Gen", length(rk_dat$t_unscaled_hi)), rep("TD2Gen", length(td_dat$t_unscaled_low)), rep("TD3.5Gen", length(td_dat$t_unscaled_hi)))
cutoffs = c(rep(cutoff7gen, length(rk_dat$t_unscaled_low)), rep(cutoff10gen, length(rk_dat$t_unscaled_hi)), rep(cutoff2gen, length(td_dat$t_unscaled_low)), rep(cutoff3.5gen, length(td_dat$t_unscaled_hi)))
kan_tas = data.frame(stacktimes, stackpops, types, cutoffs)
colnames(kan_tas) = c("Times","Pops","Type","Cutoffs")

### KANGA + DEVIL

#LGM_max = 26500 # time for maximum pos the ice sheets reached
LGM_deg = 14500 # deglaciation
IA_start = 110000 # Ice Age, not LGM... from https://en.wikipedia.org/wiki/Last_glacial_period
n <- dim(kan_tas)[1] # how many time points there are
lgTime <- seq(LGM_deg,IA_start,length.out=n) # time between ice age start and antartic ice sheets declining
lgmupper <- rep(max(kan_tas$Pops)*1.1,n) # lgm upper/lower are just things to fill in
lgmlower <- rep(0,n)
kan_tas_plot = ggplot(kan_tas, aes(x=Times, y=Pops, colour = Type ))+theme_bw()+geom_step(size=1) +
  xlab("\nYears Before Present")+ylab(bquote(''*N[e]*'\n'))+ ggtitle('Red Kangaroo and Tasmanian Devil Effective Population Estimates from PSMC') + coord_cartesian(xlim=c(20,4*10^5),ylim=c(0,10000))+
     geom_ribbon(aes(x=lgTime,ymax=lgmupper,ymin=lgmlower),colour='yellow',alpha=0.2,fill='yellow') +
  geom_vline(xintercept=Humans, colour='black') + geom_vline(aes(xintercept=cutoffs, colour=Type), linetype="dashed")

kan_tas_plot
ggsave("/home/alex/Desktop/Data/tasmanianDevil/plots/ACAD/kan_tas_nocut_nolog.pdf", kan_tas_plot, paper="a4r", width=27, height=19)

kan_tas_plot_log = kan_tas_plot + ggtitle('Red Kangaroo and Tasmanian Devil Effective Population Estimates from PSMC (Log Time)') + scale_x_continuous(trans="log10") + coord_cartesian(xlim=c(500,4*10^5),ylim=c(0,10000))
ggsave("/home/alex/Desktop/Data/tasmanianDevil/plots/ACAD/kan_tas_nocut_log.pdf", kan_tas_plot_log, paper="a4r", width=27, height=19)

### TAKING CUTOFFS INTO ACCOUNT
# types = c(rep("RK7Gen", length(rk_dat$t_unscaled_low)), rep("RK10Gen", length(rk_dat$t_unscaled_hi)), rep("TD2Gen", length(td_dat$t_unscaled_low)), rep("TD3.5Gen", length(td_dat$t_unscaled_hi)))
kan_tas_cut = kan_tas[!(kan_tas$Type=="RK7Gen" & kan_tas$Times<cutoff7gen),] # remove all lines which are RK7Gen and have a time less than the cutoff fpr 7, etc.
kan_tas_cut = kan_tas_cut[!(kan_tas_cut$Type=="RK10Gen" & kan_tas_cut$Times<cutoff10gen),]
kan_tas_cut = kan_tas_cut[!(kan_tas_cut$Type=="TD2Gen" & kan_tas_cut$Times<cutoff2gen),]
kan_tas_cut = kan_tas_cut[!(kan_tas_cut$Type=="TD3.5Gen" & kan_tas_cut$Times<cutoff3.5gen),]

n <- dim(kan_tas_cut)[1] # how many time points there are
lgTime <- seq(LGM_deg,IA_start,length.out=n) # time between ice age start and antartic ice sheets declining
lgmupper <- rep(max(kan_tas_cut$Pops)*1.1,n) # lgm upper/lower are just things to fill in
lgmlower <- rep(0,n)
kan_tas_cutoff_plot = ggplot(kan_tas_cut, aes(x=Times, y=Pops, colour = Type ))+theme_bw()+geom_step(size=1) +
  xlab("\nYears Before Present")+ylab(bquote(''*N[e]*'\n'))+ coord_cartesian(xlim=c(20,4*10^5),ylim=c(0,10000))+
  geom_ribbon(aes(x=lgTime,ymax=lgmupper,ymin=lgmlower),colour='yellow',alpha=0.2,fill='yellow') +
  geom_vline(xintercept=Humans, colour='black') #+ geom_vline(aes(xintercept=cutoffs, colour=Type), linetype="dashed")

kan_tas_cutoff_plot_log = kan_tas_cutoff_plot + ggtitle('Red Kangaroo and Tasmanian Devil With Cutoffs (Log Time)') + scale_x_continuous(trans="log10") + coord_cartesian(xlim=c(500,4*10^5),ylim=c(0,10000))
ggsave("/home/alex/Desktop/Data/tasmanianDevil/plots/ACAD/kan_tas_cut_log.pdf", kan_tas_cutoff_plot_log, paper="a4r", width=27, height=19)


### KANGA HI-RES VS LO-RES
rk_lores_dat = read.delim("/home/alex/Desktop/Data/redKangaroo/redKangarooInt4+5*3+4.txt")

rk_lores_theta0 = 0.018540
rk_lores_N0 = rk_lores_theta0 / (4*mu) / s
rk_lores_dat$t_unscaled_low = rk_lores_dat$t_k*2 * rk_lores_N0 * rk_gen_low
rk_lores_dat$t_unscaled_hi = rk_lores_dat$t_k*2 * rk_lores_N0 * rk_gen_hi
rk_lores_dat$pop_hist = rk_lores_dat$lambda_k * rk_lores_N0

stacktimes = c(rk_dat$t_unscaled_low, rk_dat$t_unscaled_hi, rk_lores_dat$t_unscaled_low, rk_lores_dat$t_unscaled_hi)
stackpops = c(rk_dat$pop_hist, rk_dat$pop_hist, rk_lores_dat$pop_hist, rk_lores_dat$pop_hist)
types = c(rep("RK7GenHiRes", length(rk_dat$t_unscaled_low)), rep("RK10GenHiRes", length(rk_dat$t_unscaled_hi)), rep("RK7GenLoRes", length(rk_lores_dat$t_unscaled_low)), rep("RK10GenLoRes", length(rk_lores_dat$t_unscaled_hi)))
cutoffs = c(rep(cutoff7gen, length(rk_dat$t_unscaled_low)), rep(cutoff10gen, length(rk_dat$t_unscaled_hi)), rep(cutoff7gen, length(rk_lores_dat$t_unscaled_low)), rep(cutoff7gen, length(rk_lores_dat$t_unscaled_hi)))
kanga_res = data.frame(stacktimes, stackpops, types, cutoffs)
colnames(kanga_res) = c("Times","Pops","Type","Cutoffs")

n <- dim(kanga_res)[1] # how many time points there are
lgTime <- seq(LGM_deg,IA_start,length.out=n) # time between ice age start and antartic ice sheets declining
lgmupper <- rep(max(kanga_res$Pops)*1.1,n) # lgm upper/lower are just things to fill in
lgmlower <- rep(0,n)
kanga_res_plot = ggplot(kanga_res, aes(x=Times, y=Pops, colour = Type ))+theme_bw()+geom_step(size=1) +
  xlab("\nYears Before Present")+ylab(bquote(''*N[e]*'\n'))+ ggtitle('Tasmanian Devil Effective Population Estimates from PSMC - Comparison of Resolution') + coord_cartesian(xlim=c(500,4*10^5),ylim=c(0,10000))+
  geom_ribbon(aes(x=lgTime,ymax=lgmupper,ymin=lgmlower),colour='yellow',alpha=0.2,fill='yellow') +
  geom_vline(xintercept=Humans, colour='black') + geom_vline(aes(xintercept=cutoffs, colour=Type), linetype="dashed")+ scale_x_continuous(trans="log10")
kanga_res_plot 
ggsave("/home/alex/Desktop/Data/tasmanianDevil/plots/ACAD/kanga_res_log.pdf", kanga_res_plot, paper="a4r", width=27, height=19)




# # bigass ggplot2 command :)
#   #kanga
# rk_td_plot <- ggplot(rk_dat,aes(x=t_unscaled_hi,y=pop_hist))+theme_bw()+geom_step(colour="red")+ ggtitle('Red Kangaroo Upper/Lower Estimates (red/orange) and Tasmanian Devil Upper/Lower Estimates (blue/green)')+coord_cartesian(xlim=c(20,4*10^5),ylim=c(0,10000))+geom_point(colour="red")+
#   xlab("\nYears Before Present")+ylab(bquote(''*N[e]*'\n'))+geom_vline(xintercept=cutoff10gen,colour='red', linetype="dashed")+
#   geom_ribbon(aes(x=lgTime,ymax=lgmupper,ymin=lgmlower),colour='purple',alpha=0.3,fill='purple') +
#   geom_step(data=rk_dat, aes(x=t_unscaled_low, y=pop_hist ), colour="orange")+geom_point(data=rk_dat, aes(x=t_unscaled_low, y=pop_hist ), colour="orange")+geom_vline(xintercept=cutoff7gen,colour='orange', linetype="dashed") +
#   # devil
#   geom_step(data=td_dat, aes(x=t_unscaled_hi, y=pop_hist ), colour="blue")+geom_point(data=td_dat, aes(x=t_unscaled_hi, y=pop_hist ), colour="blue")+geom_vline(xintercept=cutoff3.5gen,colour='blue', linetype="dashed") +
#   geom_step(data=td_dat, aes(x=t_unscaled_low, y=pop_hist ), colour="green")+geom_point(data=td_dat, aes(x=t_unscaled_low, y=pop_hist ), colour="green")+geom_vline(xintercept=cutoff2gen,colour='green', linetype="dashed")
# 
# rk_td_plot + scale_x_continuous(trans='log10')
# 
# ggsave("/home/alex/Desktop/Data/tasmanianDevil/Devil+Kanga_Int40.pdf", rk_td_plot)
# ggsave("/home/alex/Desktop/Data/tasmanianDevil/Devil+Kanga_Int40Log.pdf", rk_td_plot + scale_x_continuous(trans='log10'))
