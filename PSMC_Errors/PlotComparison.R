# THE TRUTH!
TrueTimes = 4*c(0.01,0.06,0.2,1,2) # from ms command
TruePops = 4*c(0.1,1,0.5,1,2) # dunno why you need to to this *4

# the 30 Mbp test
setwd("/home/alex/Desktop/Simulations/SplittingChromosomes/binarySplitData1Chromosome30Mbp/extractPSMCoutput/")
filelist30 = list.files(pattern = "*.txt")
# filelist30 = filelist30[-3] # remove the fullData file

# fullData30 = read.table("fullData.psmcfa.psmc.txt",header=TRUE)
#plot(log(TrueTimes),TruePops,"s",ylim=c(0,15),xlim=c(-4,3),lwd=4)
#legend("right",c("Black: full"),bty=n)

#for (infile in filelist30) {
#  infile.data = read.table(infile,header=TRUE)
#  lines(log(infile.data$t_k),infile.data$lambda_k,"s",col="red")
#}

# HalfSplitData30 = read.table(filelist30[3],header=TRUE)
# lines(log(HalfSplitData30$t_k),HalfSplitData30$lambda_k,"s",col="green",lwd=2)

# START HERE 
plot(log(TrueTimes),TruePops,"s",ylim=c(0,15),xlim=c(-4,3),lwd=5,main="30Mbp")
collist = c("blue","brown","gold","red","green","purple","cyan","orange","magenta","pink","darkolivegreen")
counter = 1
for (infile in rev(filelist30)) { # the "rev" (reverse order) is just so that the plots are a bit more legible
  infile.data = read.table(infile,header=TRUE)
  lines(log(infile.data$t_k/2),infile.data$lambda_k,"s",col=collist[counter],lwd=2.5)
  counter = counter + 1
}
legend("topleft",rev(filelist30),fill=collist,bty="n",cex=0.8) # so basically, 64 and below are fine and 128 and above are kinda stuffed.


# the 20 Mbp test
setwd("/home/alex/Desktop/Simulations/SplittingChromosomes/binarySplitData1Chromosome20Mbp/")
filelist20 = list.files(pattern = "*.txt")
# filelist20 = filelist20[-4] # remove the fullData file

#fullData20 = read.table("20Mbp_run.ms.psmcfa.psmc.txt",header=TRUE)
plot(log(TrueTimes),TruePops,"s",ylim=c(0,40),xlim=c(-3,4),lwd=5,main="20Mbp")

counter = 1
for (infile in rev(filelist20)) {
  infile.data = read.table(infile,header=TRUE)
  lines(log(infile.data$t_k/2),infile.data$lambda_k,"s",col=collist[counter],lwd=2.5)
  counter = counter + 1
}
legend("topleft",rev(filelist20),fill=collist,bty="n",cex=0.8) #

# the 40 Mbp test
setwd("/home/alex/Desktop/Simulations/SplittingChromosomes/binarySplitData1Chromosome40Mbp/")
filelist40 = list.files(pattern = "*.txt")
# filelist40 = filelist40[-7] # remove the fullData file

#fullData40 = read.table("40Mbp_run.ms.psmcfa.psmc.txt",header=TRUE)
plot(log(TrueTimes),TruePops,"s",ylim=c(0,17),xlim=c(-3,4),lwd=5,main="40Mbp")

counter = 1
for (infile in rev(filelist40)) {
  infile.data = read.table(infile,header=TRUE)
  lines(log(infile.data$t_k/2),infile.data$lambda_k,"s",col=collist[counter],lwd=2.5)
  counter = counter + 1
}
legend("topleft",rev(filelist40),fill=collist,bty="n",cex=0.8) #

# the 40 Mbp test (SECOND TIME)
setwd("/home/alex/Desktop/Simulations/SplittingChromosomes/binarySplitData1Chromosome40Mbp2/")
filelist40_2 = list.files(pattern = "*.txt")
# filelist40 = filelist40[-7] # remove the fullData file

#fullData40 = read.table("40Mbp_run.ms.psmcfa.psmc.txt",header=TRUE)
plot(log(TrueTimes),TruePops,"s",ylim=c(0,7),xlim=c(-3,3),lwd=5,main="40Mbp (2)")

counter = 1
for (infile in rev(filelist40_2)) {
  infile.data = read.table(infile,header=TRUE)
  lines(log(infile.data$t_k/2),infile.data$lambda_k,"s",col=collist[counter],lwd=2.5)
  counter = counter + 1
}
legend("topleft",rev(filelist40_2),fill=collist,bty="n",cex=0.8) #

##### test with some of shaun's stuff

