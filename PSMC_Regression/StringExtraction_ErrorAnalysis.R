library(stringr)

setwd("~/Documents/SummerScholarship/simulatedData/regressionSimulations/") 
#setwd("/home/alex/Desktop/Simulations/Regression/")

filelist = list.files(pattern = "*.txt",recursive=T)

############### TRUE POPULATION AND TIME VALUES AS SPECIFIED BY THE MS COMMANDS
ms_eN_2_vec = function(ms_input, t_or_p) {
# takes ms_input of form "-eN 0.0055 0.0832 -eN 0.0089 0.0489" copied from ms
# t_or_p can be TRUE or FALSE. If TRUE, gives times. If FALSE, gives populations.
  TimesPops = as.numeric(unlist(str_extract_all(gsub(pattern = "-eN ", replacement = "", x = ms_input), "[\\.0-9e-]+")))  # remove the -eN switch, convert to a vector
  if (t_or_p == TRUE) {
    Times = 4*TimesPops[c(TRUE,FALSE)] # extracts the odd values (i.e. the times) and scales
    return(Times)
  } # returning pops
  else {
    Pops = 4*TimesPops[c(FALSE,TRUE)] # extracts the even values (i.e. the pops) and scales
    return(Pops)    
  }
}
############### TRUE POPULATION AND TIME VALUES AS SPECIFIED BY THE MS COMMANDS extracted from a *.ms file
msFile_eN_2_vec = function(msFilePath, t_or_p) {
  # takes msFilePath as a path to a .ms file and returns:
  # t_or_p can be TRUE or FALSE. If TRUE, gives times. If FALSE, gives populations.
  ms_input=system2('head',args=paste('-n 1 ',msFilePath),stdout=TRUE)
  ms_input=paste('-eN',regmatches(ms_input, regexpr("-eN", ms_input), invert = TRUE)[[1]][2])
  ms_input=regmatches(ms_input, regexpr("-l", ms_input), invert = TRUE)[[1]][1]
  TimesPops = as.numeric(unlist(str_extract_all(gsub(pattern = "-eN ", replacement = "", x = ms_input), "[\\.0-9e-]+")))  # remove the -eN switch, convert to a vector
  if (t_or_p == TRUE) {
    Times = 4*TimesPops[c(TRUE,FALSE)] # extracts the odd values (i.e. the times) and scales
    return(Times)
  } # returning pops
  else {
    Pops = 4*TimesPops[c(FALSE,TRUE)] # extracts the even values (i.e. the pops) and scales
    return(Pops)    
  }
}
#FauxHumanMs_eN = "-eN 0.0055 0.0832 -eN 0.0089 0.0489 -eN 0.0130 0.0607 -eN 0.0177 0.1072 -eN 0.0233 0.2093 -eN 0.0299 0.3630 -eN 0.0375 0.5041 -eN 0.0465 0.5870 -eN 0.0571 0.6343 -eN 0.0695 0.6138 -eN 0.0840 0.5292 -eN 0.1010 0.4409 -eN 0.1210 0.3749 -eN 0.1444 0.3313 -eN 0.1718 0.3066 -eN 0.2040 0.2952 -eN 0.2418 0.2915 -eN 0.2860 0.2950 -eN 0.3379 0.3103 -eN 0.3988 0.3458 -eN 0.4701 0.4109 -eN 0.5538 0.5048 -eN 0.6520 0.5996 -eN 0.7671 0.6440 -eN 0.9020 0.6178 -eN 1.0603 0.5345 -eN 1.4635 1.7931" # as specified by ms
#FauxHumanTimesPops = as.numeric(unlist(str_extract_all(gsub(pattern = "-eN ", replacement = "", x = FauxHumanMs_eN), "[\\.0-9e-]+")))  # remove the -eN switch, convert to a vector
#FauxHumanTimes = 4*FauxHumanTimesPops[c(TRUE,FALSE)] # extracts the odd values (i.e. the times) and scales
#FauxHumanTimes

#FauxHumanPops = 4*FauxHumanTimesPops[c(FALSE, TRUE)] # extracts the even values (i.e. the pops) and scales

# need to define true times and populations
#msFilelist = list.files(pattern = "*ms",recursive=T)
#for(msFile in msFilelist){
#  
#}

fauxHumanMs_eN = "-eN 0.0055 0.0832 -eN 0.0089 0.0489 -eN 0.0130 0.0607 -eN 0.0177 0.1072 -eN 0.0233 0.2093 -eN 0.0299 0.3630 -eN 0.0375 0.5041 -eN 0.0465 0.5870 -eN 0.0571 0.6343 -eN 0.0695 0.6138 -eN 0.0840 0.5292 -eN 0.1010 0.4409 -eN 0.1210 0.3749 -eN 0.1444 0.3313 -eN 0.1718 0.3066 -eN 0.2040 0.2952 -eN 0.2418 0.2915 -eN 0.2860 0.2950 -eN 0.3379 0.3103 -eN 0.3988 0.3458 -eN 0.4701 0.4109 -eN 0.5538 0.5048 -eN 0.6520 0.5996 -eN 0.7671 0.6440 -eN 0.9020 0.6178 -eN 1.0603 0.5345 -eN 1.4635 1.7931"
fauxHumanTimes = ms_eN_2_vec(fauxHumanMs_eN,TRUE)
fauxHumanPops = ms_eN_2_vec(fauxHumanMs_eN,FALSE)
constantPopTimes = 4*c(0.01,1:100) # dummy times
constantPopPops = 4*rep(1,length(constantPopTimes))
psmcSim1MS_eN = "-eN 0.01 0.1 -eN 0.06 1 -eN 0.2 0.5 -eN 1 1 -eN 2 2"
psmcSim1Times = ms_eN_2_vec(psmcSim1MS_eN,TRUE)
psmcSim1Pops = ms_eN_2_vec(psmcSim1MS_eN,FALSE)
psmcSim2MS_eN = "-eN 0.1 1 -eN 0.6 4 -eN 1 1 -eN 3 2 -eN 5 1"
psmcSim2Times = ms_eN_2_vec(psmcSim2MS_eN,TRUE)
psmcSim2Pops = ms_eN_2_vec(psmcSim2MS_eN,FALSE)
trench_Ms_eN = "-eN 0.01 1 -eN 0.1 0.8 -eN 0.2 0.6 -eN 0.3 0.4 -eN 0.4 0.3 -eN 0.5 0.25 -eN 0.6 0.2 -eN 0.7 0.3 -eN 0.8 0.4 -eN 0.9 0.6 -eN 1 0.8 -eN 1.1 1"
trenchTimes = ms_eN_2_vec(trench_Ms_eN,TRUE)
trenchPops = ms_eN_2_vec(trench_Ms_eN,FALSE)

# note: we can set these here because it's the same for all of the simulations we've run so far. HOWEVER, these rates are scaled by the number of bp in the sample. this will be addressed in the for loop later
mut_rate_theta = 65130
recomb_rate_rho = 10973

##############################33
# need to do a min of maxs, max of mins to determine integration limits
logmins = rep(0,length(filelist)+1) # the plus one for the smallest (and largest) true times
logmaxs = rep(0,length(filelist)+1)
i = 1
for (infile in filelist) {
  data.infile = read.table(infile,header=TRUE)
  logmins[i] = log(min(data.infile$t_k/2)) # gotta divide by 2 to scale PSMC time output to N_0 generations. THIS IS IMPORTANT OTHERWISE THE PROGRAMS STUFF! MUST DO WHENEVER THERE IS A t_k
  logmaxs[i] = log(max(data.infile$t_k/2))
  i = i+1
}
logmins[length(filelist)+1] = log(max(min(fauxHumanTimes),min(constantPopTimes),min(psmcSim2Times),min(psmcSim1Times),min(trenchTimes))) # add true value min
logmaxs[length(filelist)+1] = log(min(max(fauxHumanTimes),max(constantPopTimes),max(psmcSim2Times),max(psmcSim1Times),max(trenchTimes))) # add true value max
logmax_of_mins = max(logmins)
logmin_of_maxs = min(logmaxs)

# need to set up integration functions
# this function will tell you PSMC's population size estimate at time "pos" given the x (time) and y (relative pop) data
eval_popsize = function(pos, x, y){ 
  xIndex = length(which(x <= pos))
  if(xIndex==0){
    print(pos)
    print(x)
    print(y)
    return(NaN)
  }
  return(y[xIndex])
}

# this function will tell you the absolute difference in between two curves A and B at position "xpos". we will later integrate this to get an area between two curves
absdifference = function(xpos,d1,d2,d3,d4){ # d1, d2 are x and y of curve A. d1, d2 are x and y of curve B. you want to compare curves A and B.
  out = NULL
  n=length(xpos)
  for (i in 1:n) {
    out[i] = abs( eval_popsize(xpos[i],d1,d2) - eval_popsize(xpos[i],d3,d4) )
  }
  return(out)
}

# EXTRACT NUMBERS BETWEEN SPECIFIC WORDS! reference: http://stackoverflow.com/questions/22924183/extract-symbols-between-two-particular-words
#numbaaaa = gsub(pattern = "(.*Split)(.*)(\\.psmcfa.*)",
#     replacement = "\\2",
#     x = filelist[1])
#numbaaaa = as.numeric(numbaaaa)
#numbaaaa

#lettaaaa = gsub(pattern = "(.*)(Mbp.*)",
#                replacement = "\\1",
#                x = filelist[1])
#lettaaaa

#### END TEST

# Initialise input variables
Bp = NULL
Int = NULL
Split = NULL
Bp_per_contig = NULL
Per_Site_Mut_Rate = NULL
Per_Site_Recomb_Rate = NULL
Pop_Dynamics_Type = NULL
Error = NULL

# files have names like "fauxHuman/fauxHumanBp100000Int10*1Split1/fauxHumanBp100000Int10*1Split1.txt"

for (infile in filelist) { # make sure I change things to Shaun's file format
  Bp_tmp = as.numeric(gsub(pattern = "(.*Bp)(.*)(Int.*)", replacement = "\\2", x = infile))
  Bp = append(Bp,Bp_tmp)
  Int = append(Int,as.numeric(gsub(pattern = "(.*Int)(.*)(\\*)(.*)(Split.*)", replacement = "\\2", x = infile))) # for something like 10*1, only extracts the 10 (we wouldn't do anything like x*y for y!= 1)
  Split_tmp = as.numeric(gsub(pattern = "(.*Split)(.*)(\\.txt)", replacement = "\\2", x = infile))
  Split = append(Split,Split_tmp)
  Bp_per_contig = append(Bp_per_contig,ceiling(Bp_tmp/Split_tmp)) # to convert from Mbp to bp)
  Pop_Dynamics_Type_tmp = gsub(pattern = "(\\/.*)(.*)(Bp.*)", replacement = "\\2", x = infile)
  Pop_Dynamics_Type = append(Pop_Dynamics_Type,Pop_Dynamics_Type_tmp)
  
  Per_Site_Mut_Rate_tmp = mut_rate_theta/Bp_tmp
  Per_Site_Mut_Rate = append(Per_Site_Mut_Rate,Per_Site_Mut_Rate_tmp)
  Per_Site_Recomb_Rate_tmp = recomb_rate_rho/Bp_tmp
  Per_Site_Recomb_Rate = append(Per_Site_Recomb_Rate,Per_Site_Mut_Rate_tmp)
  
  # error analysis
  data = read.table(infile,header=TRUE)
  #simplified if statements to concatenating Pop_Dynamics_Type_tmp (which should hopefully contain the name of the sim we're running)
  ErrorVal = integrate(absdifference, logmax_of_mins, logmin_of_maxs, d1=log(get(paste(Pop_Dynamics_Type_tmp,'Times',sep=''))), d2=get(paste(Pop_Dynamics_Type_tmp,'Pops',sep='')), d3=log(data$t_k/2), d4=data$lambda_k, subdivisions=10000)$value 
  Error = append(Error,ErrorVal)
} 
Results = data.frame(Bp,Int,Split,Bp_per_contig,Per_Site_Mut_Rate,Per_Site_Recomb_Rate,Pop_Dynamics_Type,Error)

Results


######################################################
######################################################
######################################################

# Stuff broken ... time to fix

#plot(log(fauxHumanTimes),fauxHumanPops,"s",ylim=c(0,17),xlim=c(-3,4),lwd=5,main=">_<")
#collist = c("blue","brown","gold","red","green","purple","cyan","orange","magenta","pink","darkolivegreen")

#counter = 1
#for (infile in rev(filelist)) {
#  infile.data = read.table(infile,header=TRUE)
#  lines(log(infile.data$t_k/2),infile.data$lambda_k,"s",col=collist[counter],lwd=2.5)
#  counter = counter + 1
#}
#legend("topleft",rev(filelist),fill=collist,bty="n",cex=0.8) #

