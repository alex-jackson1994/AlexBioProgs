library(stringr)

setwd("~/Documents/SummerScholarship/simulatedData/regressionSimulations/") 
#setwd("/home/alex/Desktop/Simulations/Regression/")

filelist = list.files(pattern = "*.txt",recursive=T)

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

#Automatically discover all files and assign times and pop sizes to variables
simulationNames=vector()
simulationDirs=list.dirs(recursive = FALSE)
for(sims in simulationDirs){
  firstMsFile = list.files(pattern = "*ms",path = sims,recursive=T)[1]
  msFilePath=paste(sims,firstMsFile,sep='/')
  simName=strsplit(sims,split='/')[[1]][2] #split between ./ and the the name, [[1]][2] chooses the name from the list
  simulationNames=append(simulationNames,simName)
  assign(paste(simName,'Times',sep=''),msFile_eN_2_vec(msFilePath,TRUE))
  assign(paste(simName,'Pops',sep=''),msFile_eN_2_vec(msFilePath,FALSE))
}
constantPopTimes = 4*c(0.01,1:100) # dummy times
constantPopPops = 4*rep(1,length(constantPopTimes))

# note: we can set these here because it's the same for all of the simulations we've run so far. HOWEVER, these rates are scaled by the number of bp in the sample. this will be addressed in the for loop later
mut_rate_theta = 65130
recomb_rate_rho = 10973

##############################33
# need to do a min of maxs, max of mins to determine integration limits
logmins = rep(0,length(filelist)) # the plus one for the smallest (and largest) true times
logmaxs = rep(0,length(filelist))
i = 1
for (infile in filelist) {
  data.infile = read.table(infile,header=TRUE)
  logmins[i] = log(min(data.infile$t_k/2)) # gotta divide by 2 to scale PSMC time output to N_0 generations. THIS IS IMPORTANT OTHERWISE THE PROGRAMS STUFF! MUST DO WHENEVER THERE IS A t_k
  logmaxs[i] = log(max(data.infile$t_k/2))
  i = i+1
}
for(sims in simulationNames){
  logmins=append(logmins,log(min(get(paste(sims,'Times',sep='')))))
  logmaxs=append(logmaxs,log(max(get(paste(sims,'Times',sep='')))))
}
logmax_of_mins = max(logmins)
logmin_of_maxs = min(logmaxs)

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

View(Results)


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

