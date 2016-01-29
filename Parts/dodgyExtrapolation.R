library(stringr)
library(ggplot2)
setwd("/home/alex/Desktop/Data/shaunWork/steppeBison/partAnalysis/steppeBison/")
summaryFileList = list.files(pattern = "*Summary.txt",recursive=T)
filelist = list.files(pattern = "*.txt",recursive=T)
dataFileList<-setdiff(filelist,summaryFileList)

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

findEndPoints<-function(dataFileList){
  logmins = rep(0,length(dataFileList))
  logmaxs = rep(0,length(dataFileList))
  i = 1
  for (infile in dataFileList) {
    data.infile = read.table(infile,header=TRUE)
    logmins[i] = log(min(data.infile$t_k/2)) # gotta divide by 2 to scale PSMC time output to N_0 generations. THIS IS IMPORTANT OTHERWISE THE PROGRAMS STUFF! MUST DO WHENEVER THERE IS A t_k
    logmaxs[i] = log(max(data.infile$t_k/2))
    i = i+1
  }
  logmax_of_mins = max(logmins)
  logmin_of_maxs = min(logmaxs)
  return(c(logmin_of_maxs,logmax_of_mins))
}
endPoints<-findEndPoints(dataFileList)
logmax_of_mins = endPoints[2]
logmax_of_mins = -1 #THESE ENDPOINTS WERE CHOSEN BY EYE. I plotted the data and chose endpoints that seemed to occur before it all became shit
logmin_of_maxs = endPoints[1]
logmin_of_maxs = 1

halfData<-read.table('steppeBisonProb0.5/steppeBisonProb0.5.txt',header=T)
halfTimeEstimate<-log(halfData$t_k/2)
halfPopulationEstimate<-halfData$lambda_k
results<-data.frame()
for(i in 1:length(dataFileList)){
  dataFilePath<-dataFileList[i]
  summaryFilePath<-summaryFileList[i]
  summaryData<-read.table(summaryFilePath,header=T)
  data<-read.table(dataFilePath,header=T)
  ErrorVal<-integrate(absdifference, logmax_of_mins, logmin_of_maxs, d1=fullTimeEstimate, d2=fullPopulationEstimate, d3=log(data$t_k/2), d4=data$lambda_k, subdivisions=10000)$value 
  results<-rbind(results,cbind(dataFilePath,summaryData,ErrorVal))
}

ggplot(data=results, mapping=aes(x=mean,y=ErrorVal))+geom_point(size=4)+geom_smooth()
ggplot(data=results, mapping=aes(x=totalLength,y=ErrorVal))+geom_point(size=4)+geom_smooth()


#LINEAR MODELLING ON PARTITIONS
error.lm<-lm(ErrorVal~mean+totalLength,data=results[-10,])
summary(error.lm)
#boxcox(error.lm)
