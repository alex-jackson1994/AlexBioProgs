#!/bin/bash
#This is a modified version of the simulation script that works for real data. Make sure that simulationName below is the same as the .fq.gz file you're working with
#and change all the paths to programs from my directories

# Alex update: changed just to run over a range of time intervals.
simulationName='tasmanianDevil'
#fq2psmcfa $simulationName'.fq.gz' > $simulationName'.psmcfa'
python3 /home/alex/Desktop/Data/shaunWork/myScripts/summaryStatisticsPsmcfa.py $simulationName'.psmcfa' > $simulationName'Summary.txt'
bpLength=$(tail -n 1 $simulationName'Summary.txt' | awk '{print $5}')
for timeInterval in '10*1' '15*1' '20*1' '25*1' '30*1' '35*1' '40*1'
do
	oldStringName=$simulationName'Bp'$bpLength'Int'$timeInterval
	echo $oldStringName'start'
	cp $simulationName'.psmcfa' $oldStringName'.psmcfa'
	psmc -p $timeInterval $oldStringName'.psmcfa'>$oldStringName'.psmc'
	bash /home/alex/Desktop/Data/shaunWork/myScripts/removeDataFromPSMC.sh $oldStringName'.psmc'>$oldStringName'.txt'
	python3 /home/alex/Desktop/Data/shaunWork/myScripts/summaryStatisticsPsmcfa.py $oldStringName'.psmcfa' > $oldStringName'Summary.txt'
	echo $oldStringName'end'
#	for split in '2' '4' '8' '16' '32' '64' '128' '256' '512' '1024'
#	do
#		stringName=$simulationName'Bp'$bpLength'Int'$timeInterval'Split'$split
#		echo $stringName'start'
#		python ~/Documents/SummerScholarship/myScripts/binarySplitPsmcfaPrint.py $oldStringName'.psmcfa' > $stringName'.psmcfa'
#		~/Documents/SummerScholarship/psmc-master/psmc -p $timeInterval $stringName'.psmcfa'>$stringName'.psmc'
#		bash ~/Documents/SummerScholarship/myScripts/removeDataFromPSMC.sh $stringName'.psmc'>$stringName'.txt'
#		python ~/Documents/SummerScholarship/myScripts/summaryStatisticsPsmcfa.py $stringName'.psmcfa' > $stringName'Summary.txt'
#		oldStringName=$stringName
#		echo $stringName'end'
#	done
#	for split in '1' '2' '4' '8' '16' '32' '64' '128' '256' '512' '1024'
#	do
	stringName=$simulationName'Bp'$bpLength'Int'$timeInterval
	mkdir $stringName
	mv $stringName'.psmc' $stringName'.psmcfa' $stringName'.txt' $stringName'Summary.txt' $stringName
#	done
done
mkdir $simulationName
mv $simulationName* $simulationName
7z a $simulationName
