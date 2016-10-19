#!/bin/bash
# Make training/testing folds for all datasets listed in a data file.
# USAGE: make_datafolds datalist results_top

if [ $# != 2 ]; then
	echo "USAGE: make_datafolds datalist results_top"
	echo "datalist is a 2-column file listing where to find all input feature files: name, filepath"
	echo "results_top is the top-level of an output directory where we will put cell-type-specific sub-directories."
	exit
fi

datafile=$1
restop=$2

if [ ! -e $datafile ]; then
	echo "Cannot find $datafile"
	return
fi
if [ ! -d $restop ]; then
	echo "Making output directory $restop"
	mkdir -p $restop
fi

# executables
makeRowPartitions=execs/makeRowPartitions/makeRowPartitions #/mnt/ws/sysbio/roygroup/shared/programs/makeRowPartitions/makeRowPartitions

while read trainline trainfile
do
	echo "Trainline $trainline, $trainfile" 
	if [ ! -e $trainfile ]; then
		echo "Cannot find training file $trainfile"
		exit
	fi
	# make training/testing folds
	traindata=${restop}/train_${trainline}/traindata
	if [ ! -e $traindata ]; then
		mkdir -p $traindata
	fi
	$makeRowPartitions $trainfile 10 $traindata foldsRand 1 
	echo "Produced training folds in $traindata."
done < $datafile

exit


