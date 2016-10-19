#!/bin/bash
# Collects all the AUPR values from the results for a cell line.
# You need to have run the AUCCalculator already (run_aucpr_crossval/crossline.sh).

if [ $# != 2 ]; then
	echo "USAGE: collect_aupr.sh cellname results_top"
	echo "cellname is the name of the training set for this instance (must be present in first col of datalist)"
	echo "results_top is the top-level of an output directory that we already made to store the training/testing folds. (Not unique per cell type.)"
	exit
fi

trainline=$1
restop=$2

if [ ! -d $restop ]; then
	echo "Cannot find output directory $restop"
	exit 1
fi

# locate directories
traindir=${restop}/train_${trainline}/
if [ ! -d $traindir ];
then
	echo "Cannot find $traindir"
	exit 1
fi

# make a master output file
masterfile=${traindir}/master_auc_${trainline}.txt

printf "trainline\ttestdir\tntrees\taupr\n" > $masterfile
for dir in ${traindir}/result_cv ${traindir}/test_*;
do
	dirbase=$(basename $dir)
	if [ "$dirbase" == "result_cv" ]; then
		# name for test in the output file
		dirbase="test_${trainline}_CV"
	fi
	for ntrees in 50 #1 10 50 100 500
	do
		treedir=${dir}/trees_${ntrees}
		if [ ! -d $treedir ]; then
			continue
		fi
		# get AUC file
		aucfile=${treedir}/train_${trainline}*.auc
		if [ ! -e $aucfile ]; then
			continue
		fi
		#echo $treedir
		aucpr=$(tail -n 2 $aucfile | head -n 1  | cut -f10 -d" ")
		printf "%s\t%s\t%d\t%f\n" $trainline $dirbase $ntrees $aucpr >> $masterfile

	done
done

echo $masterfile
#head $masterfile

