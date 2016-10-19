#!/bin/bash
# Get results from cross-validation and run the AUCCalculator.
# starting with one training cell line.
# USAGE: gather_predictions datalist train_cellline results_top

if [ $# != 2 ]; then
	echo "USAGE: run_auccalc.sh cellname results_top"
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

# executables
auc=execs/auc/auc.jar #/mnt/ws/sysbio/roygroup/shared/thirdparty/aupr/auc.jar


# locate internal CV
cvresult=${restop}/train_${trainline}/result_cv
if [ ! -d $cvresult ];
then
	echo "Cannot find $cvresult"
	exit 1
fi

echo $cvresult
for ntrees in 50 #1 10 50 100 500
do
	echo $ntrees
	treesult=${cvresult}/trees_${ntrees}
	echo $treesult

	if [ ! -d $treesult ]; then
		echo "Cannot find $treesult"
		continue
	fi

	if [ ! -d ${treesult}/fold_9 ]; then
		echo "Fold 9 not finished for $treesult"
		continue
	fi

	catfile=${treesult}/train_${trainline}_trees${ntrees}.txt 
	if [ -e $catfile ]; then
		rm $catfile
	fi
	for i in {0..9};
	do
		folddir=${treesult}/fold_${i}
		if [ ! -d $folddir ]; then
			echo "CANNOT FIND $folddir"
			rm $catfile
			break
		fi
		errfile=${folddir}/testset_error.txt
		cat $errfile >> $catfile
	done

	if [ ! -e $catfile ]; then
		continue
	fi
	
	# make list file out of it
	listfile=${catfile/.txt/.list}	
	aucfile=${catfile/.txt/.auc}
	cut -f 5,6 $catfile | grep -v ConfIn > $listfile


	echo "Running AUPR tool on $listfile"
	java -jar $auc $listfile list > ${aucfile}

	tail -n 2 ${aucfile}

done
exit

