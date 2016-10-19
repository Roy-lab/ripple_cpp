#!/bin/bash
# Gather cross-validation predictions for one cell line.
# USAGE: gather_predictions_crossval.sh datalist train_cellline results_top

if [ $# != 2 ]; then
	echo "USAGE: gather_predictions_crossval.sh cellname results_top"
	echo "Runs cross-validation on one cell line."
	echo "Assumes training/testing folds have already been generated under results_top/${cellname}/traindata."
	echo "cellname is the name of the training set for this instance."
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
randomForest=execs/randomForest/randomForest #/mnt/ws/sysbio/roygroup/shared/programs/regForest/randomForest/src/randomForest

if [ ! -e data/ripple_predictor_network.txt ]; then
	echo "Cannot find predictor specification file, data/ripple_predictor_network.txt"
	exit
fi

# locate training/testing folds
traindata=${restop}/train_${trainline}/traindata
echo $restop
if [ ! -d $traindata ]; then
	echo "Cannot find training fold location $traindata"
	exit 1
fi

# do internal CV
cvresult=${restop}/train_${trainline}/result_cv
if [ ! -e $cvresult ];
then
	mkdir $cvresult
fi

echo $cvresult
for ntrees in 50
do
	echo $ntrees
	treesult=${cvresult}/trees_${ntrees}
	if [ -e $treesult ]; then
		echo "We already started this -- skipping"
		continue
	fi

	if [ ! -d $treesult ]; then
		mkdir $treesult
	fi
	for i in {0..9};
	do
		folddir=${treesult}/fold_${i}
		echo $folddir
		if [ -d $folddir ]; then
			rm -r $folddir
		fi
		mkdir $folddir
		trainfold=${traindata}/trainset${i}.txt
		testfold=${traindata}/testset${i}.txt
		$randomForest -v c -m $trainfold -d $testfold -n${ntrees} -l10 -b data/ripple_predictor_network.txt -o ${folddir}  1> ${folddir}/fold_${i}.log 2> ${folddir}/fold_${i}.err
		#exit
	done
done
exit

