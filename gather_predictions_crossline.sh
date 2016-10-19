#!/bin/bash
# Gather predictions cross-cell-line for ONE test cell
# starting with one training cell line, testing on ONE other
# USAGE: gather_predictions_crossline datalist train_cellline test_cellline esults_top

if [ $# != 4 ]; then
	echo "USAGE: gather_predictions_crossline datalist train_cell test_cell results_top"
	echo "For each of 10 folds of training data for a train_cell type, makes predictions for ALL examples for the test cell type."
	echo "datalist is a 2-column file listing where to find all input feature files: name, filepath"
	echo "train_cell is the name of the training set for this instance (must be present in first col of datalist)"
	echo "test_cell is the testing cell type for this instance (must also be present in the datalist file)"
	echo "results_top is the top-level of an output directory that we already made to store the training/testing folds. (Not unique per cell type.)"
	exit
fi

datafile=$1
trainline=$2
testlineMaster=$3
restop=$4

if [ ! -e $datafile ]; then
	echo "Cannot find $datafile"
	exit 1
fi
if [ ! -d $restop ]; then
	echo "Cannot find output directory $restop"
	exit 1
fi

if [ ! -e data/ripple_predictor_network.txt ]; then
        echo "Cannot find predictor specification file, data/ripple_predictor_network.txt"
        exit
fi


# executables
randomForest=execs/randomForest/randomForest #/mnt/ws/sysbio/roygroup/shared/programs/regForest/randomForest/src/randomForest


# locate training/testing folds
traindata=${restop}/train_${trainline}/traindata
echo $restop
if [ ! -d $traindata ]; then
	echo "Cannot find training fold location $traindata"
	exit 1
fi

#echo "OTHER FILES_____"
#grep -v "^${trainline}" $datafile

while read testline testdata
do
	echo $trainline $testline

	# do cross-line predictions
	result=${restop}/train_${trainline}/test_${testline}
	if [ ! -e $result ];
	then
		mkdir $result
	fi

	echo $result
	for ntrees in 50 # 1 10 50 100 #500
	do
		echo $ntrees
		treesult=${result}/trees_${ntrees}
		if [ -d $treesult ]; then
			echo "We already started this -- skipping"
			continue
		fi

		if [ ! -d $treesult ]; then
			mkdir $treesult
		fi	
		# for each training fold
		for i in {0..9};
		do
			folddir=${treesult}/fold_${i}
			if [ -d $folddir ]; then
				rm -r $folddir
			fi
			mkdir $folddir
			trainfold=${traindata}/trainset${i}.txt
	
			#echo $randomForest -v c -m $trainfold -d $testdata -n${ntrees} -e0.1 -k1 -x1 -t0.1 -l10 -b data/ripple_predictor_network.txt -o ${folddir} 
	
			$randomForest -v c -m $trainfold -d $testdata -n${ntrees} -l10 -b data/ripple_predictor_network.txt -o ${folddir}  1> ${folddir}/fold_${i}.log 2> ${folddir}/fold_${i}.err
			echo "Done with RF"
			echo $folddir	
		done
	done
done < <(grep "^${testlineMaster}"$'\t' $datafile)
exit
