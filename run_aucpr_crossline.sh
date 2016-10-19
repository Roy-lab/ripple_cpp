#!/bin/bash
# Get results from cross-validation and run the AUCCalculator.
# starting with one training cell line.
# USAGE: gather_predictions datalist train_cellline results_top

if [ $# != 2 ]; then
	echo "USAGE: run_aucpr_crossline.sh cellname results_top"
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
mergeData=execs/mergeData/mergeData # /mnt/ws/sysbio/roygroup/shared/programs/mergedata/mergeData
collapseTP=execs/collapseTP/collapseTP #/mnt/ws/sysbio/roygroup/shared/programs/collapsetimepoints/collapseTP

# locate directories
traindir=${restop}/train_${trainline}/
if [ ! -d $traindir ];
then
	echo "Cannot find $traindir"
	exit 1
fi

# for each testing directory...
for testdir in ${traindir}/test_*
do
	testbase=$(basename $testdir)
	testline=${testbase/test_/}
	echo $testline	
	
	# for each tree setting
	for ntrees in 50 #1 10 50 100 500
	do
		echo $ntrees
		treesult=${testdir}/trees_${ntrees}
		echo $treesult
		if [ ! -d $treesult ]; then
			echo "Cannot find $treesult"
			continue
		fi
		if [ ! -d ${treesult}/fold_9 ]; then
			echo "Fold 9 not finished for $treesult"
			continue
		fi

		# locate the testset_err files 
		# extract just the relevant columns 
		# rename the headers
		
		# make file for merging columns
		colfile=${treesult}/column_merge.txt
		if [ -e $colfile ]; then
			rm $colfile
		fi
		mergefile=${treesult}/mergefile.txt
		if [ -e $mergefile ]; then
			rm $mergefile
		fi
		
		for i in {0..9};
		do
			folddir=${treesult}/fold_${i}
			if [ ! -d $folddir ]; then
				echo "CANNOT FIND $folddir"
				break
			fi
			errfile=${folddir}/testset_error.txt
			head -n 1 $errfile | cut -f 1,5,6 | sed "s/Column/Gene/; s/\t/\tFold${i}_/g;" > ${folddir}/err_for_merging.txt
				
			#cat ${folddir}/err_for_merging.txt
			tail -n +2 $errfile | cut -f 1,5,6 >> ${folddir}/err_for_merging.txt
			
			# get columns
			printf "Fold${i}_ConfIn_1\tConfIn_1\n" >> $colfile
			printf "Fold${i}_IsClass_1\tIsClass_1\n" >> $colfile
			
			echo ${folddir}/err_for_merging.txt >> $mergefile
		done
			
			
		# merge files
		mergedres=${treesult}/merged_results.txt
		$mergeData $mergefile $mergedres
		order=${treesult}/order.txt
		printf "ConfIn_1\nIsClass_1\n" > $order
		#echo $order
		
		# if any missing data, abort
		if grep -q "<nodata>" $mergedres
		then
			echo "Something wrong with test results for $testline, ${ntrees} -- found nodata in merged file."
			break
		fi
			
		# collapse to MEAN
		collapsed=${treesult}/collapsed_results.txt
			
		$collapseTP ${colfile} mean $order no ${mergedres} $collapsed
			
		# make list file out of it
		listfile=${treesult}/train_${trainline}_test_${testline}_trees${ntrees}.list
		aucfile=${listfile/.list/.auc}
		cut -f 2,3 $collapsed > $listfile

		echo $auc
		echo "Running AUPR tool on $listfile"
		java -jar $auc $listfile list &> ${aucfile}
		tail -n 2 ${aucfile}
	done # end trees
done #end testline
exit
