#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <string.h>
#include <unistd.h>
#include "gsl/gsl_randist.h"
using namespace std;

// { genename : {columnID: expression_value}
map<string,map<int,double>*> geneExpression;

vector<string> headers; // header row that matches the data columns
vector<string> samples; // Store the samples (rows) in the order in which they were read


int randseed;	// random seed (optional)

int readDataMatrix(const char*);
int generatePartition(int,const char*,int);
int generatePartition_Excl(int,const char*,int);

/*
* Generates n-fold cross validation data sets. 
* Prints both training and testing components as two files per fold.
* bool flag to randomize the columns before generating data sets.
*/
int generatePartition_Excl_TrainAndTest(int,const char*,bool doRand);
int printToFile(map<int,int>& usedInit, char* oFN, char* dFN);


int 
main(int argc, const char** argv)
{
	cout << argc << endl;
	if(argc!=5 && argc!=6)
	{
		cout <<"Usage: makeRowPartitions inputdata partitions outputdir partitiontype[foldsRand|foldsOrdered] [rand_seed]" << endl;
		return 0;
	}
	readDataMatrix(argv[1]);
	
	// user-specified seed
	if (argc==6)
	{
		randseed=atoi(argv[5]);
	} 
	else
	{
		randseed=getpid();
	}
	
	// original mode
	if(strcmp(argv[4],"foldsRand")==0)
	{
		generatePartition_Excl_TrainAndTest(atoi(argv[2]),argv[3],true);
	} 
	else if(strcmp(argv[4],"foldsOrdered")==0)
	{
		generatePartition_Excl_TrainAndTest(atoi(argv[2]),argv[3], false);
	} 
	else 
	{
		cout <<"Usage: makePartitions inputdata partitions outputdir partitiontype[foldsRand|foldsOrdered] [rand_seed]" << endl;
		cout <<"partitionsize only used by plain rand mode" <<endl;
		return 0;
	}
	return 0;
}

/**
* Reads data matrix into { geneName : {columnID, expressionValue}}
*/ 
int
readDataMatrix(const char* aFName)
{
	ifstream inFile(aFName);
	char* buffer=NULL;
	int bufflen=0;
	string buffstr;
	bool first=true; // set to false after reading header
	while(inFile.good())
	{
		getline(inFile,buffstr);
		if(buffstr.length()<=0)
		{
			continue;
		}
		if(bufflen<=buffstr.length())
		{
			bufflen=buffstr.length()+1;
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		map<int,double>* expMap=NULL;
		
		
		// read header row
		if (first)
		{
			while(tok!=NULL)
			{
				// skip the first token
				// read others as headers
				if(tokCnt>0)
				{
					string header(tok);
					headers.push_back(header);
				}
				tok=strtok(NULL,"\t");
				tokCnt++;
			}
			first=false;
		}
		else 
		{
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					string geneName(tok);
					expMap=new map<int,double>;
					geneExpression[geneName]=expMap;
					samples.push_back(geneName); // save in reading order
				}
				else
				{
					double aVal=atof(tok);
					(*expMap)[tokCnt-1]=aVal;
				}
				tok=strtok(NULL,"\t");
				tokCnt++;
			}
		}
	}
	inFile.close();
	return 0;
}

/*
* Leftover from makeColumnPartitions -- not updated
*/
int
generatePartition(int partitions, const char* outputDir,int subsetsize)
{
	cerr << "NOT IMPLEMENTED" << endl;
	return 1;
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	//int randseed=getpid(); // replaced by global var randseed
	gsl_rng_set(r,randseed);
	int size=geneExpression.begin()->second->size();
	double step=1.0/(double)size;
	char aFName[1024];
	int maxind=0;
	int minind=1000;
	for(int p=0;p<partitions;p++)
	{
		sprintf(aFName,"%s/dataindices%d.txt",outputDir,p);
		ofstream oFile(aFName);
		sprintf(aFName,"%s/dataset%d.txt",outputDir,p);
		ofstream dFile(aFName);
		map<int,int> usedInit;
		for(int i=0;i<subsetsize;i++)
		{
			double rVal=gsl_ran_flat(r,0,1);
			int rind=(int)(rVal/step);
			while(usedInit.find(rind)!=usedInit.end())
			{
				rVal=gsl_ran_flat(r,0,1);
				rind=(int)(rVal/step);
			}
			usedInit[rind]=0;
			if(rind<minind)
			{
				minind=rind;
			}
			if(rind>maxind)
			{
				maxind=rind;
			}
		}
		// print headers in the chosen order
		dFile << "Gene";
		for(map<int,int>::iterator aIter=usedInit.begin();aIter!=usedInit.end();aIter++)
		{
			dFile << "\t" << headers[aIter->first];
		}
		dFile << endl;
		
		// for each gene 
		for(map<string,map<int,double>*>::iterator gIter=geneExpression.begin();gIter!=geneExpression.end();gIter++)
		{
			dFile << gIter->first; // gene name
			map<int,double>* expMap=gIter->second;
			
			// print columns in the chosen order
			for(map<int,int>::iterator aIter=usedInit.begin();aIter!=usedInit.end();aIter++)
			{
				double aVal=(*expMap)[aIter->first];
				dFile <<"\t" << aVal;
			}
			dFile << endl;
		}
		
		// print out the random index files
		for(map<int,int>::iterator rIter=usedInit.begin();rIter!=usedInit.end();rIter++)
		{
			oFile << rIter->first << endl;
		}
		oFile.close();
		dFile.close();
		usedInit.clear();
	}
	cout <<"Maxind " << maxind << " Minind " << minind << endl;
	return 0;
}

/*
* Leftover from makeColumnPartitions -- not updated
*/
int
generatePartition_Excl(int partitions, const char* outputDir,int subsetsize)
{
	cerr << "NOT IMPLEMENTED" << endl;
	return 1;
	int size=geneExpression.begin()->second->size();
	int testsetsize=size/partitions;
	char aFName[1024];
	int maxind=0;
	int minind=1000;
	for(int p=0;p<partitions;p++)
	{
		sprintf(aFName,"%s/dataindices%d.txt",outputDir,p);
		ofstream oFile(aFName);
		sprintf(aFName,"%s/dataset%d.txt",outputDir,p);
		ofstream dFile(aFName);
		map<int,int> usedInit;
		int testSetStart=p*testsetsize;
		int testSetEnd=(p+1)*testsetsize;
		if(p==partitions-1)
		{
			testSetEnd=size;
		}
		for(int i=0;i<size;i++)
		{
			if(i>=testSetStart && i<testSetEnd)
			{
				continue;
			}
			usedInit[i];
		}
		// print headers in the chosen order
		dFile << "Gene";
		for(map<int,int>::iterator aIter=usedInit.begin();aIter!=usedInit.end();aIter++)
		{
			dFile << "\t" << headers[aIter->first];
		}
		dFile << endl;
		
		// for each gene, print values in the chosen order
		for(map<string,map<int,double>*>::iterator gIter=geneExpression.begin();gIter!=geneExpression.end();gIter++)
		{
			dFile << gIter->first;
			map<int,double>* expMap=gIter->second;
			for(map<int,int>::iterator aIter=usedInit.begin();aIter!=usedInit.end();aIter++)
			{
				double aVal=(*expMap)[aIter->first];
				dFile <<"\t" << aVal;
			}
			dFile << endl;
		}
		for(map<int,int>::iterator rIter=usedInit.begin();rIter!=usedInit.end();rIter++)
		{
			oFile << rIter->first << endl;
		}
		oFile.close();
		dFile.close();
		usedInit.clear();
	}
	//cout <<"Maxind " << maxind << " Minind " << minind << endl;
	return 0;
}

/**
* Generates partitions on the rows, prints out both the partition and its complement 
* as "testing" and "training" files.
* If doRand, will randomize rows prior to splitting into folds.
*/
int
generatePartition_Excl_TrainAndTest(int partitions, const char* outputDir, bool doRand)
{
	int size=geneExpression.size(); // number of genes
	int size2=samples.size(); // should match
	cout << "map size: " << size << ", read in " << size2 << endl;
	
	int testsetsize=size/partitions;
	char trainInFN[1024];
	char trainSetFN[1024];
	char testInFN[1024];
	char testSetFN[1024];
	
	int maxind=0;
	int minind=1000;
	
	// If we want to randomize the ROW order, we do it here.
	int* indices=new int[size];
	for (int i=0; i<size; i++)
	{
		indices[i]=i;
	}
	
	if (doRand)
	{
		gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
		// use global var randseed
		gsl_rng_set(r,randseed);
		gsl_ran_shuffle(r,indices,size,sizeof(int));
		cout << "Did randomization of samples." << endl;	
	} 
	
	for(int p=0;p<partitions;p++)
	{
		cout << "PARTITION " << p << endl;
		
		// maintain train and test indices for samples
		map<int,int> testInit;
		map<int,int> trainInit;
		
		int testSetStart=p*testsetsize;
		int testSetEnd=(p+1)*testsetsize;
		if(p==partitions-1)
		{
			testSetEnd=size;
		}
		
		// Pick train or test indices
		for(int i=0;i<size;i++)
		{
			int j=indices[i]; // we may have randomized the indices, so take from here.
			if(i>=testSetStart && i<testSetEnd)
			{
				testInit[j];
			}
			else
			{
				trainInit[j];
			}
		}
		cout << "\ttest size: " << testInit.size() << endl;
		cout << "\ttrain size: " << trainInit.size() << endl;
		
		// make output files
		sprintf(trainInFN,"%s/trainindices%d.txt",outputDir,p);
		sprintf(trainSetFN,"%s/trainset%d.txt",outputDir,p);
		sprintf(testInFN,"%s/testindices%d.txt",outputDir,p);
		sprintf(testSetFN,"%s/testset%d.txt",outputDir,p);
		
		printToFile(trainInit, trainInFN, trainSetFN);
		printToFile(testInit, testInFN, testSetFN);
	
		testInit.clear();	
		trainInit.clear();
	}
	//cout <<"Maxind " << maxind << " Minind " << minind << endl;
	
	delete [] indices;
	
	return 0;
}

/*
* Print data and indices to files.
* usedInit: indices for the set. (use key only)
* oFN: index file name
* dFN: data file name
*/
int
printToFile(map<int,int>& usedInit, char* oFN, char* dFN)
{
	// Print data file
	//cout << dFN << endl;
	ofstream dFile(dFN); 

	dFile << "Sample";
	for (int i=0; i< headers.size(); i++)
	{
		dFile << "\t" << headers[i];
	}
	dFile << endl;
		
	// now, for this partition..
	// loop over chosen gene order
	for(map<int,int>::iterator aIter=usedInit.begin();aIter!=usedInit.end();aIter++)
	{
		int geneId=aIter->first;
		string geneName=samples[geneId];
		map<int,double>* expmap=geneExpression[geneName];
		
		dFile << geneName;
		for (int j=0; j< expmap->size(); j++)
		{
			double aVal=(*expmap)[j];
			dFile <<"\t" << aVal;
		}
		dFile << endl;
	}
	dFile.close();
	
	
	//cout << oFN << endl;
	// print indices
	ofstream oFile(oFN); // index file
	for(map<int,int>::iterator rIter=usedInit.begin();rIter!=usedInit.end();rIter++)
	{
		oFile << rIter->first << endl;
	}
	oFile.close();
}


