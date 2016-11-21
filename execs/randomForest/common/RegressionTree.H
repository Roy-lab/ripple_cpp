/*The goal of this class is to store the Conditional distribution of a variable X given its neighbourhood N.*/
#ifndef _DECISIONTREE_
#define _DECISIONTREE_
#include <vector>
#include <queue>
#include "gsl/gsl_randist.h"
using namespace std;
#define PI 3.1472

#include "CommonTypes.H"

//typedef map<int,INTINTMAP*>  PARTITION;
//typedef map<int,INTINTMAP*>::iterator PARTITION_ITER;
typedef map<int,INTVECT*>  PARTITION;
typedef map<int,INTVECT*>::iterator PARTITION_ITER;
class EvidenceManager;
class Variable;
class Rule;
class RegressionTree
{
	public:
	RegressionTree();
	~RegressionTree();
	int setOutputVariable(int);
	int setEvidenceManager(EvidenceManager*);
	int setTestVariable(int);
	int setSubtreeVariable(int,const string&);
	int setDataID(int);
	int setPenalizedScore(double);
	double getPenalizedScore();

	double getTestValue();
	int getTestVariable();
	INTVECT& getDataSubset();
	INTINTMAP& getSubtreeVariables();

	int setMean(double);
	int setVariance(double);
	int setMarginalEntropy(double);

	int split(double testValue,PARTITION* p);
	int setChildParams(INTDBLMAP* childMean,INTDBLMAP* childVariance,INTDBLMAP* childEntropy);
	map<int,RegressionTree*>& getChildren();

	double getMarginalEntropy();
	RegressionTree* getChildAt(int);
	double getMarginalPDF(double);
	bool isPureNode();
	int clear();
	double getMean();
	double getVariance();
	
	// classification version -- store number of each outcome at leaf
	int addOutputValue(double output);
	
	// get probability of some output value (fraction)
	double getOutputProbability(double output);
	// get possible outcomes for this node (store in outcomes)
	int getOutputValues(vector<double>& outcomes);
	
	typedef enum
	{
		ROOT=0,
		NONLEAF=1,
		LEAF=2
	}NodeType;
	struct Pair
	{
		int ind;
		double val;
	};	

	int setType(NodeType);
	int setRNG(gsl_rng*);
	NodeType getType();
	int setParentInfo(int,int);
	int getParentInfo(int&,int&);
	int setParent(RegressionTree*);
	RegressionTree* getParent();

	int getLeafNodes(vector<RegressionTree*>&);
	int getTreeVars(INTINTMAP&);
	int showMe(string& starter,map<int,Variable*>&);
	int dumpTree(ostream&, map<int,Variable*>&,const char*);
	
	// Writes out tree to Cytoscape-friendly tab file
	// DEPRECATED -- use the one in common/Potential instead
	int dumpTreeToNetwork(ostream&, map<int,Variable*>&,const char*, map<int,int> &timesSeen);
	
	int generateRuleSet();
	vector<Rule*>& getRuleSet();
	int computeCodingLength(RegressionTree*);
	int setCodingLength(double);
	double getCodingLength();
	int genGenanatomy(ostream&,map<int,Variable*>&);
	
	// Wrapper around two functions - classification or regression
	int learn(double lambda,int leafSize,int varsToSample);
	
	// set mode of tree for learning 
	int setRegressionMode(bool doRegression);
	
	int deserialize(istream& iFile, map<int,Variable*>& varSet, RegressionTree* pp);
	int serialize(ostream& oFile, map<int,Variable*>& varSet);
	private:
	int estimateMarginal();

	// Called from learn() -- which one depends on the regressionMode flag
	int learnRegression(double lambda,int leafSize,int varsToSample);
	int learnClassification(double lambda,int leafSize,int varsToSample);

	int getDataAtLeaf(RegressionTree* currNode, vector<int>* evidVect);
	int getPartitions(int vId, double& splitValue, double marginalEntropy,INTINTMAP& dataSet, double& infoGain);
	int generateVarsTobeSampled(INTINTMAP& varSet,INTVECT& allVarSet,int varCnt);
	int getPartitions_Cached_Buggy(int vId, double& splitValue, double marginalEntropy,vector<int>& dataSet, double& infoGain);
	int getPartitions_Cached(int vId, double& splitValue, double marginalEntropy,vector<int>& dataSet, double& infoGain);
	
	// DC added - discrete target version
	int getPartitions_Discrete(int vId, double& splitValue, double marginalEntropy,vector<int>& dataSet, double& infoGain);
	
	
	int clearCache();
	int estimateMarginal(RegressionTree*);
	int sortAttrVals(DBLVECT&, INTVECT&);
	int sortAttrVals_Qsort(Pair**,int,DBLVECT&, INTVECT&);
	
	// Not used anymore?
	int getSubEntropy(double,double,double&);
	
	// Wrapper around discrete and continuous version
	int getSubEntropy(INTVECT&,int,int,double&,double&,double&);
	
	// For continuous outcomes
	int getSubEntropyContinuous(INTVECT&,int,int,double&,double&,double&);
	
	// For discrete outcomes
	int getSubEntropyDiscrete(INTVECT&,int,int,double&,double&,double&);
	int getShannonEntropy(map<double,int>& outcomes, int numItems, double& majority, double& entropy);
	int getGiniImpurity(map<double,int>& outcomes, int numItems, double& majority, double& impurity);
	
	double computeEntropyIfDeleted(Rule* r,int delMe,int&);

	int orderRulesByCoverage();

	RegressionTree* parent;
	map<int,RegressionTree*> children;
	//The set of datapoints in which the parent set configuration corresponding to the path
	//used to reach this node, holds true
	//This is for ease of computing the smaller CPTs
	//INTINTMAP dataSubset;
	vector<int> dataSubset;

	int testVarID;
	double testValue;
	
	//Need to come back and merge these two
	INTINTMAP subtreeVarIDs;
	INTSTRMAP subtreeVarIDNameMap;

	double mean;
	double variance;
	double marginalEntropy;
	
	// store the distribution of outcomes here
	// key is outcome; value is fraction of total
	map<double,double> outcomeDistro;
	
	NodeType nType;
	int parentNodeID;
	int parentBranch;
	double score;
	vector<Rule*> ruleSet;
	double codingLength;
	
	//Stuff like cache and current node set to be moved here from Potential
	queue<RegressionTree*> currentLeafNodes;
	map<int,PARTITION*> allPartitions;
	map<int,INTDBLMAP*> allMeans;
	map<int,INTDBLMAP*> allVariances;
	map<int,INTDBLMAP*> allMarginalEntropy;
	INTDBLMAP minValMap;
	INTDBLMAP maxValMap;
	int minLeafSize;
	int classVarID; // output variable ID
	EvidenceManager* evMgr;
	gsl_rng* r;
	
	// Run mode: regression (true) vs classification (false)
	bool regressionMode;
	
	// for classification mode -- hold onto the number of each kind of outcome
	map<double,int> outputValues;
};
#endif
