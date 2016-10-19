#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <map>
#include <string>
#include <vector>
using namespace std;

class Framework
{
	public:
	Framework();
	~Framework();
	int readGroupedTimepoints(const char*);
	int setTransformationType(const char*);
	int readExpressionData(const char*);
	int readOrder(const char*);
	int collapseData_Mean(map<string,double>&, map<string,double>*);
	int collapseData_Median(map<string,double>&, map<string,double>*);
	int showExpression(const char*);
	int sortAlphabetically();
	typedef enum
	{
		MEAN,
		MEDIAN
	} PredictionType;
	private:
	int readColumnHeaders(char*);
	map<string,map<string,double>*> geneExp;
	map<string,string> groupedTimepoints;
	vector<string> groupNames;
	map<int,string> headerNames;
	PredictionType ptype;
	bool alphabetSort;
};
#endif
