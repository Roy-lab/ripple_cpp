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
	int readExpressionSet(const char*);
	int mergeData(const char*);
	private:
	int readExpression(const char*);
	vector<map<string,string>*> geneExpSets;
	map<string,int> geneCnt;
	vector<int> measurements;
};
#endif
