#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>

#include "Framework.H"

Framework::Framework()
{
	alphabetSort=false;
}

Framework::~Framework()
{
}

int
Framework::readGroupedTimepoints(const char* groupTP)
{
	ifstream inFile(groupTP);
	char buffer[1024];
	map<string,int> readGroups;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string exprName;
		string groupName;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				exprName.append(tok);
			}
			else if(tokCnt==1)
			{
				groupName.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		groupedTimepoints[exprName]=groupName;
		if(readGroups.find(groupName)==readGroups.end())
		{
			readGroups[groupName]=0;
		}
	}
	inFile.close();
	return 0;
}

int
Framework::setTransformationType(const char* aType)
{
	if(strcmp(aType,"mean")==0)
	{
		ptype=Framework::MEAN;
	}
	else if(strcmp(aType,"median")==0)
	{
		ptype=Framework::MEDIAN;
	}
	return 0;
}

int
Framework::readOrder(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string name(buffer);
		groupNames.push_back(name);
	}
	inFile.close();
}

int 
Framework::readExpressionData(const char* groupTP)
{
	ifstream inFile(groupTP);
	string buffstr;
	char* buffer=NULL;
	int bufflen=0;
	int lineCnt=0;
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
				delete [] buffer;
			}
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		if(lineCnt==0)
		{
			readColumnHeaders(buffer);
		}
		else
		{
			char* tok=strtok(buffer,"\t");
			int tokCnt=0;
			string geneName;
			map<string,double> expVals;
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					geneName.append(tok);
				}
				else
				{
					double eVal=-999;
					if(strlen(tok)>0 && strstr(tok,"nodata")==NULL && strstr(tok,"no_data")==NULL)
					{
						eVal=atof(tok);
					}
					expVals[headerNames[tokCnt-1]]=eVal;
				}
				tok=strtok(NULL,"\t");
				tokCnt++;
			}
			map<string,double>* newExpVals=new map<string,double>;
			if(strcmp(geneName.c_str(),"YAL001C")==0)
			{
				cout << "Found gene " << endl;
			}
			if(ptype==Framework::MEAN)
			{
				collapseData_Mean(expVals,newExpVals);
			}
			else if(ptype==Framework::MEDIAN)
			{
				collapseData_Median(expVals,newExpVals);
			}
	//		cout <<"Processing data for " << geneName << endl;
			geneExp[geneName]=newExpVals;
			expVals.clear();
		}
		lineCnt++;
	}
	inFile.close();
	return 0;
}

int
Framework::readColumnHeaders(char* buffer)
{
	char* tok=strtok(buffer,"\t");
	int tokCnt=0;
	while(tok!=NULL)
	{
		if(tokCnt>0)
		{
			string colName(tok);
			headerNames[tokCnt-1]=colName;
		}
		tok=strtok(NULL,"\t");
		tokCnt++;
	}
	if(alphabetSort)
	{
		map<string,int> tempOrder;
		cout <<"Previous order" << endl;
		for(int i=0;i<headerNames.size();i++)
		{
			cout << headerNames[i] << endl;
			tempOrder[headerNames[i]]=0;
		}
		headerNames.clear();
		int colID=0;
		cout <<"Reordered" << endl;
		for(map<string,int>::iterator aIter=tempOrder.begin();aIter!=tempOrder.end();aIter++)
		{
			headerNames[colID]=aIter->first;
			colID++;
		}
		for(int i=0;i<headerNames.size();i++)
		{
			cout << headerNames[i] << endl;
		}
	}
	return 0;
}


int 
Framework::collapseData_Mean(map<string,double>& exprVal, map<string,double>* collapsedExprVal)
{
	double currSum=0;
	double tpCnt=0;
	string currGroupName;
	//for(map<string,double>::iterator eIter=exprVal.begin();eIter!=exprVal.end();eIter++)
	for(map<int,string>::iterator hIter=headerNames.begin();hIter!=headerNames.end();hIter++)
	{
		if(groupedTimepoints.find(hIter->second)==groupedTimepoints.end())
		{
			cout <<"No timepoint map for " << hIter->second << endl;
			continue;
		}
		double eVal=exprVal[hIter->second];
		string& tpGroup=groupedTimepoints[hIter->second];
		if(currGroupName.length()==0)
		{
			currGroupName.append(tpGroup.c_str());
		}
		else
		{
			if(strcmp(currGroupName.c_str(),tpGroup.c_str())!=0)
			{
				//Time to update the entry
				double collapsedVal=currSum/tpCnt;
				(*collapsedExprVal)[currGroupName]=collapsedVal;
				currSum=0;
				tpCnt=0;
				currGroupName.clear();
				currGroupName.append(tpGroup.c_str());
			}
		}
		//We assume eval=-999 is for missing value
		if(eVal>-999)
		{
			currSum=currSum+eVal;
			tpCnt++;
		}
	}
	if(currGroupName.length()>0 && tpCnt>0)
	{
		double collapsedVal=currSum/tpCnt;
		(*collapsedExprVal)[currGroupName]=collapsedVal;
	}
	return 0;
}

int 
Framework::collapseData_Median(map<string,double>& exprVal, map<string,double>* collapsedExprVal)
{
	vector<double> currValues;
	double tpCnt=0;
	string currGroupName;
	//for(map<string,double>::iterator eIter=exprVal.begin();eIter!=exprVal.end();eIter++)
	for(map<int,string>::iterator hIter=headerNames.begin();hIter!=headerNames.end();hIter++)
	{
		if(groupedTimepoints.find(hIter->second)==groupedTimepoints.end())
		{
			cout <<"No timepoint map for " << hIter->second << endl;
			continue;
		}
		double eVal=exprVal[hIter->second];
		string& tpGroup=groupedTimepoints[hIter->second];
		if(currGroupName.length()==0)
		{
			currGroupName.append(tpGroup.c_str());
		}
		else
		{
			double collapsedVal=0;
			if(strcmp(currGroupName.c_str(),tpGroup.c_str())!=0)
			{
				//Time to update the entry
				for(int i=0;i<currValues.size();i++)
				{
					for(int j=i+1;j<currValues.size();j++)
					{
						double aval=currValues[i];
						double bval=currValues[j];
						if(aval>bval)
						{
							currValues[i]=bval;
							currValues[j]=aval;
						}
					}
				}
				if((currValues.size()%2)==1)
				{	
					int midpoint=currValues.size()/2;
					collapsedVal=currValues[midpoint];
				}
				else
				{
					int midpoint=currValues.size()/2;
					if(midpoint==0)
					{
						cout << "Only one sample!"<< currValues.size() << endl;
					}
					collapsedVal=(currValues[midpoint]+currValues[midpoint-1])/2;
				}
				(*collapsedExprVal)[currGroupName]=collapsedVal;
				tpCnt=0;
				currValues.clear();
				currGroupName.clear();
				currGroupName.append(tpGroup.c_str());
			}
		}
		currValues.push_back(eVal);
		tpCnt++;
	}
	if(currValues.size()>0)
	{
		double collapsedVal=0;
		//Time to update the entry
		for(int i=0;i<currValues.size();i++)
		{
			for(int j=i+1;j<currValues.size();j++)
			{
				double aval=currValues[i];
				double bval=currValues[j];
				if(aval>bval)
				{
					currValues[i]=bval;
					currValues[j]=aval;
				}
			}
		}
		if((currValues.size()%2)==1)
		{	
			int midpoint=currValues.size()/2;
			collapsedVal=currValues[midpoint];
		}
		else
		{
			int midpoint=currValues.size()/2;
			collapsedVal=(currValues[midpoint]+currValues[midpoint-1])/2;
		}
		(*collapsedExprVal)[currGroupName]=collapsedVal;
	}
	return 0;
}




int
Framework::showExpression(const char* aFName)
{
	ofstream oFile(aFName);
	oFile <<"Gene";
	for(int t=0;t<groupNames.size();t++)
	{
		oFile<<"\t" << groupNames[t];
	}
	oFile << endl;
	for(map<string,map<string,double>*>::iterator eIter=geneExp.begin();eIter!=geneExp.end();eIter++)
	{
		oFile <<eIter->first;
		map<string,double>* exprVals=eIter->second;
		for(int t=0;t<groupNames.size();t++)
		{
			if(exprVals->find(groupNames[t])==exprVals->end())
			{
				cout <<"Yikes! No expression for " << groupNames[t] << " for gene " <<eIter->first << endl;
			}
			double expVal=(*exprVals)[groupNames[t]];
			oFile<< "\t" << expVal;
		}	
		oFile << endl;
	}
	oFile.close();
	return 0;
}

int
Framework::sortAlphabetically()
{
	alphabetSort=true;
	return 0;
}

int
main(int argc, const char** argv)
{
	if(argc!=7)
	{
		cout <<"Usage: collapseTimepoints timepointgroups collapsetype[mean|median] expressionorder order_input_alphebatically[yes|no] expressionlevels outputexpressionlevels" << endl;
		return 0;
	}
	Framework fw;
	fw.readGroupedTimepoints(argv[1]);
	fw.setTransformationType(argv[2]);
	fw.readOrder(argv[3]);
	if(strcmp(argv[4],"yes")==0)
	{
		fw.sortAlphabetically();
	}
	fw.readExpressionData(argv[5]);
	fw.showExpression(argv[6]);
	return 0;
}
