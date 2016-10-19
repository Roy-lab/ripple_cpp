#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include "Framework.H"

Framework::Framework()
{
}

Framework::~Framework()
{
}

int 
Framework::readExpressionSet(const char* aFName) 
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
		if(strstr(buffer,"#")!=NULL)
		{
			continue;
		}
		readExpression(buffer);
	}
	inFile.close();
	return 0;
}

int 
Framework::mergeData(const char* aFName)
{
	ofstream oFile(aFName);

	map<string,int>::iterator gIter = geneCnt.find("Gene");
	if(gIter != geneCnt.end())
	{
		oFile << gIter->first;
		for(int s=0;s<geneExpSets.size();s++)
		{
			map<string,string>* exprs=geneExpSets[s];
			int measures=measurements[s];
			if(exprs->find(gIter->first)!=exprs->end())
			{
				oFile<<"\t" <<(*exprs)[gIter->first];
			}
			else
			{
				for(int i=0;i<measures;i++)
				{
					oFile <<"\t<nodata>";
				}
			}
		}
		oFile <<endl;
	}

	for(map<string,int>::iterator gIter=geneCnt.begin();gIter!=geneCnt.end();gIter++)
	{
		if(gIter->second < geneExpSets.size())
		{
			//continue;
		}
	
		if(strcmp(gIter->first.c_str(),"Gene")==0)
		{
			continue;
		}

		oFile << gIter->first;

		for(int s=0;s<geneExpSets.size();s++)
		{
			map<string,string>* exprs=geneExpSets[s];
			int measures=measurements[s];
			if(exprs->find(gIter->first)!=exprs->end())
			{
				oFile<<"\t" <<(*exprs)[gIter->first];
			}
			else
			{
				for(int i=0;i<measures;i++)
				{
					oFile <<"\t<nodata>";
				}
			}
		}
		oFile <<endl;
	}
	oFile.close();
	return 0;
}

int 
Framework::readExpression(const char* aFName)
{	
	ifstream inFile(aFName);
	string buffstr;
	char* buffer=NULL;
	int bufflen=0;
	map<string,string>* dataSet=new map<string,string>;
	int measures=0;
	while(inFile.good())
	{
		getline(inFile,buffstr);
		if(buffstr.length()<=0)
		{
			continue;
		}
		if(buffstr.length()>=bufflen)
		{
			bufflen=buffstr.length()+1;
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		char* tok=strchr(buffer,'\t');
		if(tok!=NULL)
		{
			*tok='\0';
		}
		string gene(buffer);
		string vals(tok+1);
		if(geneCnt.find(gene)==geneCnt.end())
		{
			geneCnt[gene]=1;
		}
		else
		{
			geneCnt[gene]=geneCnt[gene]+1;
		}
		(*dataSet)[gene]=vals;
		if(dataSet->size()==1)
		{
			char* tok2=strtok(tok+1,"\t");
			int tokcnt=0;
			while(tok2!=NULL)
			{
				tok2=strtok(NULL,"\t");
				tokcnt++;
			}
			measures=tokcnt;
		}
	}
	cout <<"Read " << dataSet->size() << " total genes" << endl; 
	geneExpSets.push_back(dataSet);
	measurements.push_back(measures);
	inFile.close();
	return 0;
}

int
main(int argc, const char** argv)
{
	if(argc!=3)
	{
		cout <<"mergeData infileset outputexpr" << endl;
		return 0;
	}
	Framework fw;
	fw.readExpressionSet(argv[1]);
	fw.mergeData(argv[2]);
	return 0;
}
