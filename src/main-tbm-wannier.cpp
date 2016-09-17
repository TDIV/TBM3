/*-----------------------------------------------------------|
| Copyright (C) 2016 Yuan-Yen Tai, Hongchul Choi,            |
|                    Jian-Xin Zhu                            |
|                                                            |
| This file is distributed under the terms of the BSD        |
| Berkeley Software Distribution. See the file `LICENSE' in  |
| the root directory of the present distribution, or         |
| https://en.wikipedia.org/wiki/BSD_licenses.                |
|                                                            |
|-----------------------------------------------------------*/
//
//  main-tbm-wannier.cpp
//  TBM^3
//

#include "tbm3.hpp"

#include <iostream>
#include <random>

class TBWannierTranslator {
public:
	tbm::Parameter			parameter;
	tbm::BasisVector		basisVector;
	tbm::OrbitalProfile		orbitalProfile;
	tbm::AtomStringParser	atomParser;
	
	vector<string>			wannierBlockLines;
	
	TBWannierTranslator(string _filename){
		using namespace tbm;
		
		basisVector.clear();
		orbitalProfile.clear();
		atomParser.clear();
		
		parameter.clear();
		
		wannierBlockLines.clear();
		
		string line;
		string header = "";
		string flag = "";
		
		ifstream infile(_filename);
		if ( infile.good() ) {
			
			while ( getline(infile, line) ) {
				deleteComment(line); // Clean the commented words
				istringstream iss(line);
				iss >> header;
				
				if	( header == parameter())		{flag = header; continue;}
				if	( header == basisVector())		{flag = header; continue;}
				if	( header == orbitalProfile())	{flag = header; continue;}
				if	( header == atomParser())		{flag = header; continue;}
				if	( header == "#Wannier")			{flag = header; continue;}
				
				if	( flag == parameter())			{ parameter.append(line); continue;
				}
				if	( flag == basisVector())		{ basisVector.append(line);			continue;	}
				if	( flag == orbitalProfile())		{ orbitalProfile.append(line);		continue;	}
				if	( flag == atomParser())			{ atomParser.append(line);			continue;	}
				if	( flag == "#Wannier")			{ wannierBlockLines.push_back(line);continue;	}
			}
		}
		else{
			ErrorMessage("Error, cannot file file: "+_filename);
		}
		infile.close();
	}
	
	void	parseWannierBlock(){
		
		vector<double> cellVectorDegeneracy;
		
		map<string, vector<pair<string, string> > > cellVectorMap;
		
		unsigned flowControler = 0;
		
		for( auto & line: wannierBlockLines ){
			
			if		( flowControler == 0 )	{
				istringstream iss(line);
				double	w_var;
				
				while( iss>>w_var ){
					if( w_var < 0 )	flowControler = 1;
					else{
						cellVectorDegeneracy.push_back(w_var);
					}
				}
			}
			
			if (flowControler == 1){
				string cellVectorStr = line.substr(0, 16);
				string subAtomStr = line.substr(16, 10);
				string cvarStr = line.substr(26, 23);
				
				/* Formate output test.
				cout<< line <<endl;
				cout<< cellVectorStr <<endl;
				cout<< subAtomStr<<endl;
				cout<< cvarStr<<endl;
				cout<< cellVectorStr << subAtomStr << cvarStr <<" --- "<<endl;
				 */
				
				
			}
		}
		
		cout<<cellVectorDegeneracy.size()<<endl;
	}
};

int main(int argc, char *argv[]) {
	
	/* -------------------------------------------
	 Organizing/Collect arguments into a container.
	 ---------------------------------------------*/
	vector<string> args;
	for (int i=0; i<argc; i++) args.push_back(string(argv[i]));
	
	vector<string>	operationList;
	for(unsigned i=1 ; i<argc ; i++){ operationList.push_back(args[i]); }
	
	if( operationList.size() < 1 ){
		tbm::ErrorMessage("Error, no input file.");
	}
	
	/* -------------------------------------------
	 Feed the filename (from argument).
	 ---------------------------------------------*/
	string filename = operationList[0];
	TBWannierTranslator		wannierTranslator(filename);
	
	wannierTranslator.parseWannierBlock();
	
	return 0;
}







