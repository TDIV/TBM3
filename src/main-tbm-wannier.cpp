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

#include "header.hpp"

#include <iostream>
#include <random>

string DoubleToPMStr(double val){
	string PM = "";
	if( val >= 0 ) PM = "+";
	return PM+tbm::DoubleToStr(val);
}

class TBWannierTranslator {
public:
	tbm::Parameter			parameter;
	tbm::BasisVector		basisVector;
	
	vector<string>			wannierBlockLines;
	vector<string>			atomSetupLines;
	
	// -----------------atomName, position, label(spinless) = 1,2,3 or label(spinfull) = 1u, 2u, 3u .... 1d, 2d, 3d ....
	vector<boost::tuple<string, tbm::r_mat, string> >	subAtomOrbitalList;
	
	string filename;
	TBWannierTranslator(string _filename){
		using namespace tbm;
		
		filename = _filename;
		
		basisVector.clear();
		
		parameter.clear();
		
		wannierBlockLines.clear();
		
		string line;
		string header = "";
		string flag = "";
		
		ifstream infile(_filename);
		if ( infile.good() )	{
			
			while ( getline(infile, line) ) {
				deleteComment(line); // Clean the commented words
				istringstream iss(line);
				iss >> header;
				
				if	( header == parameter())		{flag = header; continue;}
				if	( header == basisVector())		{flag = header; continue;}
				if	( header == "#AtomSetup")		{flag = header; continue;}
				if	( header == "#Wannier")			{flag = header; continue;}
				
				if	( flag == parameter())			{ parameter.append(line);			continue; }
				if	( flag == basisVector())		{ basisVector.append(line);			continue;	}
				if	( flag == "#AtomSetup")			{ atomSetupLines.push_back(line);	continue;	}
				if	( flag == "#Wannier")			{ wannierBlockLines.push_back(line);continue;	}
			}
		}
		else					{
			ErrorMessage("Error, input fiole not found: "+_filename);
		}
		infile.close();
		
		
		vector<string> orbitalProfile;
		vector<string> atomProfile;
		
		// Construct the subAtomOrbitalList.
		for( unsigned ii=0; ii < atomSetupLines.size() ; ii++){
			auto line = atomSetupLines[ii];
			auto parser = split(line, ">");
			
			if( parser.size() != 3 ) break;
			
			// Setup atom name
			removeSpace( parser[0] );
			string orbN = IntToStr(ii+1);
			string atomName = parser[0]+"-"+ orbN;
			
			atomProfile.push_back(fformat(IntToStr(ii+1))+" "+parser[1]);
			orbitalProfile.push_back(fformat(atomName)+" "+parser[2]);
			
			// Setup atom position
			tbm::r_mat pos(1,3);
			auto posStrList = split(parser[1], " ");
			for( unsigned p=0 ; p<posStrList.size(); p++){
				pos[p] = StrToDouble(posStrList[p]);
			}
			
			// Loop through orbital profile and arange in : subAtomOrbitalList.
			auto orbitalStrList = split(parser[2], " ");
			if( parameter.STR("spin", "on") == "off" ){
				for( unsigned i=1 ; i<=orbitalStrList.size() ; i++){
					string orbN = IntToStr(i);
					subAtomOrbitalList.push_back(boost::make_tuple( atomName, pos, orbN ));
				}
			}
			else
			if( parameter.STR("spin", "on") == "on" ){
				for( unsigned i=1 ; i<=orbitalStrList.size() ; i++){
					string orbN = IntToStr(i)+"u";
					subAtomOrbitalList.push_back(boost::make_tuple( atomName, pos, orbN ));
				}
				for( unsigned i=1 ; i<=orbitalStrList.size() ; i++){
					string orbN = IntToStr(i)+"d";
					subAtomOrbitalList.push_back(boost::make_tuple( atomName, pos, orbN ));
				}
			}
		}
		
		ofstream outfile(filename+".lat");
		outfile<<basisVector.getFileString();
		
		outfile<<"#OrbitalProfile"<<endl;
		for(auto & line: orbitalProfile){
			outfile<<line<<endl;
		}
		outfile<<endl;
		
		outfile<<"#Atoms"<<endl;
		for(auto & line: atomProfile){
			outfile<<line<<endl;
		}
		outfile.close();
	}
	
	void	parseWannierBlock(){
		using namespace tbm;
		
		auto avec = basisVector.getAVec();
		
		auto spin = parameter.STR("spin");
		auto bondN = parameter.VEC("bondN", tbm::xvec(20, 20, 20));
		double cvarCutoff = abs( parameter.VAR("precision", 0.00001).real() );
		cout<<"Converting with .."<<endl;
		cout<<"spin = "<<spin<<endl;
		cout<<"bondN = "<<bondN<<endl;
		cout<<"precision = "<<cvarCutoff<<endl<<endl;
		
		// ********************************
		// ***Converting parameter*********
		// ********************************
		vector<double> cellVectorDegeneracy;
		
		//-- cell N ----------------<PairIndex, cvar  , degeneracy>
		map<string, vector<boost::tuple<string, string, unsigned> > > cellVectorMap;
		set<string>	cellVecPool;
		
		unsigned flowControler = 0;
		// ********************************
	
		ofstream outfile(filename+".lat.tbm");
		outfile<<"#Hamiltonian"<<endl;
		for( auto & line: wannierBlockLines ){
			
			if		( flowControler == 0 )	{
				istringstream iss(line);
				double	w_var;
				
				while( iss>>w_var ){
					if( w_var > 0 and flowControler == 0){
						cellVectorDegeneracy.push_back(w_var);
					}
					else{
						flowControler = 1;
						break;
					}
				}
			}
			
			if (flowControler == 1){
				string cellVectorStr = line.substr(0, 16);
				
				if( cellVecPool.find(cellVectorStr) == cellVecPool.end()){
					cellVecPool.insert(cellVectorStr);
				}
				double degeneracy = cellVectorDegeneracy[cellVecPool.size()-1];
				
				istringstream iss(line);
				double N0, N1, N2, varReal, varImag;
				unsigned atomOrb_I, atomOrb_J;
				iss>> N0 >> N1 >> N2 >> atomOrb_I >> atomOrb_J >> varReal >> varImag ;
				
				if( abs(N0) > abs(bondN[0].real()) ) continue;
				if( abs(N1) > abs(bondN[1].real()) ) continue;
				if( abs(N2) > abs(bondN[2].real()) ) continue;
				
				varReal = varReal / degeneracy;
				varImag = varImag / degeneracy;
				
				// ****** Get sub atoms pair
				auto subAtomI = subAtomOrbitalList[atomOrb_I-1];
				auto subAtomJ = subAtomOrbitalList[atomOrb_J-1];
				
				// ****** Get sub atoms pair name
				auto subAtomI_name = subAtomI.get<0>();
				auto subAtomJ_name = subAtomJ.get<0>();
				
				// ****** Get bondIJ of sub atoms pair
				r_mat posI(1,3), posJ(1,3);
				posI = subAtomI.get<1>();
				posJ = N0 * avec[0] + N1 * avec[1] + N2 * avec[2] + subAtomJ.get<1>();
				auto bondIJ = posJ - posI;
				
				// ****** Get orbitals(spins) of sub atoms pair
				auto orbI = subAtomI.get<2>();
				auto orbJ = subAtomJ.get<2>();
				
				zvar val = varReal + Im*varImag;
				
				if( abs(val) > abs(cvarCutoff) ){
					if( parameter.STR("spin") == "off"){
						outfile<< "hopping  >  ";
					}
					if( parameter.STR("spin") == "on"){
						outfile<< "bond  >  ";
					}
					
					string bondOperation =
					subAtomI_name+":"+subAtomJ_name+":"
						+DoubleToPMStr( bondIJ[0] )
						+DoubleToPMStr( bondIJ[1] )
						+DoubleToPMStr( bondIJ[2] );
					
					string orbitalOperation = orbI+":"+orbJ;
					
					outfile << fformat(bondOperation, 44)<<" "
							<<fformat(orbitalOperation,7)<<"  >  "
							<<"("<<varReal<<","<<varImag<<")"
							<<endl;
				}
				
			}
		}
		
		outfile.close();
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



