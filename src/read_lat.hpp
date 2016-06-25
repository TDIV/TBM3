/*-----------------------------------------------------------|
 | Copyright (C) 2016 Yuan-Yen Tai, Hongchul Choi,           |
 |                    Jian-Xin Zhu                           |
 |                                                           |
 | This file is distributed under the terms of the GNU       |
 | General Public License. See the file `LICENSE' in         |
 | the root directory of the present distribution, or        |
 | http://www.gnu.org/copyleft/gpl.txt .                     |
 |                                                           |
 |-----------------------------------------------------------*/

//
//  read_lat.hpp
//  TBM^3
//
//  Created by Yuan Yen Tai on 7/02/15.
//

// The ParserBase will be used as a framework to construct for any read-in input header.
class ParserBase{
public:
	ParserBase(string _keyStr): keyString(_keyStr) {}
	string keyString;
	string operator()	(){
		return keyString;
	}
	
};

/*######################################################
 The following class will be read from xxx.lat file.
 It's full construction is located at Lat.open("xxx.lat");
 ######################################################*/
// Read in the parameter section
class Parameter:			public ParserBase{
	map<string, string>	strPool;
	map<string, x_var>	varPool;
	map<string, r_mat>	vectorPool;
	vector<string>		inputLines;
public:
	Parameter(): tbm::ParserBase("#Parameters"){ }
	
	void	append(string line){
		inputLines.push_back(line);
		
		auto lineParser = split(line, "=");
		if ( lineParser.size() == 2) {
			string parameterName =lineParser[0];
			removeSpace(parameterName);
			
			string parameterLine = lineParser[1];
			
			auto lineParser2nd = split(parameterLine, ",");
			
			if (lineParser2nd.size() == 1) {
				
				string singleParameter = lineParser2nd[0];
				removeSpace(singleParameter);
				
				if (singleParameter[0] == '\"' and singleParameter[singleParameter.size()-1] == '\"') {
					// The string type goes here.
					singleParameter.erase(0,1);
					singleParameter.pop_back();
					strPool[parameterName] = singleParameter;
				}
				else{
					// The variable type goes here.
					varPool[parameterName] = StrToComplex( singleParameter );
				}
			}
			else if (lineParser2nd.size() > 1) {
				// The complex vector type goes here
				
				r_mat	varVec(1, lineParser2nd.size()); // Make a vector of real-matrix type
				for( unsigned i=0 ; i<lineParser2nd.size() ; i++){
					varVec[i] = StrToDouble(lineParser2nd[i]);
				}
				vectorPool[parameterName] = varVec;
			}
		}
	}
	
	string &	STR(string key)						{
		if (strPool.find(key) == strPool.end()) {
			string errorStr = "Error, string type parameter: \""+key+"\" not defined.";
			ErrorMessage(errorStr);
		}
		return strPool[key];
	}
	string &	STR(string key, string defaultValue){
		if (strPool.find(key) != strPool.end()) {
			return strPool[key];
		}
		strPool[key] = defaultValue;
		return strPool[key];
	}
	x_var &		VAR(string key)						{
		if (varPool.find(key) == varPool.end()) {
			string errorStr = "Error, parameter: \""+key+"\" not defined.";
			ErrorMessage(errorStr);
		}
		return varPool[key];
	}
	x_var &		VAR(string key, x_var defaultValue)	{
		if (varPool.find(key) != varPool.end()) {
			return varPool[key];
		}
		varPool[key] = defaultValue;
		return varPool[key];
	}
	r_mat & 	VEC(string key)						{
		r_mat defaultVal(1,1);
		return VEC(key, defaultVal);
	}
	r_mat &		VEC(string key, r_mat defaultValue)	{
		if (vectorPool.find(key) != vectorPool.end()) {
			return vectorPool[key];
		}
		vectorPool[key] = defaultValue;
		return vectorPool[key];
	}
	
	string	getFileString()	{
		string	fileString = keyString +'\n';
		for( auto & line: inputLines){ fileString += line+'\n'; }
		fileString += '\n';
		return fileString;
	}
	void	clear()			{
		strPool.clear();
		varPool.clear();
		vectorPool.clear();
	}
};

// Read in the KSymmetryPoint section
class KSymmetryPoint:		public ParserBase{
public:
	vector<pair<string, r_mat> > kSymmPointList;
	
	KSymmetryPoint(): ParserBase("#KPointPath")	{ }
	void	append(string line)	{
		auto lineParser = split(line, " ");
		if( lineParser.size() == 4){
			r_mat	kpoint(1,3);
			string	label = lineParser[0];
			kpoint[0]	= StrToDouble(lineParser[1]);
			kpoint[1]	= StrToDouble(lineParser[2]);
			kpoint[2]	= StrToDouble(lineParser[3]);
			kSymmPointList.push_back(make_pair(label, kpoint));
		}
	}
	
	string	getFileString()		{
		
		string	fileString = keyString +'\n';
		for(auto & kpoint: kSymmPointList){
			fileString += fformat(kpoint.first,5)+" ";
			for(unsigned i=0 ; i<kpoint.second.size() ; i++){
				fileString += fformat(kpoint.second[i],5)+" ";
			}
			fileString += '\n';
		}
		fileString += '\n';
		return fileString;
	}
	void	clear()				{
		kSymmPointList.clear();
	}
};

// Read in the BasisVector section
class BasisVector:			public ParserBase{
	vector<r_mat>	A;
	r_var	abs(r_mat r){
		r_var sum = 0;
		for(unsigned i=0 ; i<r.size() ; i++){
			sum += r[i]*r[i];
		}
		return sqrt(sum);
		
	}
	r_var findMinmalRepeatanceOf(r_var N){
		vector<r_var>	radiusList;
		
		r_var minValue = 10000;
		if( A.size() == 3){
			for(int i=-1 ; i<=1 ; i++)
			for(int j=-1 ; j<=1 ; j++)
			for(int k=-1 ; k<=1 ; k++){
				if( i!=0 or j!=0 or k!=0){
					r_var value = abs( i*N*A[0] +j*N*A[1] +k*N*A[2]);
					if( minValue > value  ) minValue = value;
				}
			}
		}
		else if( A.size() == 2){
			for(int i=-1 ; i<=1 ; i++)
			for(int j=-1 ; j<=1 ; j++){
				if( i!=0 or j!=0){
					r_var value = abs( i*N*A[0] +j*N*A[1]);
					if( minValue > value  ) minValue = value;
				}
			}
		}
		else if( A.size() == 1){
			minValue = abs( N*A[0] );
		}

		return minValue;
	}
	vector<string>		inputLines;
public:
	BasisVector():ParserBase("#BasisVector"){  }
	
	void			append(string line)		{
		inputLines.push_back(line);
		
		auto lineParser = split(line, " ");
		
		if( lineParser.size() == 3){
			r_mat	a_vec(1, lineParser.size()); // Make a vector of real-matrix type
			for( unsigned i=0 ; i<lineParser.size() ; i++){
				a_vec[i] = StrToDouble(lineParser[i]);
			}
			A.push_back(a_vec);
		}
	}
	vector<r_mat>	getAVec()	const		{
		return A;
	}
	vector<r_mat>	getBVec()	const		{
		vector<r_mat> B; //Reciprocal lattice vector

		if (A.size() == 3) { // 3D
			auto tmpA = A;
			double AVectorVolume = cdot(tmpA[0], curl(tmpA[1],tmpA[2]));
			r_mat	b0=curl(tmpA[1], tmpA[2]); b0=b0*(2*pi/AVectorVolume);
			r_mat	b1=curl(tmpA[2], tmpA[0]); b1=b1*(2*pi/AVectorVolume);
			r_mat	b2=curl(tmpA[0], tmpA[1]); b2=b2*(2*pi/AVectorVolume);
			B.push_back(b0);
			B.push_back(b1);
			B.push_back(b2);
		}
		else if (A.size() == 2) { // 2D
			r_mat aa=curl(A[0], A[1]);
			aa = aa*(1/cdot(aa, aa));
			auto tmpA = A;
			tmpA.push_back(aa);
			
			double AVectorVolume = cdot(tmpA[0], curl(tmpA[1],tmpA[2]));
			r_mat	b0=curl(tmpA[1], tmpA[2]); b0=b0*(2*pi/AVectorVolume);
			r_mat	b1=curl(tmpA[2], tmpA[0]); b1=b1*(2*pi/AVectorVolume);
			r_mat	b2=curl(tmpA[0], tmpA[1]); b2=b2*(2*pi/AVectorVolume);
			B.push_back(b0);
			B.push_back(b1);
		}
		else if (A.size() == 1) { // 1D
			auto tmpA = A;
			double AVectorVolume = cdot(tmpA[0], tmpA[0]);
			r_mat	b0=tmpA[0]*(2*pi/AVectorVolume);
			B.push_back(b0);
		}
		
		return B;
	}
	
	unsigned	minRepeatForRadius(double radius)	{
		unsigned N=1;
		
		while( findMinmalRepeatanceOf(N) < radius){ N++; }
		return N;
	}
	string		getFileString(unsigned N1=1, unsigned N2=1, unsigned N3=1)	{
		string	fileString = keyString +'\n';
		//for( auto & line: inputLines){ fileString += line+'\n'; }
		//fileString += '\n';
		vector<unsigned> N;
		N.push_back(N1);
		N.push_back(N2);
		N.push_back(N3);
		for( unsigned i=0 ; i<A.size() && i<N.size() ; i++){
			auto exA = A[i]*N[i];
			for( unsigned ii=0 ; ii<exA.size() ; ii++){
				fileString += fformat(exA[ii],15)+" ";
			}
			fileString += '\n';
		}
		fileString += '\n';
		return fileString;
	}
	void		clear()								{
		A.clear();
	}
};

// Read in the BondVector section
class BondVector:			public ParserBase{
	map<unsigned, vector<r_mat> >	bondMap;
	
public:
	BondVector(): ParserBase("#BondVector"){ }
	void	append(unsigned mapKey, string line){
		
		auto it = bondMap.find(mapKey);
		if( it == bondMap.end() ){
			bondMap[mapKey] = vector<r_mat>();
		}
		
		auto lineParser = split(line, " ");
		
		if( lineParser.size() == 3 and bondMap[mapKey].size() < 3){
			r_mat	a_vec(1, lineParser.size()); // Make a vector of real-matrix type
			for( unsigned i=0 ; i<lineParser.size() ; i++){
				a_vec[i] = StrToDouble(lineParser[i]);
			}
			bondMap[mapKey].push_back(a_vec);
		}
	}
	vector<r_mat>	getBond(unsigned mapKey)	{
	auto it = bondMap.find(mapKey);
		if( it == bondMap.end()){
			ErrorMessage("Error, cannot find bondVector definition of: #"+IntToStr(mapKey));
		}
		return bondMap[mapKey];
	}
	
	string	getFileString(){
		string fileString = "";
		for( auto & it: bondMap){
			fileString += keyString +" "+IntToStr(it.first)+'\n';
			for( auto & vec: it.second){
				for( unsigned i=0 ; i<vec.size() ; i++){
					fileString += fformat(vec[i], 15)+" ";
				}
				fileString += '\n';
			}
			fileString += '\n';
		}
		return fileString;
	}
	void	clear(){
		bondMap.clear();
	}
};

// Read in the OrbitalProfile section
class OrbitalProfile:		public ParserBase{
	vector<deque<string> >	orbitalList;
	
	set<string>	validOrbital;
public:
	OrbitalProfile(): ParserBase("#OrbitalProfile"){
		validOrbital.clear();
		validOrbital.insert("s");
		validOrbital.insert("px");
		validOrbital.insert("py");
		validOrbital.insert("pz");
		validOrbital.insert("dxy");
		validOrbital.insert("dxz");
		validOrbital.insert("dyz");
		validOrbital.insert("dx2-y2");
		validOrbital.insert("dz2");
		
		//vector<string> vacancy;
		//vacancy.push_back("VC");
		//orbitalList.push_back(vacancy);
	}
	void	append(string line){
		auto lineParser = split(line, " ");
		
		if( lineParser.size() > 1){
			for( unsigned i=1; i<lineParser.size() ; i++){
				if( validOrbital.find(lineParser[i]) == validOrbital.end()){
					string ErrorMsg = "Error, "+lineParser[i]+" is not a valid orbital type.";
					ErrorMessage(ErrorMsg);
				}
			}
		}
		if( lineParser.size() > 0 ) orbitalList.push_back(lineParser);
	}
	
	vector<deque<string> > &	getOrbitalList()	{
		return	orbitalList;
	}
	
	string	getFileString(){
		
		string	fileString = keyString +'\n';
		for(auto & orbLine: orbitalList){
			for(auto & elem: orbLine){ fileString += fformat(elem, 10); }
			fileString += '\n';
		}
		fileString += '\n';
		return fileString;
	}
	
	void	clear(){
		orbitalList.clear();
	}
};

// Read in the Atoms section
class AtomStringParser:		public ParserBase{
public:
	vector<pair<unsigned, r_mat> >	atomInfoList;
	
	AtomStringParser(): ParserBase("#Atoms"){ }
	
	void	append(string line){
		auto lineParser = split(line, " ");
		if( lineParser.size() == 4){
			r_mat	pos(1,3);
			unsigned	index	= StrToInt(lineParser[0]);
			pos[0]	= StrToDouble(lineParser[1]);
			pos[1]	= StrToDouble(lineParser[2]);
			pos[2]	= StrToDouble(lineParser[3]);
			atomInfoList.push_back(make_pair(index, pos));
		}
	}
	void	clear(){
		atomInfoList.clear();
	}
};

/*######################################################
 The following class will be read from xxx.lat.tbm file.
 It's full construction is located at Lat.open("xxx.lat");
 ######################################################*/
// Read in the #CoreCharge section
class CoreCharge:	public ParserBase{
	map<string, r_var> coreChargeMap;	// Store the input from "xxx.lat.tbm".
public:
	CoreCharge(): ParserBase("#CoreCharge"){ }
	
	void	append(string line)	{
		auto parser = split(line, ">");
		if( parser.size() == 2){
			removeSpace(parser[0]);
			coreChargeMap[parser[0]] = StrToDouble(parser[1]);
		}
	}
	r_var	getCharge(string atomName){
		if( coreChargeMap.find(atomName) == coreChargeMap.end() ) return 0;
		return coreChargeMap[atomName];
	}
	void	clear()				{
		coreChargeMap.clear();
	}
};

// Read in the #Init section
class InitOrder:	public ParserBase{
public:
	vector<pair<string, x_mat> > orderOperationList;	// Store the input from "xxx.lat.tbm".
	
	InitOrder(): ParserBase("#Init"){ }
	
	void	append(string line)	{
		auto parser = split(line, ">");
		if( parser.size() == 2){
			orderOperationList.push_back(
				make_pair(parser[0], StrToXVec(parser[1]))
			);
		}
	}
	void	clear()				{
		orderOperationList.clear();
	}
};

// Read in the #Hamiltonian section
class HamiltonianParser:	public ParserBase{
public:
	vector< boost::tuple<string, string, string> >	hamOperationList;	// Store the input from "xxx.lat.tbm".
	map< string, vector<pair<string, string> > >	hamOperationMap;	// Store the input from "xxx.lat.tbm".
	
	HamiltonianParser(): ParserBase("#Hamiltonian"){ }
	
	void	append(string line)	{
		auto parser = split(line, ">");
		if( parser.size() == 3){
			
			removeSpace(parser[0]);
			replaceAll(parser[1], "\t", " ");
			removeSpace(parser[2]);
			hamOperationList.push_back( boost::make_tuple(parser[0], parser[1], parser[2]));
			
			if( hamOperationMap.find(parser[0]) == hamOperationMap.end() ){
				hamOperationMap[parser[0]] =  vector<pair<string, string> >();
				hamOperationMap[parser[0]].push_back(make_pair(parser[1], parser[2]));
			}
			else{
				hamOperationMap[parser[0]].push_back(make_pair(parser[1], parser[2]));
			}
		}
	}
	vector<pair<string, string> > & getOperationList(string strKey){
		if( hamOperationMap.find(strKey) == hamOperationMap.end() ){
			hamOperationMap[strKey] =  vector<pair<string, string> >();
		}
		return hamOperationMap[strKey];
	}
	
	void	clear()				{
		hamOperationList.clear();
	}
};














