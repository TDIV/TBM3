/*-----------------------------------------------------------|
 | Copyright (C) 2016 Yuan-Yen Tai, Hongchul Choi,            |
 |                    Jian-Xin Zhu                            |
 |                                                            |
 | This file is distributed under the terms of the GNU        |
 | General Public License. See the file `LICENSE' in          |
 | the root directory of the present distribution, or         |
 | http://www.gnu.org/copyleft/gpl.txt .                      |
 |                                                            |
 |-----------------------------------------------------------*/

//
//  wannier90_constructor.cpp
//  TBM^3
//
//  Created by Yuan Yen Tai on 5/20/16.
//
//

class HoppingBaseConstructor {
protected:
	string		name;
	
	unsigned	maxOrbital;
	
	// The formate of two bondings: {+1+0+0: (-1.0, 0.0)}
	vector<pair<string, x_var> > hoppingList;
	
public:
	
	HoppingBaseConstructor(){ }
	HoppingBaseConstructor(string _name, const Lattice & Lat, string filename) {
		
		name = _name;
		maxOrbital = 0;
		
		ifstream infile(filename.c_str());
		
		// --------------------
		if (infile.is_open()) {
			
			string line;
			
			while ( getline(infile, line)) {
				
				auto wordInLine = split(line," ");
				
				if (wordInLine.size() == 7) { // Analyze if the condition matched.
					
					string bondKey = "";
					for (int i=0 ; i<3 ; i++){
						if (wordInLine[i][0] != '-'){
							wordInLine[i] = "+" + wordInLine[i];
						}
						bondKey += wordInLine[i];
					}
					
					auto orbital_I = StrToInt(wordInLine[3]);
					auto orbital_J = StrToInt(wordInLine[4]);
					
					if (orbital_I > maxOrbital) { maxOrbital = orbital_I; }
					if (orbital_J > maxOrbital) { maxOrbital = orbital_J; }
					
					auto bondOperation = bondKey + ":" + wordInLine[3] +":" + wordInLine[4];
					
					
					// First, check if the bondKey exist in the Lat (Lattice) formate.
					// Then, assign the correspond hopping with the key 'bondOperation'.
					if (Lat.hasBondKey(bondKey)) {
						r_var realPart = StrToDouble(wordInLine[5]);
						r_var imagPart = StrToDouble(wordInLine[6]);
						auto hopVar = x_var(realPart,imagPart);
						hoppingList.push_back(pair<string, x_var>(bondOperation, hopVar) );
					}
				}
				else{
					string warningStr= "Warning, formate not matched for w90:"+line+" of file:"+filename+"!";
					cout<<warningStr<<endl;
				}
			}
		}
		else {
			// Faild to open the file
			string errorStr = "Error, faild to open the file:"+filename+" !";
			cout<<errorStr<<endl;
			throw errorStr;
		}
		
		infile.close();
	}
	
	string				getName()				const	{ return name; }
	
	// Access to the hopping list.
	unsigned			size()					const	{return hoppingList.size();}
	pair<string, x_var> operator[](unsigned i)	const	{return hoppingList[i];}
	
	HoppingBaseConstructor & operator= (const HoppingBaseConstructor & hop) {
		name = hop.name;
		maxOrbital = hop.maxOrbital;
		hoppingList = hop.hoppingList;
		return *this;
	}
};


class HoppingOrderParameter: public OrderParameter {
private:
	map<string, HoppingBaseConstructor>	hoppingBaseMap;
	
public:
	HoppingOrderParameter(Lattice & _Lat): OrderParameter(_Lat) {}
	
	void	appendHoppingBase(string name, string filename){
		hoppingBaseMap[name] = HoppingBaseConstructor(name, Lat, filename);
	}
	void	constructHoppingOrder(){
		// Clean all previous storages, if any.
		clear();
		
		cout<<"Warning, the HoppingOrderParameter only takes a single-phase input formate."<<endl;
		cout<<"The heterostructure hopping scheme should be implemented in the future."<<endl;
		cout<<endl;
		
		if (hoppingBaseMap.size() == 0) {
			string ErrorMsg = "Error, there is no 'HoppingBaseConstructor' inside the map, please 'appendHoppingBase' before calling this method (constructHoppingOrder).";
			cout<<ErrorMsg<<endl;
			throw ErrorMsg;
		}
		
		auto pair_iteration = Lat.pair_iteration();
		
		// Begin to construct the hopping terms
		for (int ii=0 ; ii<pair_iteration.size() ; ii++) {
			auto pit = pair_iteration[ii];
			
			auto atomI = pit.AtomI;
			auto atomJ = pit.AtomJ;
			
			// isHeterostructureHopping will be used later to construct the hopping list.
			bool	isSinglePhaseHopping = atomI.Name() == atomJ.Name();
			
			string	pairName = atomI.Name()+"-"+atomJ.Name();
			
			r_var	aSmallValue = 0.000000001;
			
			// Detect that if the pairName is inside the hoppingBaseMap.
			if (hoppingBaseMap.find(pairName) != hoppingBaseMap.end()) {
				
				// Select the founded hoppingBase
				auto & hoppingBase = hoppingBaseMap[pairName];
				
				for ( int hopIndex = 0 ; hopIndex < hoppingBase.size(); hopIndex++ ) {
					auto hopElement = hoppingBase[hopIndex];
					
					auto tmpStrParser = split( hopElement.first , ":" );
					
					if (isSinglePhaseHopping) {
						// ----------------------------------------------------
						// The single phase hopping operations goes here.
						// ----------------------------------------------------
						if (abs(hopElement.second.real()) > aSmallValue or abs(hopElement.second.imag()) > aSmallValue ) {
							order(pit, atomI.Name()+":"+atomJ.Name()+":"+tmpStrParser[0]+" hop."+tmpStrParser[1]+"."+tmpStrParser[2]) = hopElement.second;
						}
					}
					else {
						
						// ****************************************************
						// ****************************************************
						// The heterostructure hopping will be implement later.
						// ****************************************************
						// ****************************************************
					}
				}
			}
		}// End of construction of the hopping terms
		
		save("hopping");
	}
};


class	HoppingMatrixElement{
public:
	int		I;
	int		J;
	x_var	val;
	r_mat	bondVector;
	string	spaceDegreeLabel;
	
	HoppingMatrixElement(){ }
	HoppingMatrixElement(int i, int j, x_var _val, r_mat & _bondVecor, string spaceDegree){
		I=i;
		J=j;
		val=_val;
		bondVector = _bondVecor;
		spaceDegreeLabel = spaceDegree;
	}
	void set(int i, int j, x_var _val, r_mat & _bondVecor, string spaceDegree){
		I=i;
		J=j;
		val=_val;
		bondVector = _bondVecor;
		spaceDegreeLabel = spaceDegree;
	}
	
	
};

class HoppingModelConstructor{
private:
	Lattice & latticeRef;
	
protected:
	HoppingOrderParameter			hoppingOrder;
	vector<HoppingMatrixElement>	hoppingMatrixElementList;
	
public:
	HoppingModelConstructor(Lattice & _Lat): latticeRef(_Lat), hoppingOrder(_Lat) { }
	
	void appendHoppingBase(string name, string filename) { hoppingOrder.appendHoppingBase(name, filename); }
	
	void initHoppingTerms() {
		hoppingOrder.constructHoppingOrder();
		
		for (unsigned ii = 0 ; ii < hoppingOrder.pairListSize() ; ii++) {
			auto pairOperation = hoppingOrder.getPairOperationAt(ii);
			appendHoppingMatrixElement(pairOperation);
		}
	}
	
	void constructHoppingHamiltonian(x_mat & Ham, r_mat k_space){
		
		for( int i = 0 ; i < hoppingMatrixElementList.size() ; i++ ){
			auto indexI = hoppingMatrixElementList[i].I;
			auto indexJ = hoppingMatrixElementList[i].J;
			auto spaceDegree = hoppingMatrixElementList[i].spaceDegreeLabel;
			auto val = hoppingMatrixElementList[i].val;
			auto bvec = hoppingMatrixElementList[i].bondVector;
			
			x_var phase = k_space[0]*bvec[0]+ k_space[1]*bvec[1]+ k_space[2]*bvec[2];
			auto	exp_IJ=exp(-Im*phase);
			
			if (spaceDegree == "")	{ Ham(indexI, indexJ) += val*exp_IJ; }
			if (spaceDegree == "u")	{ Ham(indexI, indexJ) += val*exp_IJ; }
			if (spaceDegree == "d")	{ Ham(indexI, indexJ) += val*exp_IJ; }
		}
		
	}
	
private:
	void appendHoppingMatrixElement(PairOperationUnit & pairOpt){
		
		AtomPair & ap = pairOpt.pit;
		//cout<<ap.bond<<":"<<ap.bvec[0]<<" "<<ap.bvec[1]<<" "<<ap.bvec[2]<<endl;
		
		if (ap.AtomI.getSpinDegree() != ap.AtomJ.getSpinDegree() ) {
			string ErrorMsg = "Error, the spin-degree should be the same for construct a general Hopping model.\n";
			ErrorMsg += "Please check the lattice input file: "+ latticeRef.filename + ".";
			cout<<ErrorMsg<<endl;
			throw ErrorMsg;
		}
		string spinDegree = ap.AtomI.getSpinDegree();
		
		
		auto optSegmentA	= split( pairOpt.opt, " ");
		auto optSegmentA_B	= split( optSegmentA[1], ".");
		
		string orbital_I = optSegmentA_B[1];
		string orbital_J = optSegmentA_B[2];
		
		
		//HoppingMatrixElement hopElem;
		// ------------------------------------------------------
		// Normal space
		// ------------------------------------------------------
		int indexI = 0;
		int indexJ = 0;
		if( latticeRef.getSymmetry() == NORMAL ) {
			if(spinDegree == ""){	//
				indexI = ap.AtomI.getOrbitalIndex(orbital_I);
				indexJ = ap.AtomI.getOrbitalIndex(orbital_J);
				hoppingMatrixElementList.push_back(
												   HoppingMatrixElement(indexI, indexJ, pairOpt.val, ap.bvec, "")
												   );
			}
			else if(spinDegree == "s"){
				indexI = ap.AtomI.getOrbitalIndex(orbital_I+"u");
				indexJ = ap.AtomI.getOrbitalIndex(orbital_J+"u");
				hoppingMatrixElementList.push_back(
												   HoppingMatrixElement(indexI, indexJ, pairOpt.val, ap.bvec, "u")
												   );
				
				indexI = ap.AtomI.getOrbitalIndex(orbital_I+"d");
				indexJ = ap.AtomI.getOrbitalIndex(orbital_J+"d");
				hoppingMatrixElementList.push_back(
												   HoppingMatrixElement(indexI, indexJ, pairOpt.val, ap.bvec, "d")
												   );
			}
		}
		
		// ------------------------------------------------------
		// Nambu space
		// ------------------------------------------------------
		if( latticeRef.getSymmetry() == NAMBU ) {
			if(spinDegree == "s"){
				indexI = ap.AtomI.getOrbitalIndex(orbital_I+"Au");
				indexJ = ap.AtomI.getOrbitalIndex(orbital_J+"Au");
				hoppingMatrixElementList.push_back(
												   HoppingMatrixElement(indexI, indexJ, pairOpt.val, ap.bvec, "Au")
												   );
				
				indexI = ap.AtomI.getOrbitalIndex(orbital_I+"Bd");
				indexJ = ap.AtomI.getOrbitalIndex(orbital_J+"Bd");
				hoppingMatrixElementList.push_back(
												   HoppingMatrixElement(indexI, indexJ,-pairOpt.val, ap.bvec, "Bd") // negative value for particle-hole symmetry.
												   );
			}
		}
		
		// ------------------------------------------------------
		// Extended Nambu space
		// ------------------------------------------------------
		if( latticeRef.getSymmetry() == EXNAMBU ) {
			if(spinDegree == ""){
				indexI = ap.AtomI.getOrbitalIndex(orbital_I+"A");
				indexJ = ap.AtomI.getOrbitalIndex(orbital_J+"A");
				hoppingMatrixElementList.push_back(
												   HoppingMatrixElement(indexI, indexJ, pairOpt.val, ap.bvec, "A")
												   );
				
				indexI = ap.AtomI.getOrbitalIndex(orbital_I+"B");
				indexJ = ap.AtomI.getOrbitalIndex(orbital_J+"B");
				hoppingMatrixElementList.push_back(
												   HoppingMatrixElement(indexI, indexJ,-pairOpt.val, ap.bvec, "B") // negative value for particle-hole symmetry.
												   );
			}
			if(spinDegree == "s"){
				indexI = ap.AtomI.getOrbitalIndex(orbital_I+"Au");
				indexJ = ap.AtomI.getOrbitalIndex(orbital_J+"Au");
				hoppingMatrixElementList.push_back(
												   HoppingMatrixElement(indexI, indexJ, pairOpt.val, ap.bvec, "Au")
												   );
				
				indexI = ap.AtomI.getOrbitalIndex(orbital_I+"Bu");
				indexJ = ap.AtomI.getOrbitalIndex(orbital_J+"Bu");
				hoppingMatrixElementList.push_back(
												   HoppingMatrixElement(indexI, indexJ,-pairOpt.val, ap.bvec, "Bu") // negative value for particle-hole symmetry.
												   );
				
				indexI = ap.AtomI.getOrbitalIndex(orbital_I+"Ad");
				indexJ = ap.AtomI.getOrbitalIndex(orbital_J+"Ad");
				hoppingMatrixElementList.push_back(
												   HoppingMatrixElement(indexI, indexJ, pairOpt.val, ap.bvec, "Ad")
												   );
				
				indexI = ap.AtomI.getOrbitalIndex(orbital_I+"Bd");
				indexJ = ap.AtomI.getOrbitalIndex(orbital_J+"Bd");
				hoppingMatrixElementList.push_back(
												   HoppingMatrixElement(indexI, indexJ,-pairOpt.val, ap.bvec, "Bd") // negative value for particle-hole symmetry.
												   );
			}
		}
		
	}
	
};


























