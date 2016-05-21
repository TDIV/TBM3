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

class HoppingConstructor {
protected:
	string		name;
	
	unsigned	maxOrbital;
	
	// The formate of two bondings: {+1+0+0: (-1.0, 0.0)}
	vector<pair<string, x_var> > hoppingList;
	
public:

	HoppingConstructor(){}
	HoppingConstructor(string _name, const Lattice & Lat, string filename) {
		
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
	
	HoppingConstructor & operator= (const HoppingConstructor & hop) {
		name = hop.name;
		maxOrbital = hop.maxOrbital;
		hoppingList = hop.hoppingList;
		return *this;
	}
};




class HoppingOrderParameter: OrderParameter {
private:
	map<string, HoppingConstructor>	hoppingSelectionMap;
	
public:
	HoppingOrderParameter(Lattice Lat): OrderParameter::OrderParameter(Lat) {}
	
	void	appendHoppingStructure(string name, HoppingConstructor hop){
		hoppingSelectionMap[name] = hop;
	}
	
	void	constructHoppingOrder(){
		// Clean all previous storages, if any.
		clear();
		
		if (hoppingSelectionMap.size() == 0) {
			string ErrorMsg = "Error, there is no 'HoppingConstructor' inside the map, please 'appendHoppingStructure' before calling this method (constructHoppingOrder).";
			cout<<ErrorMsg<<endl;
			throw ErrorMsg;
		}
		
		auto pair_iteration = Lat.pair_iteration();
		
		// Begin to construct the hopping terms from wan
		for (int ii=0 ; ii<pair_iteration.size() ; ii++) {
			auto pit = pair_iteration[ii];
			
			auto atomI = pit.AtomI;
			auto atomJ = pit.AtomJ;
			
			// isHeterostructureHopping will be use latter to construct the hopping list.
			bool	isHeterostructureHopping = atomI.Name() == atomJ.Name();
			
			string	pairName = atomI.Name()+"-"+atomJ.Name();
			
			// Detect that if the pairName is inside the hoppingSelectionMap.
			if (hoppingSelectionMap.find(pairName) != hoppingSelectionMap.end()) {
				
				auto & selectedHopConstructor = hoppingSelectionMap[pairName];
				for ( int hopIndex = 0 ; hopIndex < selectedHopConstructor.size(); hopIndex++ ) {
					auto hopElement = selectedHopConstructor[hopIndex];
					
					auto tmpStrParser = split( hopElement.first , ":" );
				}
			}
		}
	}
	
};


























