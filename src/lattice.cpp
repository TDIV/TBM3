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
//  lattice.cpp
//  TBM^3
//
//  Created by Yuan Yen Tai on 7/02/15.
//


/* -------------------------------------------------------------------
 The BondNameVectorPair class will be used to be the element of bonding map.
 -------------------------------------------------------------------*/
class BondNameVectorPair{
public:
	string	bondName;
	r_mat	bondVector;
	BondNameVectorPair(string name, r_mat vec){
		bondName = name;
		bondVector = vec;
	}
};
/* -------------------------------------------------------------------
 The 'Lattice' class create the Lattice structure from the Atom class of atom.cpp.
 It can construct the basic bonding relations from the input file.
 This is convenient to construct the near by neighbors for the model (e.g. hopping, pairing).
 However, it cannot create the bonding relations 'beyond the input file.
 A complete bonding relation can get from the next class : LatticeBond.
 
 The 'Lattice' class will create two important products:
	1. atomList:		This is the place holder for all the atoms.
	2. atomPairList:	This store all the relations of neighboring atoms.
 Out side the 'Lattice' class, the 'TBModel' will use these two components to construct the Hamiltonian.
 -------------------------------------------------------------------*/
class Lattice{
private:
	H_SYMMETRY					symmetry;		//The required model symmetry : Normal, Nambu, EXNambu

	int							_index_size;
	
    vector<string>				strBasisVector;	//Basis vector: a1, a2, a3
    vector<string>      		subAtom;		//Sub atoms in a unitcell
    unsigned            		N1, N2, N3;		//Supercell Dim

    vector<string>				bondStringList;	/*
	The bond basis store the information of the basic bond label and vector information
		in the string formate (bondStringList.size() should be 3):
		example:
		bondStringList[0] = "x	0.5		0.0		0.0"
		bondStringList[1] = "y	0.0		0.5		0.0"
		bondStringList[2] = "z	0.0		0.0		0.5" */
	
	/* ----------------------------------------------------------------------------------
	 latticeBondIndexMapList[ AtomIndex_I ] [ bondName ] -> AtomIndex_J
	 ----------------------------------------------------------------------------------*/
	vector< map<string, int> >	latticeBondIndexMapList;
	
	/* ----------------------------------------------------------------------------------
	   Create the index of each atom in the lattice, also construct the atom-pair map list
	 ----------------------------------------------------------------------------------*/
	vector<Atom>		atomList;	 // Place holder for all the atoms
	vector<AtomPair>	atomPairList;// Pair up two neighboring atoms as : "AtomI - bondName - AtomJ"
	void				createIndex	()			{
		_index_size=0;
		int Index=-1;
		
		for (unsigned i=0 ; i<atomList.size() ; i++) {
			
			atomList[i].createIndexLabel(symmetry);
			for (unsigned j=0; j<atomList[i].indexLabel.size(); j++) {
			
				if (atomList[i].Name()!="VA" or atomList[i].Name()!="VC" or atomList[i].Name()!="BD") {
					
					Index++;
					string label=atomList[i].indexLabel[j];
					atomList[i].indexMap[label]=Index;
				
				} else {
					atomList[i].indexLabel.clear();
				}
			}
		}
		_index_size=Index+1;
		
		// Construct the pair iteration list
		for (unsigned i =0 ; i<atomList.size(); i++) {
			
			auto AtomI = operator()(i);
			
			auto bkey=getBondKeysFromAtomI(i);
			
			for (unsigned j=0; j<bkey.size(); j++) {
				bondKeyVectorMap[bkey[j]] = createBondVectorFromKey(bkey[j]);
				
				auto AtomJ = operator()(i,bkey[j]);
				
				AtomPair ap;
				ap.AtomI = AtomI;
				ap.AtomJ = AtomJ;
				ap.bond	= bkey[j];
				ap.bvec	= BVEC(bkey[j]).get();
				atomPairList.push_back( ap );
				
			}
			
		}
	
	}
	
	/* ----------------------------------------------------------------------------------
	   Create the Key-Vector mapping for each representing bond
	 ----------------------------------------------------------------------------------*/
	map<string, r_mat>	bondKeyVectorMap;	// The map : { "+1+0+0" : [[ 1, 0, 0 ]] }
	vector<string>		getBondKeysFromAtomI(unsigned int ii)	{
		vector<string> vstrs;
		for(auto imap: latticeBondIndexMapList[ii]) {
			vstrs.push_back(imap.first);
		}
		return vstrs;
	}
	
	/* ----------------------------------------------------------------------------------
	   Create the bondNameVector list, a name-vector mapping will be given inside the map
	 ----------------------------------------------------------------------------------*/
	vector<BondNameVectorPair>	bondNameVectorList;
	r_mat				createBondVectorFromKey(string bond)	{
		/* Convert the bonding information from the input 'bond' and
			return the correspond bond vector of it
			example:	+1+0+0 ---> [[ 1, 0, 0 ]]
						+1.. ---> [[ 1, 0, 0 ]]
		 */
		r_mat vec(3,1);
		vec(0,0)=0;
		vec(1,0)=0;
		vec(2,0)=0;
		
		// Analyze the bond component and split it into : 1. arrayNxyzStr, 2. arraySignStr.
		auto tmpBnd = bond;
		replaceAll(tmpBnd, "+", " ");
		replaceAll(tmpBnd, "-", " ");
		auto arrayNxyzStr = split(tmpBnd, " ");
		
		vector<char> arraySignStr;
		for (int i=0 ; i<bond.size(); i++) {
			if (bond[i] == '+' or bond[i] == '-') {
				arraySignStr.push_back(bond[i]);
			}
		}
		
		if (arrayNxyzStr.size() == arraySignStr.size()) { // Test for valid formate of 'bond'.
		
			// Construct the bond_vector
			for (unsigned i=0 ; i<arraySignStr.size(); i++) {
				if ( i< bondNameVectorList.size() ){
					switch (arraySignStr[i]) {
						case '+':
							vec = vec + StrToInt(arrayNxyzStr[i]) * bondNameVectorList[i].bondVector;
							break;
						case '-':
							vec = vec - StrToInt(arrayNxyzStr[i]) * bondNameVectorList[i].bondVector;
							break;
					}
				} else {
					string errorMsg = "Error, cannot convert the formate of bond vector: "+bond;
					cout<<errorMsg<<endl;
					throw errorMsg;
				}
			}
		} else {
			string errorMsg = "Error, cannot convert the formate of bond vector: "+bond;
			cout<<errorMsg<<endl;
			throw errorMsg;
		}
		
		return vec;
	}
	void				createBondVectorList()					{ // This method will be called once from the 'open' method.
		bondNameVectorList.clear();
		
		for (unsigned i=0 ; i<bondStringList.size(); i++) { // Construct the bond_vector
			
			istringstream iss(bondStringList[i]);
			string name;
			double x,y,z;
			iss>>name>>x>>y>>z;
			r_mat vec(3,1);
			vec[0]=x; vec[1]=y; vec[2]=z;
			
			bondNameVectorList.push_back(BondNameVectorPair(name, vec));
		}
	}

public:
	map<string, string> bondStringMap;
	string	translateBondString(string bondStr)		const		{
		string retBondStr = bondStr;
		
		// The first layer translation that could be defined by the user.
		auto it = bondStringMap.find(bondStr);
		if (it != bondStringMap.end()) { retBondStr = it->second; }
		
		// The second layer translation to replace the dot, '.', symbol to "+0".
		replaceAll(retBondStr, ".", "+0");
		
		return retBondStr;
	}
	
	string	filename;
	
	Lattice			()											{ _index_size=0; }
	Lattice			(string _filename, H_SYMMETRY sym)			{ open(_filename, sym); }
	~Lattice		()											{
		strBasisVector.clear();
		subAtom.clear();
		bondStringList.clear();
		atomList.clear();
		latticeBondIndexMapList.clear();
		bondNameVectorList.clear();
	}
	bool			open(string _filename, H_SYMMETRY _symmetry){
		_index_size = 0;
		filename = _filename;
		symmetry = _symmetry;
		
		strBasisVector.clear();
		subAtom.clear();
		bondStringList.clear();
		atomList.clear();
		latticeBondIndexMapList.clear();
		bondNameVectorList.clear();
		
		N1=0; N2=0; N3=0;
		
		vector<string>		lattice_neighbor;	//The neighboring relations
		
		string filename_lif = filename+".lif";
        ifstream infile(filename_lif.c_str());

		
		// ------ Read in the file and construct the parts
        if ( infile.is_open() ) {
			int flag=-1;
            string line;
            while ( getline(infile, line) ) {
                istringstream iss(line);
                string head;
                
                iss>>head;
                if (head=="")   {continue;}
				if (head=="#0") { flag=0; continue;}
                if (head=="#1") { flag=1; continue;}
                if (head=="#2") { flag=2; continue;}
                if (head=="#3") { flag=3; continue;}
                if (head=="#4") { flag=4; continue;}
                if (head=="#5") { flag=5; continue;}
                if (head=="#6") { flag=6; continue;}
                if (head=="#7") { flag=7; continue;}
                if (head=="#8") { flag=8; continue;}
                
                if (flag==1) { strBasisVector.push_back(line); }
                
				if (flag==2) { subAtom.push_back(line); }
				
                if (flag==3) {
                    istringstream iss2(line);
                    unsigned X=0,Y=0,Z=0;
                    iss2>>X>>Y>>Z;
                    N1=X; N2=Y; N3=Z;
                }
                
                if (flag==4) { bondStringList.push_back(line); }
                
                if (flag==7) {atomList.push_back( Atom(line) );}
				
				if (flag==8) {lattice_neighbor.push_back( line );}
            }
			
            infile.close();
			
			// ------ Assign atom_info --------------
			for (unsigned i=0 ; i<atomList.size() ; i++) {
				atomList[i].setAtomInfo(
									subAtom[atomList[i].subIndex()]
									);
			}
			
			// ------ Construct the latticeBondIndexMapList structure
			if ( lattice_neighbor.size() != 0 )
			if ( lattice_neighbor.size() == atomList.size() )
			for (unsigned i=0 ; i<lattice_neighbor.size() ; i++) {
				
				auto content = split(lattice_neighbor[i], " ");
				
				map<string, int> site_map;
				
				if (content.size() > 1)
				for (unsigned i=1 ; i<content.size(); i++) {
					auto sm=split(content[i],":");
					site_map[sm[0]]= StrToInt(sm[1]);
				}
				latticeBondIndexMapList.push_back(site_map);
			}

		createBondVectorList();
		createIndex();
        }
		else{
			string ErrorMsg = "Error, cannot find correspond .lif file for: "+filename+".";
			cout<<ErrorMsg<<endl;
			throw ErrorMsg;
			return false;
		}
		
		return true;
	}

	H_SYMMETRY		getSymmetry	()					const		{ return symmetry; }

	int				index_size	()					const		{return _index_size;}
	int				lattice_size()					const		{return atomList.size();}

	Atom			operator()	(unsigned int ii)				{ return atomList[ii]; }
	Atom			operator()	(unsigned int ii, string bond)	{
		if (bond=="+0+0+0"){
			return atomList[ii];
		}
		int jj=latticeBondIndexMapList[ii][bond];
		return atomList[jj];
	}
	Atom			operator()	(Atom si, string bond)			{
		return this->operator()(si.atomIndex(), bond);
	}
	r_mat			BVEC		(string bond)					{
		return bondKeyVectorMap[bond];
	}
	vector<Atom>	site_iteration	()				const		{ return atomList; }
	vector<AtomPair>pair_iteration	()				const		{ return atomPairList; }
	vector<string>	getBondVectorStr()				const		{
		/*This method will be depracated in the future.*/
		return bondStringList;
	}
	bool			hasBondKey		(string bond)	const		{
		auto it = bondKeyVectorMap.find(bond);
		if (it == bondKeyVectorMap.end()) { return false; }
		return true;
	}
	
	r_mat			basis_vec()						const		{
		r_mat vec(3,3);
		
		vector<int> N;
		N.push_back(N1);
		N.push_back(N2);
		N.push_back(N3);
		for (unsigned i=0; i<strBasisVector.size(); i++) {
			istringstream iss(strBasisVector[i]);
			string name;
			double x,y,z;
			iss>>name>>x>>y>>z;
			
			vec(i,0)=N[i]*x;
			vec(i,1)=N[i]*y;
			vec(i,2)=N[i]*z;
		}
	
		return vec;
	}
	r_mat			get_reciprocal()				const		{
		
		r_mat A(3,3); //Basis vector
		r_mat B(3,3); //Reciprocal lattice vector
		A=basis_vec().get();
		
		// Condition flags
		bool vflag0=true, vflag1=true, vflag2=true;
		if ( abs(A(0,0))==0 and abs(A(0,1))==0 and abs(A(0,2))==0) { vflag0=false; }
		if ( abs(A(1,0))==0 and abs(A(1,1))==0 and abs(A(1,2))==0) { vflag1=false; }
		if ( abs(A(2,0))==0 and abs(A(2,1))==0 and abs(A(2,2))==0) { vflag2=false; }
		
		if ( !vflag0 and !vflag1 and !vflag2) { cout<<"run time error: Basis vector not setted."<<endl; }
		
		// Conditions for 2-D lattice setup
		if ( !vflag0 and vflag1 and vflag2 ) {
			r_mat aa=curl(A.row(1), A.row(2));
			aa*(1/cdot(aa, aa));
			A(0,0) = aa[0]; A(0,1) = aa[1]; A(0,2) = aa[2];
		}
		
		if ( vflag0 and !vflag1 and vflag2 ) {
			r_mat aa=curl(A.row(2), A.row(0));
			aa*(1/cdot(aa, aa));
			A(1,0) = aa[0]; A(1,1) = aa[1]; A(1,2) = aa[2];
		}
		
		if ( vflag0 and vflag1 and !vflag2 ) {
			r_mat aa=curl(A.row(0), A.row(1));
			aa*(1/cdot(aa, aa));
			A(2,0) = aa[0]; A(2,1) = aa[1]; A(2,2) = aa[2];
		}
		
		
		double AVectorVolume = cdot(A.row(0), curl(A.row(1),A.row(2)));
		r_mat	b0=curl(A.row(1), A.row(2)); b0=b0*(2*pi/AVectorVolume);
		r_mat	b1=curl(A.row(2), A.row(0)); b1=b1*(2*pi/AVectorVolume);
		r_mat	b2=curl(A.row(0), A.row(1)); b2=b2*(2*pi/AVectorVolume);
		
		for (unsigned i=0; i<B.rows(); i++) {
			B(0,i)=b0[i]; B(1,i)=b1[i]; B(2,i)=b2[i];
		}
		
		
		// Conditions for 2-D lattice setup
		if ( !vflag0 and vflag1 and vflag2 ) { B(0,0) = 0; B(0,1) = 0; B(0,2) = 0; }
		if ( vflag0 and !vflag1 and vflag2 ) { B(1,0) = 0; B(1,1) = 0; B(1,2) = 0; }
		if ( vflag0 and vflag1 and !vflag2 ) { B(2,0) = 0; B(2,1) = 0; B(2,2) = 0; }
		
		return B;
	}
};

// -------------------------------------------------------------------
// LatticeBond is a complimentary part for the Lattice.
// It can create any bonding relations that if the shortest bonding relations are defined.
// To use this class, one need to define a "cut-off" for the maxima bond distance.
// Note that, up to now, this is only designed for the CUBIC-PEROVSKYTE system.
// More generaled version should be developed in the future.
// -------------------------------------------------------------------
class LatticeBond{
private:
	Lattice				Lat;
	
	vector<string>		bondStringList;
	map<string, r_mat>	bondKeyVectorMap;
	
	r_mat	get_vec(string bond_label){
		r_mat ret_vec(1,3);
		
		for (unsigned ii=0; ii<bondStringList.size(); ii++) {
			unsigned p_bond_count = WordCount(bond_label, "+"+bondStringList[ii]);
			unsigned m_bond_count = WordCount(bond_label, "-"+bondStringList[ii]);
			
			ret_vec = ret_vec + p_bond_count*1.0*bondKeyVectorMap[bondStringList[ii]];
			ret_vec = ret_vec - m_bond_count*1.0*bondKeyVectorMap[bondStringList[ii]];
		}
		return ret_vec;
	}
	
	unsigned	max_bond_level;
	double		cutoff_distanc;
	
	vector<pair<string, r_mat> >	atom_bond;
public:
	LatticeBond(Lattice _lat, unsigned _max_bond_level, double _cutoff_distance){
		Lat = _lat;
		max_bond_level = _max_bond_level;
		cutoff_distanc = _cutoff_distance;
		
		auto bbasis = Lat.getBondVectorStr();
		
		for (unsigned i =0 ; i<bbasis.size(); i++)			{
			auto word = split(bbasis[i], " "); /*
							Split each line lf the bondKeyVectorMap[i] = "x	0.5		0.0		0.0"
						*/
			
			r_mat	b_vec(1,3);
			auto b_name = word[0];
			b_vec[0] = StrToDouble(word[1]);
			b_vec[1] = StrToDouble(word[2]);
			b_vec[2] = StrToDouble(word[3]);

			bondKeyVectorMap[b_name] = b_vec;
			bondStringList.push_back(b_name);
		}
		
		// Construct the expended_bondKeyVectorMap for the full pack 3D bond-label and vec.
		vector<deque<string> > expended_bondKeyVectorMap(bondStringList.size());
		for (unsigned ii=0; ii<bondStringList.size(); ii++) {
			
			string m_bond="-"+bondStringList[ii]; // The minus sign for the bond basis.
			string m_bond_acc="";
			
			string p_bond="+"+bondStringList[ii]; // The plus  sign for the bond basis.
			string p_bond_acc="";
			
			for (unsigned i_lev=0; i_lev<max_bond_level; i_lev++){
				m_bond_acc+=m_bond+",";
				expended_bondKeyVectorMap[ii].push_front(m_bond_acc);
			}
			
			expended_bondKeyVectorMap[ii].push_back("o,");

			for (unsigned i_lev=0; i_lev<max_bond_level; i_lev++){
				p_bond_acc+=p_bond+",";
				expended_bondKeyVectorMap[ii].push_back(p_bond_acc);
			}
		}

		// Construct a full pack of 3D bond_label with vectors (that is : atom_bond)
		for (unsigned i0=0; i0<expended_bondKeyVectorMap[0].size(); i0++)
		for (unsigned i1=0; i1<expended_bondKeyVectorMap[1].size(); i1++)
		for (unsigned i2=0; i2<expended_bondKeyVectorMap[2].size(); i2++) {
			auto total_bond_label = expended_bondKeyVectorMap[0][i0]+expended_bondKeyVectorMap[1][i1]+expended_bondKeyVectorMap[2][i2];
			auto vec = Lat.BVEC(total_bond_label);

			r_mat total_bond_vec(1,3);
			total_bond_vec[0] = vec[0];
			total_bond_vec[1] = vec[1];
			total_bond_vec[2] = vec[2];
			
			double vec_magnitude = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
			
			if (vec_magnitude < _cutoff_distance) {
				atom_bond.push_back(make_pair(total_bond_label, total_bond_vec));
			}
		}
		
	}
	
	vector<pair<Atom, r_mat> > atom_bondings(Atom SiteI)	{
		
		vector<pair<Atom, r_mat> > ret_atom_bondings;
		
		for (unsigned ii=0; ii<atom_bond.size(); ii++) {
			auto bond_iter = split(atom_bond[ii].first, ",");
			
			Atom SiteJ = Lat(SiteI, "o");
			
			for (unsigned i=0; i<bond_iter.size(); i++){
				SiteJ = Lat(SiteJ, bond_iter[i]);
			}
			ret_atom_bondings.push_back(make_pair(SiteJ, atom_bond[ii].second));
		}

		return ret_atom_bondings;
	}
	
};






















