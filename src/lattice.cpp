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
 The Lattice class create the Lattice structure from the Atom class of atom.cpp.
 It can construct the basic bonding relations from the input file.
 This is convenient to construct the near by neighbors for the model (e.g. hopping, pairing).
 However, it cannot create the bonding relations 'beyond the input file.
 A complete bonding relation can get from the next class : LatticeBond.
 -------------------------------------------------------------------*/
class BondVector{
public:
	string	bndName;
	r_mat	bndVec;
	BondVector(string name, r_mat vec){
		bndName = name;
		bndVec = vec;
	}
};

class Lattice{
private:
	H_SYMMETRY					symmetry;    //The required model symmetry : Normal, Nambu, EXNambu

    vector<string>				basis_vector;//Basis vector for the unitcell (Will be constructed from entry level 1).
    vector<string>      		sub_atom;    //Sub atoms in a unitcell
    unsigned            		N1, N2, N3;  //Supercell Dim

    vector<string>				bond_basis;			 /*
	The bond basis store the information of the basic bond label and vector information
		in the string formate (bond_basis.size() should be 3):
		example:
		bond_basis[0] = "x	0.5		0.0		0.0"
		bond_basis[1] = "y	0.0		0.5		0.0"
		bond_basis[2] = "z	0.0		0.0		0.5" */
	vector<BondVector>			bondVector;
    vector<string>				pair_operation_basis;/*
	The pair_operation_basis could be the same as the bond_basis, if it is not defined in the lattice input file. 
													*/
	
    vector<string>      		bonding;            //Bonding relations

	vector<Atom>				lattice_atom;		//Lattice atoms
	vector< map<string, int> >	lattice_map;		//Atom to Atom Neighbor relation map
	vector<AtomPair>			vec_pair;

	int							_index_size;

	void			createIndex	(H_SYMMETRY sym)				{
		_index_size=0;
		int Index=-1;
		
		for (unsigned i=0 ; i<lattice_atom.size() ; i++) {
			
			lattice_atom[i].createIndexLabel();
			for (unsigned j=0; j<lattice_atom[i].index_label.size(); j++) {
			
				if (lattice_atom[i].Name()!="VA" or lattice_atom[i].Name()!="VC" or lattice_atom[i].Name()!="BD") {
					
					Index++;
					string label=lattice_atom[i].index_label[j];
					lattice_atom[i].index[label]=Index;
				
				} else {
					lattice_atom[i].index_label.clear();
				}
			}
		}
		_index_size=Index+1;
		
		// Construct the pair iteration list
		for (unsigned i =0 ; i<lattice_atom.size(); i++) {
			
			auto AtomI = operator()(i);
			
			auto bkey=bond_key(i);
			
			for (unsigned j=0; j<bkey.size(); j++) {
				auto AtomJ = operator()(i,bkey[j]);
				
				AtomPair ap;
				ap.AtomI = AtomI;
				ap.AtomJ = AtomJ;
				ap.bond	= bkey[j];
				ap.bvec	= BVEC(bkey[j]).get();
				vec_pair.push_back( ap );
			}
			
		}
	
	}
	vector<string>	bond_key	(unsigned int ii)				{
		vector<string> vstrs;
		for(auto imap: lattice_map[ii]) {
			vstrs.push_back(imap.first);
		}
		
		return vstrs;
	}
	void			createBondVector(){ // This method will be called once from the 'open' method.
		bondVector.clear();
		
		for (unsigned i=0 ; i<bond_basis.size(); i++) { // Construct the bond_vector
			
			istringstream iss(bond_basis[i]);
			string name;
			double x,y,z;
			iss>>name>>x>>y>>z;
			r_mat vec(3,1);
			vec[0]=x; vec[1]=y; vec[2]=z;
			
			bondVector.push_back(BondVector(name, vec));
		}
	}
	
public:
	string filename;

public:
	Lattice			()											{ _index_size=0; }
	Lattice			(string _filename, H_SYMMETRY sym)			{
		_index_size = 0;
		filename = _filename;
		symmetry = sym;
		open(_filename+".lif", sym);
	}
	~Lattice		()											{
		basis_vector.clear();
		sub_atom.clear();
		bond_basis.clear();
		bonding.clear();
		lattice_atom.clear();
		lattice_map.clear();
	}
	bool			open(string filename, H_SYMMETRY sym)		{
		
		N1=0; N2=0; N3=0;
		basis_vector.clear();
		sub_atom.clear();
		bond_basis.clear();
		bonding.clear();
		lattice_atom.clear();
		lattice_map.clear();
		
		vector<string>		lattice_neighbor;	//The neighboring relations
		
        ifstream infile(filename.c_str());

		
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
                
                if (flag==1) { basis_vector.push_back(line); }
                
				if (flag==2) { sub_atom.push_back(line); }
				
                if (flag==3) {
                    istringstream iss2(line);
                    unsigned X=0,Y=0,Z=0;
                    iss2>>X>>Y>>Z;
                    N1=X; N2=Y; N3=Z;
                }
                
                if (flag==4) { bond_basis.push_back(line); }
                
                //if (flag==5) {
                //    istringstream iss2(line);
                //    string bds;
                //    while (iss2>>bds) { bonding.push_back(bds); }
                //}
                
                if (flag==7) {lattice_atom.push_back( Atom(line) );}
				
				if (flag==8) {lattice_neighbor.push_back( line );}
            }
			
			if (pair_operation_basis.size() == 0) {
				pair_operation_basis = bond_basis;
			}
            
            infile.close();
			
			// ------ Assign atom_info --------------
			for (unsigned i=0 ; i<lattice_atom.size() ; i++) {
				lattice_atom[i].setAtomInfo(
									sub_atom[lattice_atom[i].sub_index()]
									);
			}
			
			// ------ Construct the lattice_map structure
			if ( lattice_neighbor.size() != 0 )
			if ( lattice_neighbor.size() == lattice_atom.size() )
			for (unsigned i=0 ; i<lattice_neighbor.size() ; i++) {
				
				auto content = split(lattice_neighbor[i], " ");
				
				map<string, int> site_map;
				
				if (content.size() > 1)
				for (unsigned i=1 ; i<content.size(); i++) {
					auto sm=split(content[i],":");
					site_map[sm[0]]= StrToInt(sm[1]);
				}
				lattice_map.push_back(site_map);
			}

			createBondVector();
		createIndex(sym);
        }
		else{
			return false;
		}
		
		
		return true;
	}
	
	string pairOperationTranslator(string opt){
		replaceAll(opt, ".", "+0");
		return opt;
	}

	H_SYMMETRY		getSymmetry	()								{ return symmetry; }
	void			PrintIndex	()								{
		//**** Print ****
		cout<<endl<<endl;
		cout<<"Index size: "<<index_size()<<endl;
		for (unsigned i=0 ; i<lattice_atom.size() ; i++) {
			cout<<lattice_atom[i].Name()<<" ";
			for (unsigned j=0; j<lattice_atom[i].index_label.size(); j++) {
				string lab=lattice_atom[i].index_label[j];
				cout<<lab<<":"<<lattice_atom[i].index[lab]<<" ";
			}
			cout<<endl;
		}
	}
	int				index_size	()								{return _index_size; }
	int				lattice_size()								{ return lattice_atom.size(); }
	int				count_atoms	(string atom_list)				{
		
		int total = 0;
		auto atom_names = split(atom_list, " ");
		for (unsigned i=0 ; i<lattice_atom.size(); i++) {
			auto Name = lattice_atom[i].Name();
			for (unsigned j=0; j<atom_names.size(); j++) {
				if (Name == atom_names[j]) {
					total++;
					break;
				}
			}
		}
		return total;
	}

	Atom			operator()	(unsigned int ii)				{ return lattice_atom[ii]; }
	Atom			operator()	(unsigned int ii, string bnd)	{
		if (bnd=="o"){
			return lattice_atom[ii];
		}
		int jj=lattice_map[ii][bnd];
		return lattice_atom[jj];
	}
	Atom			operator()	(Atom si, string bnd)			{
		if (bnd=="o"){
			return lattice_atom[si.AtomIndex()];
		}
		int jj=lattice_map[si.AtomIndex()][bnd];
		return lattice_atom[jj];
	}
	r_mat			BVEC		(string bnd)					{
		/* Convert the bonding information from the input 'bnd' and
			return the correspond bond vector of it
			example:	+1+0+0 ---> [[ 1, 0, 0 ]]
						+1.. ---> [[ 1, 0, 0 ]]
		 */
		r_mat vec(3,1);
		vec(0,0)=0;
		vec(1,0)=0;
		vec(2,0)=0;
		
		// Analyze the bnd component and split it into : 1. arrayNxyzStr, 2. arraySignStr.
		auto tmpBnd = bnd;
		replaceAll(tmpBnd, "+", " ");
		replaceAll(tmpBnd, "-", " ");
		auto arrayNxyzStr = split(tmpBnd, " ");
		vector<char> arraySignStr;
		for (int i=0 ; i<bnd.size(); i++) {
			if (bnd[i] == '+' or bnd[i] == '-') {
				arraySignStr.push_back(bnd[i]);
			}
		}
		
		if (arrayNxyzStr.size() == arraySignStr.size()) { // Test for valid formate of 'bnd'.
		
			// Construct the bond_vector
			for (unsigned i=0 ; i<arraySignStr.size(); i++) {
				if ( i< bondVector.size() ){
					switch (arraySignStr[i]) {
						case '+':
							vec = vec + bondVector[i].bndVec;
							break;
						case '-':
							vec = vec - bondVector[i].bndVec;
							break;
					}
				} else {
					string errorMsg = "Error, cannot convert the formate of bond vector: "+bnd;
					cout<<errorMsg<<endl;
					throw errorMsg;
				}
			}
		} else {
			string errorMsg = "Error, cannot convert the formate of bond vector: "+bnd;
			cout<<errorMsg<<endl;
			throw errorMsg;
		}
		
		return vec;
	}
	vector<Atom>	site_iteration	()							{ return lattice_atom; }
	vector<AtomPair>pair_iteration	()							{ return vec_pair; }
	vector<string>	get_bonding		()							{return bonding;}
	vector<string>	get_bond_basis	()							{return bond_basis;}
	
	r_mat			basis_vec()									{
		r_mat vec(3,3);
		
		vector<int> N;
		N.push_back(N1);
		N.push_back(N2);
		N.push_back(N3);
		for (unsigned i=0; i<basis_vector.size(); i++) {
			istringstream iss(basis_vector[i]);
			string name;
			double x,y,z;
			iss>>name>>x>>y>>z;
			
			vec(i,0)=N[i]*x;
			vec(i,1)=N[i]*y;
			vec(i,2)=N[i]*z;
		}
	
		return vec;
	}
	r_mat			get_reciprocal()							{
		
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
	
	vector<string>		bond_basis_str;
	map<string, r_mat>	bond_basis;
	
	r_mat	get_vec(string bond_label){
		r_mat ret_vec(1,3);
		
		for (unsigned ii=0; ii<bond_basis_str.size(); ii++) {
			unsigned p_bond_count = WordCount(bond_label, "+"+bond_basis_str[ii]);
			unsigned m_bond_count = WordCount(bond_label, "-"+bond_basis_str[ii]);
			
			ret_vec = ret_vec + p_bond_count*1.0*bond_basis[bond_basis_str[ii]];
			ret_vec = ret_vec - m_bond_count*1.0*bond_basis[bond_basis_str[ii]];
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
		
		auto bbasis = Lat.get_bond_basis();
		
		for (unsigned i =0 ; i<bbasis.size(); i++)			{
			auto word = split(bbasis[i], " "); /*
							Split each line lf the bond_basis[i] = "x	0.5		0.0		0.0"
						*/
			
			r_mat	b_vec(1,3);
			auto b_name = word[0];
			b_vec[0] = StrToDouble(word[1]);
			b_vec[1] = StrToDouble(word[2]);
			b_vec[2] = StrToDouble(word[3]);

			bond_basis[b_name] = b_vec;
			bond_basis_str.push_back(b_name);
		}
		
		// Construct the expended_bond_basis for the full pack 3D bond-label and vec.
		vector<deque<string> > expended_bond_basis(bond_basis_str.size());
		for (unsigned ii=0; ii<bond_basis_str.size(); ii++) {
			
			string m_bnd="-"+bond_basis_str[ii]; // The minus sign for the bond basis.
			string m_bnd_acc="";
			
			string p_bnd="+"+bond_basis_str[ii]; // The plus  sign for the bond basis.
			string p_bnd_acc="";
			
			for (unsigned i_lev=0; i_lev<max_bond_level; i_lev++){
				m_bnd_acc+=m_bnd+",";
				expended_bond_basis[ii].push_front(m_bnd_acc);
			}
			
			expended_bond_basis[ii].push_back("o,");

			for (unsigned i_lev=0; i_lev<max_bond_level; i_lev++){
				p_bnd_acc+=p_bnd+",";
				expended_bond_basis[ii].push_back(p_bnd_acc);
			}
		}

		// Construct a full pack of 3D bond_label with vectors (that is : atom_bond)
		for (unsigned i0=0; i0<expended_bond_basis[0].size(); i0++)
		for (unsigned i1=0; i1<expended_bond_basis[1].size(); i1++)
		for (unsigned i2=0; i2<expended_bond_basis[2].size(); i2++) {
			auto total_bond_label = expended_bond_basis[0][i0]+expended_bond_basis[1][i1]+expended_bond_basis[2][i2];
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






















