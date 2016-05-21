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
//  atom.cpp
//  TBM^3
//
//  Created by Yuan Yen Tai on 7/02/15.
//

/*
   The library is used to create a multi orbital lattice interface for different tight-binding models.
*/


enum H_SYMMETRY { NORMAL, NAMBU, EXNAMBU };

class Atom{
private:
    unsigned			unitcell;   //
    unsigned			atom_index; //
	string				atom_name;  //example: Fe-1
	string				atom_info;  //atom_info comes from the sub_atom,
									//example:  Fe   2s   0.5  0.5   0.5
	
	unsigned			orbitalNumber;
	string				spinDegree;
	
	pair<int, string>	getOrbitalInfo()	const	{
					
		auto	ss=split(atom_info, " ");
		string	degree_of_freedom = ss[1];
		auto	label = analyze_atom_label(degree_of_freedom);
		int		num_of_orbital	= StrToInt(label[0]);
		string	sub_degree		= label[1];
		return	make_pair(num_of_orbital, sub_degree);
	}
public:
    r_mat				pos;		//The atom position
	vector<string>		indexLabel;	//The labels for the index, example: 1u 1d 2u 2d
	map<string, int>	indexMap;	//The map of the orbital index, ex: {1u : 0}
	vector<pair<double,double> > LDOS;//To storage the temporarily calculated LDOS.

	Atom	()					{
		unitcell=0;
		atom_index=0;
		atom_name="";
		pos = r_mat(3,1).get();
		orbitalNumber = 0;
		spinDegree = "";
	}
	Atom	(const Atom &at)	{
		*this=at;
	}	// The copy constructor
	Atom	(string line)		{
        istringstream iss(line);
        double x,y,z;
        pos = r_mat(3,1).get();
        
        iss >> unitcell >> atom_index >> atom_name >> x >> y >> z;
        
        pos[0]=x;
        pos[1]=y;
        pos[2]=z;
    }	// The constructor from each Atom instance with a string input
	~Atom	()					{
		indexLabel.clear();
		indexMap.clear();
	}
	
	Atom&		operator=(const Atom &at)			{
		if ( this != &at ) {
			unitcell	=at.unitcell;
			atom_index	=at.atom_index;
			atom_name	=at.atom_name;
			pos			=at.pos;
			atom_info	=at.atom_info;
			indexLabel	=at.indexLabel;
			indexMap	=at.indexMap;
			LDOS		=at.LDOS;
			orbitalNumber = at.orbitalNumber;
			spinDegree	= at.spinDegree;
		}
		return *this;
	}// The copy operation
	
	string		Info		()				const	{
		ostringstream oss;
		oss<<endl;
		oss<<"Atom Name  : "<< Name() <<endl;
		oss<<"Atom Index : "<< atom_index <<endl;
		oss<<"Labels     : ";
		for (unsigned i=0; i<indexLabel.size(); i++) {
			oss<<indexLabel[i]<<" ";
		}
		oss<<endl;
		oss<<"Position   : ("<<pos.const_index(0)<<", "<<pos.const_index(1)<<", "<<pos.const_index(2)<<")"<<endl;
		oss<<"In unitcell: "<<unitcell<<endl;
		return oss.str();
	}
	string		Name		()				const	{
		auto ss=split(atom_name,"-");
		return ss[0];
	}// Get the Atom name from "atom_name"
	string		subName		()				const	{
		return atom_name;
	}// Get the Atom SubName from "atom_info"
	int			subIndex	()				const	{
		auto ss=split(atom_name,"-");
		return StrToInt(ss[1]);
	}// Get the Atom sub_index from "atom_name"
	unsigned	unitcellIndex()				const	{
		return unitcell;
	}
	bool		hasIndex	(string ss)		const	{
		auto it=indexMap.find(ss);
		if (it!=indexMap.end()) return true;
		
		return false;
	}
	
	void		setAtomInfo	(string ss)				{
		atom_info = ss;
		auto info = getOrbitalInfo();
		orbitalNumber = info.first;
		spinDegree = info.second;
	}
	
	unsigned	atomIndex	()				const	{ return atom_index; }
	unsigned	getOrbitalNumber()					{ return orbitalNumber; }
	string		getSpinDegree()						{ return spinDegree; }

	int			operator[]	(string label)	const	{
		if (hasIndex(label)) { return indexMap.find(label)->second; }
		return -1;
	}
	
	void createIndexLabel(H_SYMMETRY symmetry = NORMAL){
		auto ss=split(atom_info, " ");
		
		indexLabel.clear();
		
		if (ss.size()>2) {
			
			int		num_of_orbital	= getOrbitalNumber();
			string	sub_degree		= getSpinDegree();
			
			vector<string>	sub_degree_label;
			
			switch( symmetry ) {
			case NORMAL:
				if (sub_degree == "") {		// Normal space
					sub_degree_label.push_back("");
				}
				else if (sub_degree == "s"){// Spin space
					sub_degree_label.push_back("u");
					sub_degree_label.push_back("d");
				}
				break;
				
			case NAMBU:
				if (sub_degree == "N"){		// Nambu space without spin
					sub_degree_label.push_back("A");
					sub_degree_label.push_back("B");
				}
				else if (sub_degree == "Ns"){// Nambu space with spin
					sub_degree_label.push_back("Au");
					sub_degree_label.push_back("Bd");
				}
				break;
				
			case EXNAMBU:
				if (sub_degree == "") {
					cout<<endl;
					cout<<"Error: A spinless atom cannot be assigned with 'Extended Nambu' space."<<endl;
					cout<<"Program stoped!"<<endl;
					throw "Program stop";
				}
				else if (sub_degree == "Es"){	// Extended Nambu space with spin
					sub_degree_label.push_back("Au");
					sub_degree_label.push_back("Ad");
					sub_degree_label.push_back("Bu");
					sub_degree_label.push_back("Bd");
				}
				break;
			}
			
			// Create the index_label  ... ex: for a 2-orbital system
			//		Normal case	          : 1,  2
			//		Nambu wo spin         : 1A, 1B, 2A, 2B
			//		Spin case	          : 1u, 1d, 2u, 2d
			//		Nambu wi spin         : 1Au,1Bd,2Au,2Bd
			//		Extended Nambu wi spin: 1Au,1Ad,1Bu,1Bd,2Au,2Ad,2Bu,2Bd
			for (int i=0; i<num_of_orbital; i++) {
				for (int j=0; j<sub_degree_label.size(); j++){
					ostringstream oss;
					oss<<i+1<<sub_degree_label[j];
					indexLabel.push_back(oss.str());
				}
			}
		}
	} // Create "index_label"
};

class AtomPair{
public:
	Atom	AtomI;
	Atom	AtomJ;
	string	bond;
	r_mat	bvec;
	
	AtomPair(){
		AtomI	= Atom();
		AtomJ	= Atom();
		bond	= "";
		bvec	= r_mat();
	}

	AtomPair(const AtomPair& ap){
		AtomI	= ap.AtomI;
		AtomJ	= ap.AtomJ;
		bond	= ap.bond;
		bvec	= ap.bvec;
	}
	
	AtomPair& operator=(const AtomPair &ap){
		if ( this != &ap ) {
			AtomI	= ap.AtomI;
			AtomJ	= ap.AtomJ;
			bond	= ap.bond;
			bvec	= ap.bvec;
		}
		return *this;
	}
};

