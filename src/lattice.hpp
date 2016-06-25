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
//  lattice.hpp
//  TBM^3
//

// Declare the RTree object.
typedef bg::model::point<r_var, 3, bg::cs::cartesian> point;
typedef point bondVec;
typedef std::pair<point, unsigned> value;
typedef bg::model::box<point> box;


//----------------------
class Lattice{
	unsigned index_size;
public:
	Lattice(){}
	Lattice(string _filename):atomIndex(-1),  filename(_filename){
		open(filename);
		createAtomList();
	}
	Lattice(string _filename, string spin, string space):atomIndex(-1),  filename(_filename){
		open(filename);
		parameter.STR("spin") = spin;
		parameter.STR("space") = space;
		createAtomList();
	}

	// Read from xxx.lat
	Parameter			parameter;
	BasisVector			basisVector;
	BondVector			bondVector;
	OrbitalProfile		orbitalProfile;
	AtomStringParser	atomParser;
	KSymmetryPoint		kSymmPointParser;
	
	// Read from xxx.lat.tbm
	InitOrder			initOrder;
	HamiltonianParser	hamParser;
	CoreCharge			coreCharge;
	
	void			open(string _filename)	{
		parameter.clear();
		basisVector.clear();
		orbitalProfile.clear();
		atomParser.clear();
		
		extendedAtomList.clear();
		rtree.clear();
		
        ifstream infile(_filename);
		if ( infile.is_open() ) {
			string flag = "";
			string sub_flag = "";
            string line;
			
            while ( getline(infile, line) ) {
				deleteComment(line); // Clean the commented words
                istringstream iss(line);
				
				auto lineParser = split(line, " ");
				
				if (lineParser.size() > 0){
					if (lineParser[0] == parameter())			{flag = lineParser[0]; continue;}
					if (lineParser[0] == basisVector())			{flag = lineParser[0]; continue;}
					if (lineParser[0] == orbitalProfile())		{flag = lineParser[0]; continue;}
					if (lineParser[0] == atomParser())			{flag = lineParser[0]; continue;}
					if (lineParser[0] == kSymmPointParser())	{flag = lineParser[0]; continue;}
					if (lineParser[0] == bondVector())		{
						flag = lineParser[0];
						if( lineParser.size() == 1){
							sub_flag = "0";
						}
						if( lineParser.size() >= 2){
							sub_flag = lineParser[1];
						}
						continue;
					}
			
					if (flag == parameter())		{ parameter.append(line);						continue;	}
					if (flag == basisVector())		{ basisVector.append(line);						continue;	}
					if (flag == orbitalProfile())	{ orbitalProfile.append(line);					continue;	}
					if (flag == atomParser())		{ atomParser.append(line);						continue;	}
					if (flag == kSymmPointParser())	{ kSymmPointParser.append(line);				continue;	}
					if (flag == bondVector())		{ bondVector.append(StrToInt(sub_flag) ,line);	continue;	}
				}
			}
		}
		infile.close();
		
		infile.open(_filename+".tbm");
		if ( infile.is_open() ) {
			string line;
			
			string flag = "";
			string header = "";
			
            while ( getline(infile, line) ) {
				deleteComment(line); // Clean the commented words
				istringstream iss(line);
				iss >> header;
				if	( header == coreCharge())	{flag = header;	continue;}
				if	( header == initOrder())	{flag = header;	continue;}
				if	( header == hamParser())	{flag = header;	continue;}
				
				if	( flag == coreCharge())		{ coreCharge.append(line);	continue; }
				if	( flag == initOrder())		{ initOrder.append(line);	continue; }
				if	( flag == hamParser())		{ hamParser.append(line);	continue; }
			}
		}
		infile.close();
		
	}
	size_t			latticeSize()			{ return atomList.size(); }
	size_t			indexSize()				{ return index_size; }
	
	Atom			getAtom(unsigned index)					{
		if( index >= atomList.size() ){
			string ErrorMsg = "Error, index out of range.";
			ErrorMessage(ErrorMsg);
		}
		return atomList[index];
	}
	AtomPair		getPair(unsigned index, r_mat bond)		{
		Atom atomI = getAtom(index);
		r_mat posJ = atomI.pos + bond;
		point p = point(posJ[0], posJ[1], posJ[2]);
		
		vector<value>	result_s;
		rtree.query(bgi::nearest(p, 1), std::back_inserter(result_s));
		
		unsigned exIndex = result_s[0].second;
		Atom atomJ = extendedAtomList[exIndex];
		
		return AtomPair(atomI, atomJ, bond);
	}
	AtomPair		getPair(unsigned index, string bondStr)	{
		return getPair(index, vec(bondStr));
	}
	pair<Atom, vector<Atom> > getBox(unsigned index)		{
		Atom atomI = getAtom(index);
		r_mat pos = atomI.pos;
		
		vector<value> result_s;
		r_var radius = abs(parameter.VAR("bondRadius", 1).real());
		box query_box(point(pos[0] - radius, pos[1] - radius, pos[2] - radius),
					  point(pos[0] + radius, pos[1] + radius, pos[2] + radius));
		
		rtree.query(bgi::intersects(query_box), std::back_inserter(result_s));
		
		vector<Atom> atomJList;
		for(unsigned i=0 ; i<result_s.size(); i++){
			auto indexJ = result_s[i].second;
			Atom atomJ = extendedAtomList[indexJ];
			
			r_mat distanceIJ = atomJ.pos - atomI.pos;
			if( sqrt(cdot(distanceIJ, distanceIJ))<= radius )
				atomJList.push_back(atomJ);
		}
		return make_pair(atomI, atomJList);
	}
	
	Atom			getAtom()								{
		if( !(atomIndex >= 0 and atomIndex < atomList.size()) ){
			string ErrorMsg = "Error, index out of range.";
			ErrorMessage(ErrorMsg);
		}
		return atomList[atomIndex];
	}
	AtomPair		getPair(r_mat bond)						{
		Atom atomI = getAtom();
		r_mat posJ = atomI.pos + bond;
		point p = point(posJ[0], posJ[1], posJ[2]);
		
		vector<value>	result_s;
		rtree.query(bgi::nearest(p, 1), std::back_inserter(result_s));
		
		unsigned exIndex = result_s[0].second;
		Atom atomJ = extendedAtomList[exIndex];
		
		return AtomPair(atomI, atomJ, bond);
	}
	AtomPair		getPair(string bondStr)					{
		return getPair(vec(bondStr));
	}
	pair<Atom, vector<Atom> > getBox()						{
		Atom atomI = getAtom();
		r_mat pos = atomI.pos;
		
		vector<value> result_s;
		r_var radius = abs(parameter.VAR("bondRadius", 1).real());
		box query_box(point(pos[0] - radius, pos[1] - radius, pos[2] - radius),
					  point(pos[0] + radius, pos[1] + radius, pos[2] + radius));
		
		rtree.query(bgi::intersects(query_box), std::back_inserter(result_s));
		
		vector<Atom> atomJList;
		for(unsigned i=0 ; i<result_s.size(); i++){
			auto indexJ = result_s[i].second;
			Atom atomJ = extendedAtomList[indexJ];
			
			r_mat distanceIJ = atomJ.pos - atomI.pos;
			if( sqrt(cdot(distanceIJ, distanceIJ))<= radius )
				atomJList.push_back(atomJ);
		}
		return make_pair(atomI, atomJList);
	}
	
	r_mat			vec(string line)	{
		// line = "+1+0+1"		-> use Cartesian coordinate
		// line = "+1+0+1#"		-> use BondVector 0
		// line = "+1+0+1#0"	-> use BondVector 0
		// line = "+1+0+1#1"	-> use BondVector 1
		// line = "+1+0+1#2"	-> use BondVector 2
		r_mat retVec(1,3);
		
		auto lineParser = split(line, "#");
		
		int CoordinateSelect = -1;
		int countSharpSimbol = WordCount("#", line);
		if( countSharpSimbol == 1){
			if( lineParser.size() == 2)
				CoordinateSelect = StrToInt(lineParser[1]);
			else if(lineParser.size() == 1){
				CoordinateSelect = 0;
			}
		}
		
		auto vecParser = splitByIntNumber(lineParser[0]);
		if( vecParser.size() == 6){
			if( CoordinateSelect == -1){ // Cartesian coordinate
				retVec = tbm::vec(
							 PMStr(vecParser[0])*StrToDouble(vecParser[1]),
							 PMStr(vecParser[2])*StrToDouble(vecParser[3]),
							 PMStr(vecParser[4])*StrToDouble(vecParser[5])
							 );
			}
			else{
				auto bvec = bondVector.getBond(CoordinateSelect);
				retVec =PMStr(vecParser[0]) * StrToDouble(vecParser[1]) * bvec[0] +
						PMStr(vecParser[2]) * StrToDouble(vecParser[3]) * bvec[1] +
						PMStr(vecParser[4]) * StrToDouble(vecParser[5]) * bvec[2] ;
			}
		}
		
		return retVec;
	}
	H_SPACE			HSpace()			{return h_space;}
	string			FileName()			{return filename;}
	
	vector<Atom>	getAtomList(){ return atomList; }
	
	bool			iterate()			{
		atomIndex++;
		
		if( atomIndex == latticeSize())	{
			atomIndex = -1;
			return false;
		}
		return true;
	}
	
private:
	
	void createAtomList(){
		Atom::totalIndexSize = 0;
		index_size = 0;
		
		auto	spinDegree	= parameter.STR("spin", "on");
		auto	hSpace		= parameter.STR("space", "normal");
		auto	coordinate	= parameter.STR("coordinate", "cartesian");
		
		if		( hSpace == "normal")	{ h_space = NORMAL;}
		else if	( hSpace == "nambu")	{ h_space = NAMBU;}
		else if	( hSpace == "exnambu")	{ h_space = EXNAMBU;}
		
		auto	atomInfoList = atomParser.atomInfoList;
		auto	orbitalList	= orbitalProfile.getOrbitalList();
		
		// -----------------------------------------------------------------------
		// Generate the index of each atom.
		// -----------------------------------------------------------------------
		for( unsigned i=0 ; i<atomInfoList.size() ; i++){
			int			orbitalLabelIndex = atomInfoList[i].first - 1; // The info list starting with 1 (not 0).
			r_mat		atomPos = atomInfoList[i].second;
			
			if (orbitalLabelIndex >= 0){
				auto	orbitalLabel = orbitalList[orbitalLabelIndex];
				
				Atom	tmpAtom;
				tmpAtom.orbitalOriginIndex = orbitalLabelIndex + 1;
				tmpAtom.atomIndex = i;
				tmpAtom.pos = atomPos;
				tmpAtom.createOrbitalContent(orbitalLabel, spinDegree, h_space);
				atomList.push_back(tmpAtom);
			}
		}
		
		index_size = Atom::totalIndexSize;
		
		updateExtendedAtomList();

	}
	
	void updateExtendedAtomList(){
		
		extendedAtomList.clear();
		
		// -----------------------------------------------------------------------
		// Create the extendedAtomList from the atomList by a given bondRadius.
		// -----------------------------------------------------------------------
		auto bondRadius = parameter.VAR("bondRadius", 1);
		r_var N = basisVector.minRepeatForRadius(bondRadius.real());
		auto AVec = basisVector.getAVec();
		
		if		( AVec.size() == 3){
			for( r_var i = -N ; i<=N ; i+=1)
			for( r_var j = -N ; j<=N ; j+=1)
			for( r_var k = -N ; k<=N ; k+=1){
				for(unsigned index=0 ; index<atomList.size() ; index++){
					Atom tmpAtom;
					tmpAtom = atomList[index];
					tmpAtom.pos = tmpAtom.pos+ i*AVec[0] +j*AVec[1] +k*AVec[2];
					extendedAtomList.push_back(tmpAtom);
				}
			}
		}
		else if	( AVec.size() == 2){
			for( r_var i = -N ; i<=N ; i+=1)
			for( r_var j = -N ; j<=N ; j+=1){
				for(unsigned index=0 ; index<atomList.size() ; index++){
					Atom tmpAtom;
					tmpAtom = atomList[index];
					tmpAtom.pos = tmpAtom.pos+ i*AVec[0] +j*AVec[1];
					extendedAtomList.push_back(tmpAtom);
				}
			}
		}
		else if	( AVec.size() == 1){
			for( r_var i = -N ; i<=N ; i+=1){
				for(unsigned index=0 ; index<atomList.size() ; index++){
					Atom tmpAtom;
					tmpAtom = atomList[index];
					tmpAtom.pos = tmpAtom.pos+ i*AVec[0];
					extendedAtomList.push_back(tmpAtom);
				}
			}
		}
		
		// -----------------------------------------------------------------------
		// Create the rtree for atom position search.
		// -----------------------------------------------------------------------
		rtree.clear();
		for( unsigned exIndex=0 ; exIndex<extendedAtomList.size() ; exIndex++){
			auto pos =	extendedAtomList[exIndex].pos;
			point p = point( pos[0], pos[1], pos[2]);
			rtree.insert(std::make_pair(p, exIndex));
		}
	}

	bgi::rtree< value, bgi::quadratic<16> > rtree;
	int				atomIndex;
	vector<Atom>	atomList;
	vector<Atom>	extendedAtomList;
	
	string	filename;
	H_SPACE	h_space;
	
};




















