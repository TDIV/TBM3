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
//  atom.hpp
//  TBM^3
//

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

enum H_SPACE { NORMAL, NAMBU, EXNAMBU };

//----------------------
class Atom{
private:
	map<string, unsigned>	orbitalMap;
	deque<string>			orbitalLabel;
	
	map<string, size_t>		indexMap;	//The map of the orbital index.
	vector<pair<string, unsigned> > indexList;
	
	void	appendIndexLabel(string key){
		indexMap[key] = indexList.size();
		indexList.push_back(make_pair(key,totalIndexSize));
		totalIndexSize++;
	}

public:
	string				spinDegree;
	H_SPACE				h_space;

	unsigned			orbitalOriginIndex;
	unsigned			atomIndex;
	string				atomName;  //example: Fe
	r_mat				pos;
	
	static	unsigned	totalIndexSize;
	
	Atom()	{ }
	~Atom()	{
		orbitalMap.clear();
		orbitalLabel.clear();
		indexMap.clear();
		indexList.clear();
	}
	void	createOrbitalContent(deque<string> _orbitalLabel, string _spinDegree, H_SPACE _h_space){
		spinDegree = _spinDegree;
		h_space = _h_space;
		
		if( _orbitalLabel.size() > 0){
			atomName = _orbitalLabel[0];
			_orbitalLabel.pop_front();
			orbitalLabel = _orbitalLabel;
		}
		for( unsigned ii = 0; ii<orbitalLabel.size() ; ii++){
			orbitalMap[orbitalLabel[ii]] = ii;
		}
		
		for(auto & label: orbitalLabel){
			
			switch( _h_space){
			case NORMAL:
				if( spinDegree == "off"){
					appendIndexLabel(label+".n");
				}
				else if( spinDegree == "on"){
					appendIndexLabel(label+".u");
					appendIndexLabel(label+".d");
				}
				break;
					
			case NAMBU:
				if( spinDegree == "off"){
					appendIndexLabel(label+".A");
					appendIndexLabel(label+".B");
				}
				else if( spinDegree == "on"){
					appendIndexLabel(label+".Au");
					appendIndexLabel(label+".Bd");
				}
				break;
					
			case EXNAMBU:
				if( spinDegree == "off"){
					string ErrorMsg = "Error: A spinless atom cannot be assigned with 'Extended Nambu' space.";
					ErrorMessage(ErrorMsg);
				}
				else if( spinDegree == "on"){
					appendIndexLabel(label+".Au");
					appendIndexLabel(label+".Ad");
					appendIndexLabel(label+".Bu");
					appendIndexLabel(label+".Bd");
				}
				break;
			}
		}
	}
	
	unsigned			index(string strKey)		{
		string	indexKey = strKey;
		auto	firstParser = split(strKey, ".");
		
		if (firstParser.size() == 1) { // Means the user is using "1u", "1d" style output.
			auto secondParser = splitByIntNumber(strKey);
			if( secondParser.size() == 1){
				auto ii = StrToInt(secondParser[0])-1;
				auto orbLabel = orbitalLabel[ii];
				auto spaceLabel = "n";
				indexKey = orbLabel+"."+spaceLabel;
			}
			else if( secondParser.size() == 2){
				auto ii = StrToInt(secondParser[0])-1;
				auto orbLabel = orbitalLabel[ii];
				auto spaceLabel = secondParser[1];
				indexKey = orbLabel+"."+spaceLabel;
			}
		}
		
		if( indexMap.find(indexKey) == indexMap.end()){
			string ErrorMsg = "Error, : cannot find the orbital space key of:'"+strKey+"'.";
			ErrorMessage(ErrorMsg);
		}
		
		return indexList[ indexMap[indexKey] ].second;
	}
	bool				hasOrbital(string orb)		{
		if( IsIntStr(orb)){
			unsigned orbIndex = StrToInt(orb)-1;
			return ( orbIndex < orbitalLabel.size());
		}
		return (orbitalMap.find(orb) != orbitalMap.end());
	}
	string				getOrbitalNumber(string orb){
		if( IsIntStr(orb))	{ return orb; }
		else				{ return IntToStr(orbitalMap[orb]+1); }
	}
	
	vector<pair<string, unsigned> > &	allIndexList	()					{ return indexList; }
	vector<pair<string, unsigned> >		orbitalIndexList(string orb)		{
		vector<pair<string, unsigned> > subIndexList;
		if(!hasOrbital(orb)){ return subIndexList; }
		
		string orbStr = orb;
		if( IsIntStr(orb) ) orbStr = orbitalLabel[StrToInt(orb)-1];
		
		string indexKey = "";
		switch( h_space){
		case NORMAL:
			if		( spinDegree == "off"){
				subIndexList.push_back(make_pair("n", index(orbStr+".n")));
			}
			else if	( spinDegree == "on"){
				subIndexList.push_back(make_pair("u", index(orbStr+".u")));
				subIndexList.push_back(make_pair("d", index(orbStr+".d")));
			}
			break;
				
		case NAMBU:
			if		( spinDegree == "off"){
				subIndexList.push_back(make_pair("A", index(orbStr+".A")));
				subIndexList.push_back(make_pair("B", index(orbStr+".B")));
			}
			else if	( spinDegree == "on"){
				subIndexList.push_back(make_pair("Au", index(orbStr+".Au")));
				subIndexList.push_back(make_pair("Bd", index(orbStr+".Bd")));
			}
			break;
				
		case EXNAMBU:
			//if		( spinDegree == "off"){
			//}
			if	( spinDegree == "on"){
				subIndexList.push_back(make_pair("Au", index(orbStr+".Au")));
				subIndexList.push_back(make_pair("Ad", index(orbStr+".Ad")));
				subIndexList.push_back(make_pair("Bu", index(orbStr+".Bu")));
				subIndexList.push_back(make_pair("Bd", index(orbStr+".Bd")));
			}
			break;
		}
		
		return subIndexList;
	}
	vector<pair<string, unsigned> >		spinIndexList	(string orbSpin)	{
		vector<pair<string, unsigned> > subIndexList;
		
		auto parser = split(orbSpin, ".");
		if( parser.size() != 2) parser = splitByIntNumber(orbSpin);
		if( parser.size() != 2) {
			ErrorMessage("Error, wrong spin-orbital formate:"+orbSpin);
		}
		
		if(!hasOrbital(parser[0])){ return subIndexList; }
		
		string orbStr = parser[0];
		if( IsIntStr(orbStr) ) orbStr = orbitalLabel[StrToInt(orbStr)-1];
		
		string indexKey = "";
		switch( h_space){
		case NORMAL:
			subIndexList.push_back(make_pair(parser[1], index(orbStr+"."+parser[1])));
			break;
				
		case NAMBU:
			if( parser[1] == "u") subIndexList.push_back(make_pair("A", index(orbStr+".Au")));
			if( parser[1] == "d") subIndexList.push_back(make_pair("B", index(orbStr+".Bd")));
			break;
				
		case EXNAMBU:
			if( parser[1] == "u"){
				subIndexList.push_back(make_pair("A", index(orbStr+".Au")));
				subIndexList.push_back(make_pair("B", index(orbStr+".Bu")));
			}
			if( parser[1] == "d"){
				subIndexList.push_back(make_pair("A", index(orbStr+".Ad")));
				subIndexList.push_back(make_pair("B", index(orbStr+".Bd")));
			}
			break;
		}
		return subIndexList;
	}
	
	void				printIndexLabel()			{
		for(unsigned i=0 ; i<indexList.size() ; i++){
			cout<<indexList[i].first<<endl;
		}
	}
	Atom &				operator=(const Atom &at)	{
		orbitalOriginIndex = at.orbitalOriginIndex;
		orbitalMap		= at.orbitalMap;
		orbitalLabel	= at.orbitalLabel;
		indexMap		= at.indexMap;
		indexList		= at.indexList;
		
		spinDegree		= at.spinDegree;
		h_space			= at.h_space;
		atomIndex		= at.atomIndex;
		atomName		= at.atomName;
		pos				= at.pos;
		return *this;
	}
	string				posToStr()					{
		/*
		 Return a string that contain the orbitalOriginIndex and position of the atom, example:
		 1		1.0		1.0		1.0
		 ...
		 */
		string retStr = "";
		retStr += fformat(orbitalOriginIndex,4)+" ";
		for(unsigned i=0 ; i<pos.size() ; i++){
			retStr += fformat(pos[i],15)+" ";
		}
		return retStr;
	}
};

unsigned Atom::totalIndexSize = 0;

class AtomPair{
public:
	Atom atomI;
	Atom atomJ;
	r_mat inputBond;
	
	AtomPair(){ }
	AtomPair(Atom _atomI, Atom _atomJ, r_mat _inputBond):inputBond(_inputBond), atomI(_atomI), atomJ(_atomJ){
	}
	r_mat bondIJ(){
		r_mat bondDelta(1,3);
		for( unsigned i=0 ; i<bondDelta.size() ; i++)
			bondDelta[i] = atomJ.pos[i] - atomI.pos[i];
		return bondDelta;
	}
	bool withinRange(r_var valid_radius = 0.001){
		r_mat distance(1,3);
		r_mat outputBond = bondIJ();
		distance = outputBond - inputBond;
		
		return sqrt(cdot(distance, distance)) < valid_radius;
	}
	
	AtomPair & operator= (const AtomPair & ap){
		if( this != &ap){
			atomI = ap.atomI;
			atomJ = ap.atomJ;
			inputBond = ap.inputBond;
		}
		return *this;
	}
};










