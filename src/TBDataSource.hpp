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
//  TBDataSource.hpp
//  TBM^3
//

struct	MatrixElement{
public:
	int		I,J;
	x_var	val;
	r_mat	bondVector;
	
	MatrixElement(){ }
	MatrixElement(int i, int j, x_var _val, r_mat _bondVector){
		I=i;
		J=j;
		val=_val;
		bondVector = _bondVector;
	}
	MatrixElement& operator()(int i, int j, x_var _val, r_mat _bondVector){
		I=i;
		J=j;
		val=_val;
		bondVector = _bondVector;
		return *this;
	}
	
};

struct	EigVec{
	string	message;
	r_mat	kPoint;
	r_mat	eigenValue;
	x_mat	eigenVector;
};

/*
 This structure will be passed through all future models.
 It can be used to construct all matrix element, for example:
 	orbital		> Fe 1		> 1		%Create the on-site energy for the first orbital of Fe.
 	hundSpin	> Fe 1		> 1,0,0 %Create a classical spin (x-direction) coupled to orbital 1.
 
	site		> Fe 1d:2u	> 1	%Create the on-site correlation for orbital-1 spin-dn and orbital-2 spin-up.
	siteHc		> Fe 1d:2u	> 1	% Hermitian conjugate operation.
 
	hopping		> Fe:Fe:+1+0+0	1:2	> 1		% Hopping terms.
	hoppingHc	> Fe:Fe:+1+0+0	1:2	> 1		% Hermitian conjgate of Hopping terms.
	bond		> Fe:Fe:+1+0+0	1u:2d	>	1	%
	bondHc		> Fe:Fe:+1+0+0	1u:2d	>	1	%
 
	pairingS	> Fe 1:1			> 1		% On-site Singlet pairing
	pairingS	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site Singlet pairing
	pairingU	> Fe 1:1			> 1		% On-site Triplet pairing for 'up' spin
	pairingU	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site Triplet pairing for 'up' spin
	pairingD	> Fe 1:1			> 1		% On-site Triplet pairing for 'dn' spin
	pairingD	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site Triplet pairing for 'dn' spin
	pairingT	> Fe 1:1			> 1		% On-site Triplet pairing
	pairingT	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site Triplet pairing
 */
class TBDataSource{
public:
	/////////////////////////////////////////////////
	//   The constructor                            /
	/////////////////////////////////////////////////
	TBDataSource(Lattice & _lat):Lat(_lat), order(_lat), order_old(_lat)	{	}
	~TBDataSource(){
		KEigenValVec.clear();
	}
	
	Lattice &				Lat;
	OrderParameter			order;
	OrderParameter			order_old;
	
	vector<MatrixElement>	hamElementList;	// A list to store all the matrix elements
	x_mat					Ham;			// The body of the Hamiltonian (matrix).
	vector<EigVec>			KEigenValVec;	// Calculated K-space dependent eigen value and vectore.
	//EigVec					tmpEigValVec;	// A temporary storage of eigen value vector for a specific K-point.
	r_var					maxE,	minE;
	
	map<string, double>		energyMap;

	void initHam	()		{
		Ham = x_mat(Lat.indexSize(), Lat.indexSize());
		initOrder();
	}
	void initOrder	()		{
		order.clear();
		order.load();
	}

	/*------------------------------------------------
	 Using these methods to construct the "hamElementList".
	 -------------------------------------------------*/

	void		addHundSpin		(string opt, deque<string> optList, deque<string> varList)	{
		/* Site operation: "Fe 1".	(operate for only a single orbital) */
		
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage("Error, 'HundSpin' operation:\n"+opt+"\nCannot be applied for a spin-independent Hamiltonian.");
		}
		if( optList.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto atomI = Lat.getAtom();
		if( atomI.atomName != optList[0]){ return; } // Name does not match.
		if(!atomI.hasOrbital(optList[1])){ return; } // No such orbital.
		
		r_var Jh = 1.0;
		x_mat xvec = parseSiteString(atomI, varList[0]);
		for( unsigned i=1 ; i<varList.size() ; i++){
			Jh = Jh * parseSiteString(atomI, varList[i])[0].real();
		}
		
		normalizeXVec(xvec);
		x_mat sVec = xvec * Jh;
	
		if( sVec.size() != 3 ){
			ErrorMessage("Error, operation:\n"+opt+".\n Has been applied with a wrong spin-variable with size:"+IntToStr(sVec.size()));
		}
		
		string orbital = optList[1];
		r_var Sx =-sVec[0].real();
		r_var Sy =-sVec[1].real();
		r_var Sz =-sVec[2].real();
		switch( Lat.HSpace() ){
			case NORMAL: {
				unsigned indexNu = atomI.index(orbital+"u");
				unsigned indexNd = atomI.index(orbital+"d");
				
				hamElementList.push_back(MatrixElement(indexNu, indexNd, Sx-Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexNd, indexNu, Sx+Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexNu, indexNu, Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexNd, indexNd,-Sz			, vec(0,0,0)));
				break;
			}
			case NAMBU: {
				unsigned indexAu = atomI.index(orbital+"Au");
				unsigned indexBd = atomI.index(orbital+"Bd");
				hamElementList.push_back(MatrixElement(indexAu, indexAu, Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBd, Sz			, vec(0,0,0)));
				if( Sx > 0.0001 or Sy > 0.0001){
					ErrorMessage("Error, spin-flip operation:\n"+
							 opt+
							 ".\n Cannot be applied for a spin-depende Nambu space Hamiltonian.");
				}
				break;
			}
			case EXNAMBU: {
				unsigned indexAu = atomI.index(orbital+"Au");
				unsigned indexAd = atomI.index(orbital+"Ad");
				unsigned indexBu = atomI.index(orbital+"Bu");
				unsigned indexBd = atomI.index(orbital+"Bd");
				
				hamElementList.push_back(MatrixElement(indexAu, indexAd, Sx-Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexAd, indexAu, Sx+Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexAu, indexAu, Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexAd, indexAd,-Sz			, vec(0,0,0)));
				
				hamElementList.push_back(MatrixElement(indexBu, indexBd,-Sx+Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBu,-Sx-Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBu, indexBu,-Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBd, Sz			, vec(0,0,0)));
				break;
			}
		}
	}
	void		addOrbitalEng	(string opt, deque<string> optList, deque<string> varList)	{
		/*
		 Site operation: "Fe 1".	(operate for only a single orbital)
		 */
		if( optList.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto atomI = Lat.getAtom();
		if( atomI.atomName != optList[0]){ return; } // Name does not match.
		if(!atomI.hasOrbital(optList[1])){ return; } // No such orbital.
		
		x_mat sVec = parseSiteString(atomI, varList[0]);
		r_var val = sVec[0].real();
		for( unsigned i=1 ; i<varList.size() ; i++){
			val = val * parseSiteString(atomI, varList[i])[0].real();
		}
		
		auto subIndex = atomI.orbitalIndexList(optList[1]);
		for(auto & iter: subIndex){
			if(iter.first[0] == 'n' or iter.first[0] == 'u' or iter.first[0] == 'd'){
				hamElementList.push_back(MatrixElement(iter.second, iter.second, val, vec(0,0,0)));
			}
			if(iter.first[0] == 'A'){
				hamElementList.push_back(MatrixElement(iter.second, iter.second, val, vec(0,0,0)));
			}
			if(iter.first[0] == 'B'){
				hamElementList.push_back(MatrixElement(iter.second, iter.second,-val, vec(0,0,0)));
			}
		}
	}

	void		addSiteCouple	(string opt, deque<string> optList, deque<string> varList)	{
		/*
		 Site operation: "Fe 1u:1u".	(operate for only a single orbital)
		 Site operation: "Fe 1u:2d".	(operate between orbitals and spin)
		 */
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage("Error, operation:\n"+opt+".\n Cannot be applied for a spin-independent Hamiltonian.");
		}
		if( optList.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(optList[0], ":");
		auto secondSec = split(optList[1], ":");
		if( firstSec.size() != 1 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		if( secondSec[0][secondSec[0].size()-1] != secondSec[1][secondSec[1].size()-1] and
		   Lat.HSpace() == NAMBU) {
			ErrorMessage("Error, operation:\n"+opt+".\n Spin-flip term cannot be applied under 'nambu' space.\n Please set the 'space' to 'normal' or 'exnambu'.");
		}
		
		auto atomI = Lat.getAtom();
		if( optList[0] != atomI.atomName )	return;
		
		x_var val = parseSiteString(atomI, varList[0])[0];
		for( unsigned i=1 ; i<varList.size() ; i++){
			val = val * parseSiteString(atomI, varList[i])[0];
		}
		
		auto spinIndexList_i = atomI.spinIndexList(secondSec[0]);
		auto spinIndexList_j = atomI.spinIndexList(secondSec[1]);
		
		if(spinIndexList_i.size() <= 0)
			ErrorMessage("Error, in the operation: \n'"+opt+
						 "'\n atom:'"+atomI.atomName+"' don't have any orbital to operate.");
		if(spinIndexList_j.size() <= 0)
			ErrorMessage("Error, in the operation: \n'"+opt+
						 "'\n atom:'"+atomI.atomName+"' don't have any orbital to operate.");
		
		if(spinIndexList_i.size() != spinIndexList_j.size()){
			ErrorMessage("Error, operation:\n"+opt+".\n the 'i' and 'j' spin-operation-list does not match eatch other.");
		}
		
		for( unsigned ii=0 ; ii<spinIndexList_i.size() ; ii++){
			auto & label_i = spinIndexList_i[ii].first;
			auto & index_i = spinIndexList_i[ii].second;
			auto & index_j = spinIndexList_j[ii].second;
			
			if( label_i == "n" or label_i == "u" or label_i == "d" or label_i == "A")
				hamElementList.push_back(MatrixElement(index_i, index_j, val, vec(0,0,0)));
			if( label_i == "B")
				hamElementList.push_back(MatrixElement(index_i, index_j,-val, vec(0,0,0)));
		}
	}
	void		addSiteCoupleHc	(string opt, deque<string> optList, deque<string> varList)	{
		/*
		 Site operation: "Fe 1u:1u".	(operate for only a single orbital)
		 Site operation: "Fe 1u:2d".	(operate between orbitals and spin)
		 */
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage("Error, operation:\n"+opt+".\n Cannot be applied for a spin-independent Hamiltonian.");
		}
		if( optList.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(optList[0], ":");
		auto secondSec = split(optList[1], ":");
		if( firstSec.size() != 1 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		if( secondSec[0][secondSec[0].size()-1] != secondSec[1][secondSec[1].size()-1] and
		   Lat.HSpace() == NAMBU) {
			ErrorMessage("Error, operation:\n"+opt+".\n Spin-flip term cannot be applied under 'nambu' space.\n Please set the 'space' to 'normal' or 'exnambu'.");
		}
		
		auto atomI = Lat.getAtom();
		if( optList[0] != atomI.atomName )	return;
		
		x_mat sVec = parseSiteString(atomI, varList[0]);
		x_var val = sVec[0];
		for( unsigned i=1 ; i<varList.size() ; i++){
			val = val * parseSiteString(atomI, varList[i])[0];
		}
		
		auto spinIndexList_i = atomI.spinIndexList(secondSec[0]);
		auto spinIndexList_j = atomI.spinIndexList(secondSec[1]);
		
		if(spinIndexList_i.size() <= 0)
			ErrorMessage("Error, in the operation: \n'"+opt+
						 "'\n atom:'"+atomI.atomName+"' don't have any orbital to operate.");
		if(spinIndexList_j.size() <= 0)
			ErrorMessage("Error, in the operation: \n'"+opt+
						 "'\n atom:'"+atomI.atomName+"' don't have any orbital to operate.");
		
		if(spinIndexList_i.size() != spinIndexList_j.size()){
			ErrorMessage("Error, operation:\n"+opt+".\n the 'i' and 'j' spin-operation-list does not match eatch other.");
		}
		
		for( unsigned ii=0 ; ii<spinIndexList_i.size() ; ii++){
			auto & label_i = spinIndexList_i[ii].first;
			auto & index_i = spinIndexList_i[ii].second;
			auto & index_j = spinIndexList_j[ii].second;
			
			if( label_i == "n" or label_i == "u" or label_i == "d" or label_i == "A"){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(index_j, index_i,conj(val), vec(0,0,0)));
			}
			if( label_i == "B"){
				hamElementList.push_back(MatrixElement(index_i, index_j,-val, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(index_j, index_i,-conj(val), vec(0,0,0)));
			}
		}
	}
	
	void		addHoppingInt	(string opt, deque<string> optList, deque<string> varList)	{
		/*
		 Bond operation: "Fe:O:+1+0+0 1:2".	(operate for only a single orbital)
		 */
		if( optList.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(optList[0], ":");
		auto secondSec = split(optList[1], ":");
		if( firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto pair = Lat.getPair(firstSec[2]);
		auto bond = pair.bondIJ();
		if( pair.atomI.atomName != firstSec[0])	return;
		if( pair.atomJ.atomName != firstSec[1])	return;
		if( !pair.withinRange() )				return;
		
		if(!pair.atomI.hasOrbital(secondSec[0]))return;
		if(!pair.atomJ.hasOrbital(secondSec[1]))return;
		
		x_mat sVec = parseBondString(pair, varList[0]);
		x_var val = sVec[0];
		for( unsigned i=1 ; i<varList.size() ; i++){
			val = val * parseBondString(pair, varList[i])[0];
		}
		
		auto subIndexI = pair.atomI.orbitalIndexList(secondSec[0]);
		auto subIndexJ = pair.atomJ.orbitalIndexList(secondSec[1]);
		for( unsigned ii = 0; ii<subIndexI.size() ; ii++){
			auto & iterI = subIndexI[ii];
			auto & iterJ = subIndexJ[ii];
			char & tmpChar = iterI.first[0];
			if(tmpChar == 'n' or tmpChar == 'u' or tmpChar == 'd'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));
			}
			if(tmpChar == 'A'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));
			}
			if(tmpChar == 'B'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,-val, pair.bondIJ()));
			}
		}
	}
	void		addHoppingIntHc	(string opt, deque<string> optList, deque<string> varList)	{
		/*
		 Bond operation: "Fe:O:+1+0+0 1:2".	(operate for only a single orbital)
		 */
		if( optList.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(optList[0], ":");
		auto secondSec = split(optList[1], ":");
		if( firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto pair = Lat.getPair(firstSec[2]);
		if( pair.atomI.atomName != firstSec[0])	return;
		if( pair.atomJ.atomName != firstSec[1])	return;
		if( !pair.withinRange() )				return;
		
		if(!pair.atomI.hasOrbital(secondSec[0]))return;
		if(!pair.atomJ.hasOrbital(secondSec[1]))return;
		
		x_mat sVec = parseBondString(pair, varList[0]);
		x_var val = sVec[0];
		for( unsigned i=1 ; i<varList.size() ; i++){
			val = val * parseBondString(pair, varList[i])[0];
		}
		
		auto subIndexI = pair.atomI.orbitalIndexList(secondSec[0]);
		auto subIndexJ = pair.atomJ.orbitalIndexList(secondSec[1]);
		for( unsigned ii = 0; ii<subIndexI.size() ; ii++){
			auto & iterI = subIndexI[ii];
			auto & iterJ = subIndexJ[ii];
			char & tmpChar = iterI.first[0];
			
			if(tmpChar == 'n' or tmpChar == 'u' or tmpChar == 'd'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val,			pair.bondIJ()));
				hamElementList.push_back(MatrixElement(iterJ.second, iterI.second, conj(val),	-pair.bondIJ()));
			}
			if(tmpChar == 'A'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val,			pair.bondIJ()));
				hamElementList.push_back(MatrixElement(iterJ.second, iterI.second, conj(val),	-pair.bondIJ()));
			}
			if(tmpChar == 'B'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,-val,			pair.bondIJ()));
				hamElementList.push_back(MatrixElement(iterJ.second, iterI.second,conj(-val),	-pair.bondIJ()));
			}
		}
	}
	
	void		addBondCouple	(string opt, deque<string> optList, deque<string> varList)	{
		/*
		 Bond operation: "Fe:Fe:+1+0+0 1u:1u".	(operate for only a single orbital)
		 Bond operation: "Fe:O:+1+0+0# 1u:2d".	(operate between orbitals and spin)
		 */
		if( optList.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(optList[0], ":");
		auto secondSec = split(optList[1], ":");
		if( firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto pair = Lat.getPair(firstSec[2]);
		if( pair.atomI.atomName != firstSec[0])	return;
		if( pair.atomJ.atomName != firstSec[1])	return;
		if( !pair.withinRange() )				return;
		
		//if(!pair.atomI.hasOrbital(secondSec[0]))return;
		//if(!pair.atomJ.hasOrbital(secondSec[1]))return;
		
		x_mat sVec = parseBondString(pair, varList[0]);
		x_var val = sVec[0];
		for( unsigned i=1 ; i<varList.size() ; i++){
			val = val * parseBondString(pair, varList[i])[0];
		}
		
		auto spinIndexList_i = pair.atomI.spinIndexList(secondSec[0]);
		auto spinIndexList_j = pair.atomJ.spinIndexList(secondSec[1]);
		
		if(spinIndexList_i.size() <= 0)
			ErrorMessage("Error, in the operation: \n'"+opt+
						 "'\n atom:'"+pair.atomI.atomName+"' don't have any orbital to operate.");
		if(spinIndexList_j.size() <= 0)
			ErrorMessage("Error, in the operation: \n'"+opt+
						 "'\n atom:'"+pair.atomJ.atomName+"' don't have any orbital to operate.");
		
		if(spinIndexList_i.size() != spinIndexList_j.size()){
			ErrorMessage("Error, operation:\n"+opt+".\n The 'i' and 'j' spin-operation-list does not match eatch other.");
		}
		
		for( unsigned ii=0 ; ii<spinIndexList_i.size() ; ii++){
			auto & tmpChar = spinIndexList_i[ii].first[0];
			
			auto & index_i = spinIndexList_i[ii].second;
			auto & index_j = spinIndexList_j[ii].second;
			
			if(tmpChar == 'n' or tmpChar == 'u' or tmpChar == 'd'){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, pair.bondIJ()));
			}
			if(tmpChar == 'A'){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, pair.bondIJ()));
			}
			if(tmpChar == 'B'){
				hamElementList.push_back(MatrixElement(index_i, index_j,-val, pair.bondIJ()));
			}
		}
	}
	void		addBondCoupleHc	(string opt, deque<string> optList, deque<string> varList)	{
		/*
		 Bond operation: "Fe:Fe:+1+0+0 1u:1u".	(operate for only a single orbital)
		 Bond operation: "Fe:O:+1+0+0# 1u:2d".	(operate between orbitals and spin)
		 */
		if( optList.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(optList[0], ":");
		auto secondSec = split(optList[1], ":");
		if( firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto pair = Lat.getPair(firstSec[2]);
		if( pair.atomI.atomName != firstSec[0])	return;
		if( pair.atomJ.atomName != firstSec[1])	return;
		if( !pair.withinRange() )				return;
		
		//if(!pair.atomI.hasOrbital(secondSec[0]))return;
		//if(!pair.atomJ.hasOrbital(secondSec[1]))return;
		
		x_mat sVec = parseBondString(pair, varList[0]);
		x_var val = sVec[0];
		for( unsigned i=1 ; i<varList.size() ; i++){
			val = val * parseBondString(pair, varList[i])[0];
		}
		
		auto spinIndexList_i = pair.atomI.spinIndexList(secondSec[0]);
		auto spinIndexList_j = pair.atomJ.spinIndexList(secondSec[1]);
		
		if(spinIndexList_i.size() <= 0)
			ErrorMessage("Error, in the operation: \n'"+opt+
						 "'\n atom:'"+pair.atomI.atomName+"' don't have any orbital to operate.");
		if(spinIndexList_j.size() <= 0)
			ErrorMessage("Error, in the operation: \n'"+opt+
						 "'\n atom:'"+pair.atomJ.atomName+"' don't have any orbital to operate.");
		
		if(spinIndexList_i.size() != spinIndexList_j.size()){
			ErrorMessage("Error, operation:\n"+opt+".\n The 'i' and 'j' spin-operation-list does not match eatch other.");
		}
		
		for( unsigned ii=0 ; ii<spinIndexList_i.size() ; ii++){
			auto & tmpChar = spinIndexList_i[ii].first[0];
			
			auto & index_i = spinIndexList_i[ii].second;
			auto & index_j = spinIndexList_j[ii].second;
			
			if(tmpChar == 'n' or tmpChar == 'u' or tmpChar == 'd'){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, pair.bondIJ()));
				hamElementList.push_back(MatrixElement(index_j, index_i, conj(val), -pair.bondIJ()));
			}
			if(tmpChar == 'A'){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, pair.bondIJ()));
				hamElementList.push_back(MatrixElement(index_j, index_i, conj(val), -pair.bondIJ()));
			}
			if(tmpChar == 'B'){
				hamElementList.push_back(MatrixElement(index_i, index_j,-val, pair.bondIJ()));
				hamElementList.push_back(MatrixElement(index_j, index_i, conj(-val),-pair.bondIJ()));
			}
		}
	}
	
	void		addScreenCoulomb(string opt, deque<string> optList, deque<string> varList)	{
		
		double radius = Lat.parameter.VAR("bondRadius", 0).real();
		
		if( optList.size()	!= 2)	{ ErrorMessage("Error, not a valid operation:\n "+opt); }
		if( optList[1][0]	!= '~')	{ ErrorMessage("Error, not a valid operation:\n "+opt); }
		if( varList.size()	< 2)	{ ErrorMessage("Error, not enough arguments:\n "+opt); }
		
		
		
		optList[1][0] = ' ';
		double tmpRadius = StrToDouble(optList[1]);
		if( tmpRadius > radius){
			 ErrorMessage("Error, in operation:\n "+opt
						  +"\n the given radius,~"+DoubleToStr(tmpRadius)+", "
						  +"is oarger the 'bondRadius(="+DoubleToStr(radius)+")'.");
		}
		radius = tmpRadius;
		
		auto atomI = Lat.getAtom();
		if( atomI.atomName != optList[0]){ return; } // Name does not match.
		
		auto orderKey = varList[0];
		r_var alpha = parseSiteString(atomI, varList[1])[0].real();
		for( unsigned i=2 ; i<varList.size() ; i++){
			alpha = alpha * parseSiteString(atomI, varList[i])[0].real();
		}

		//*****************************************************************
		//*        Query the neighbor atoms inside a sphere.              *
		auto rboxAtoms = Lat.getRBox(radius);
		//*                                                               *
		//*****************************************************************
		
		//*****************************************************************
		// Sum ove the coulomb term with the given charges (Z_i + Den_i).
		//*****************************************************************
		double sumScreenCoulomb = 0;
		for( auto & atomJ: rboxAtoms.second ){
			auto orderJ = order.findOrder(atomJ, orderKey);
			
			double Charge = -Lat.coreCharge.getCharge(atomJ.atomName);
			if( orderJ.first ){ Charge += orderJ.second[0].real(); }
			
			r_mat vecIJ(1,3);
			for(unsigned i=0 ; i<vecIJ.size() ; i++)
				vecIJ[i] = atomJ.pos[i] - atomI.pos[i];
			
			double distIJ = sqrt(cdot(vecIJ,vecIJ));
			
			if( distIJ > 0)
			sumScreenCoulomb += alpha * Charge * exp(-distIJ/radius) / distIJ;
		}
		
		//*****************************************************************
		// Attribute the surrounding Mean-field Coulomb potential to the Hamiltonian.
		//*****************************************************************
		x_mat coulomb(1,1);
		coulomb[0] = sumScreenCoulomb;
		order.setNew(atomI.atomIndex, "@:coulomb", coulomb);
		
		auto & orbitalIndex_I = atomI.allIndexList();
		for( auto & index_I: orbitalIndex_I){
			
			char ph_char_I = index_I.first[index_I.first.size()-2];
			if( ph_char_I != 'B' ){
				hamElementList.push_back(MatrixElement(index_I.second, index_I.second, sumScreenCoulomb,	vec(0,0,0))); continue;
			}
			else{
				hamElementList.push_back(MatrixElement(index_I.second, index_I.second,-sumScreenCoulomb,	vec(0,0,0))); continue;
			}
		}
	}
	
	void		addPairingS		(string opt, deque<string> optList, deque<string> varList)	{
		/*
		 Bond operation: "Fe	1:2".		(site only operation)
		 Bond operation: "Fe:O:+1+0+0 1:2".	(bond operation)
		 */
		if( Lat.HSpace() == NORMAL){
			 ErrorMessage("Error, pairing operation:\n"+opt+".\n Cannot be applited with normal space.");
		}
		if( Lat.parameter.STR("spin") == "off" and Lat.HSpace() == NAMBU){
			 ErrorMessage("Error, singlet pairing operation:\n"+opt+".\n Cannot be applited with spinless-nambu space.");
		}
		
		if( optList.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec	= split(optList[0], ":");
		auto secondSec	= split(optList[1], ":");
		if( firstSec.size() != 1 and firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size()!= 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto atomI = Lat.getAtom();
		AtomPair pair(atomI, atomI, vec(0,0,0));
		if( firstSec.size() == 3 ) { pair = Lat.getPair(firstSec[2]); }
		if( firstSec.size() == 1 ) { firstSec.push_back(firstSec[0]); }
		if( firstSec.size() == 3 and !pair.withinRange() )				return;
		if( pair.atomI.atomName != firstSec[0])	return;
		if( pair.atomJ.atomName != firstSec[1])	return;
		if(!pair.atomI.hasOrbital(secondSec[0]))return;
		if(!pair.atomJ.hasOrbital(secondSec[1]))return;
		
		x_var val = 0;
		if( firstSec.size() == 2 ){ // site only pairing
			val = parseSiteString(atomI, varList[0])[0];
			for( unsigned i=1 ; i<varList.size() ; i++){
				val = val * parseSiteString(atomI, varList[i])[0];
			}
	
		}
		if( firstSec.size() == 3 ){
			val = parseBondString(pair, varList[0])[0];
			for( unsigned i=1 ; i<varList.size() ; i++){
				val = val * parseBondString(pair, varList[i])[0];
			}
		}
		
		auto subIndexI = pair.atomI.orbitalIndexList(secondSec[0]);
		auto subIndexJ = pair.atomJ.orbitalIndexList(secondSec[1]);
		for( unsigned ii = 0; ii<subIndexI.size() ; ii++)
		for( unsigned jj = 0; jj<subIndexI.size() ; jj++){
			auto & iterI = subIndexI[ii];
			auto & iterJ = subIndexJ[jj];
			
			if(		(iterI.first == "Au" and iterJ.first == "Bd")
			   or	(iterI.first == "Ad" and iterJ.first == "Bu")){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));
				continue;
			}
			if(		(iterI.first == "Bu" and iterJ.first == "Ad")
			   or	(iterI.first == "Bd" and iterJ.first == "Au")){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val), -pair.bondIJ()));
				continue;
			}
		}
	}
	void		addPairingU		(string opt, deque<string> optList, deque<string> varList)	{
		/*
		 Bond operation: "Fe	1:2".		(site only operation)
		 Bond operation: "Fe:O:+1+0+0 1:2".	(bond operation)
		 */
		if( Lat.HSpace() == NORMAL){
			 ErrorMessage("Error, pairing operation:\n"+opt+".\n Cannot be applited with normal space.");
		}
		if( Lat.parameter.STR("spin") == "off" and Lat.HSpace() == NAMBU){
			 ErrorMessage("Error, singlet pairing operation:\n"+opt+".\n Cannot be applited with spinless-nambu space.");
		}
		
		if( optList.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec	= split(optList[0], ":");
		auto secondSec	= split(optList[1], ":");
		if( firstSec.size() != 1 and firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size()!= 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto atomI = Lat.getAtom();
		AtomPair pair(atomI, atomI, vec(0,0,0));
		if( firstSec.size() == 3 ) { pair = Lat.getPair(firstSec[2]); }
		if( firstSec.size() == 1 ) { firstSec.push_back(firstSec[0]); }
		if( firstSec.size() == 3 and !pair.withinRange() )				return;
		if( pair.atomI.atomName != firstSec[0])	return;
		if( pair.atomJ.atomName != firstSec[1])	return;
		if(!pair.atomI.hasOrbital(secondSec[0]))return;
		if(!pair.atomJ.hasOrbital(secondSec[1]))return;
		
		x_var val = 0;
		if( firstSec.size() == 2 ){ // site only pairing
			val = parseSiteString(atomI, varList[0])[0];
			for( unsigned i=1 ; i<varList.size() ; i++){
				val = val * parseSiteString(atomI, varList[i])[0];
			}
	
		}
		if( firstSec.size() == 3 ){
			val = parseBondString(pair, varList[0])[0];
			for( unsigned i=1 ; i<varList.size() ; i++){
				val = val * parseBondString(pair, varList[i])[0];
			}
		}
		
		auto subIndexI = pair.atomI.orbitalIndexList(secondSec[0]);
		auto subIndexJ = pair.atomJ.orbitalIndexList(secondSec[1]);
		for( unsigned ii = 0; ii<subIndexI.size() ; ii++)
		for( unsigned jj = 0; jj<subIndexI.size() ; jj++){
			auto & iterI = subIndexI[ii];
			auto & iterJ = subIndexJ[jj];
			
			if(	(iterI.first == "Au" and iterJ.first == "Bu")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));		continue;}
			if(	(iterI.first == "Bu" and iterJ.first == "Au")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val), -pair.bondIJ()));	continue;}
		}
	}
	void		addPairingD		(string opt, deque<string> optList, deque<string> varList)	{
		/*
		 Bond operation: "Fe	1:2".		(site only operation)
		 Bond operation: "Fe:O:+1+0+0 1:2".	(bond operation)
		 */
		if( Lat.HSpace() == NORMAL){
			 ErrorMessage("Error, pairing operation:\n"+opt+".\n Cannot be applited with normal space.");
		}
		if( Lat.parameter.STR("spin") == "off" and Lat.HSpace() == NAMBU){
			 ErrorMessage("Error, singlet pairing operation:\n"+opt+".\n Cannot be applited with spinless-nambu space.");
		}
		
		if( optList.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec	= split(optList[0], ":");
		auto secondSec	= split(optList[1], ":");
		if( firstSec.size() != 1 and firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size()!= 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto atomI = Lat.getAtom();
		AtomPair pair(atomI, atomI, vec(0,0,0));
		if( firstSec.size() == 3 ) { pair = Lat.getPair(firstSec[2]); }
		if( firstSec.size() == 1 ) { firstSec.push_back(firstSec[0]); }
		if( firstSec.size() == 3 and !pair.withinRange() )				return;
		if( pair.atomI.atomName != firstSec[0])	return;
		if( pair.atomJ.atomName != firstSec[1])	return;
		if(!pair.atomI.hasOrbital(secondSec[0]))return;
		if(!pair.atomJ.hasOrbital(secondSec[1]))return;
		
		x_var val = 0;
		if( firstSec.size() == 2 ){ // site only pairing
			val = parseSiteString(atomI, varList[0])[0];
			for( unsigned i=1 ; i<varList.size() ; i++){
				val = val * parseSiteString(atomI, varList[i])[0];
			}
	
		}
		if( firstSec.size() == 3 ){
			val = parseBondString(pair, varList[0])[0];
			for( unsigned i=1 ; i<varList.size() ; i++){
				val = val * parseBondString(pair, varList[i])[0];
			}
		}
		
		auto subIndexI = pair.atomI.orbitalIndexList(secondSec[0]);
		auto subIndexJ = pair.atomJ.orbitalIndexList(secondSec[1]);
		for( unsigned ii = 0; ii<subIndexI.size() ; ii++)
		for( unsigned jj = 0; jj<subIndexI.size() ; jj++){
			auto & iterI = subIndexI[ii];
			auto & iterJ = subIndexJ[jj];
			
			if(	(iterI.first == "Ad" and iterJ.first == "Bd")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));		continue;}
			if(	(iterI.first == "Bd" and iterJ.first == "Ad")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val), -pair.bondIJ()));	continue;}
		}
	}
	void		addPairingT		(string opt, deque<string> optList, deque<string> varList)	{
		/*
		 Bond operation: "Fe	1:2".		(site only operation)
		 Bond operation: "Fe:O:+1+0+0 1:2".	(bond operation)
		 */
		if( Lat.HSpace() == NORMAL){
			 ErrorMessage("Error, pairing operation:\n"+opt+".\n Cannot be applited with normal space.");
		}
		if( Lat.parameter.STR("spin") == "off" and Lat.HSpace() == NAMBU){
			 ErrorMessage("Error, singlet pairing operation:\n"+opt+".\n Cannot be applited with spinless-nambu space.");
		}
		
		if( optList.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec	= split(optList[0], ":");
		auto secondSec	= split(optList[1], ":");
		if( firstSec.size() != 1 and firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size()!= 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto atomI = Lat.getAtom();
		AtomPair pair(atomI, atomI, vec(0,0,0));
		if( firstSec.size() == 3 ) { pair = Lat.getPair(firstSec[2]); }
		if( firstSec.size() == 1 ) { firstSec.push_back(firstSec[0]); }
		if( firstSec.size() == 3 and !pair.withinRange() )				return;
		if( pair.atomI.atomName != firstSec[0])	return;
		if( pair.atomJ.atomName != firstSec[1])	return;
		if(!pair.atomI.hasOrbital(secondSec[0]))return;
		if(!pair.atomJ.hasOrbital(secondSec[1]))return;
		
		x_var val = 0;
		if( firstSec.size() == 2 ){ // site only pairing
			val = parseSiteString(atomI, varList[0])[0];
			for( unsigned i=1 ; i<varList.size() ; i++){
				val = val * parseSiteString(atomI, varList[i])[0];
			}
	
		}
		if( firstSec.size() == 3 ){
			val = parseBondString(pair, varList[0])[0];
			for( unsigned i=1 ; i<varList.size() ; i++){
				val = val * parseBondString(pair, varList[i])[0];
			}
		}
		auto subIndexI = pair.atomI.orbitalIndexList(secondSec[0]);
		auto subIndexJ = pair.atomJ.orbitalIndexList(secondSec[1]);
		for( unsigned ii = 0; ii<subIndexI.size() ; ii++)
		for( unsigned jj = 0; jj<subIndexI.size() ; jj++){
			auto & iterI = subIndexI[ii];
			auto & iterJ = subIndexJ[jj];
			
			if(		(iterI.first == "A"  and iterJ.first == "B" )){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,	val,	pair.bondIJ())); continue;}
			if(		(iterI.first == "B"  and iterJ.first == "A" )){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val),-pair.bondIJ())); continue;}
			if(		(iterI.first == "Au" and iterJ.first == "Bu")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,	val,	pair.bondIJ())); continue;}
			if(		(iterI.first == "Bu" and iterJ.first == "Au")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val),-pair.bondIJ())); continue;}
			if(		(iterI.first == "Ad" and iterJ.first == "Bd")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,	val,	pair.bondIJ())); continue;}
			if(		(iterI.first == "Bd" and iterJ.first == "Ad")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val),-pair.bondIJ())); continue;}
		}
	}

	void		addFieldB		(string opt, deque<string> optList, deque<string> varList)	{ // svar === " @:cspin * Jse "
		// *******************************
		// fieldB	> Fe   > [0,0,1] * B
		// fieldB	> Fe 1 > [0,0,1] * B
		// fieldB	> Fe 2 > [0,0,1] * B
		// *******************************
		
		if( optList.size() > 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( optList.size() == 1){
			// This part will not be effective in the quantum Hamiltonian, but goes to the Classical spin part.
			// That handeling the following operations.
			// fieldB	> Fe   > [0,0,1] * B (ignored)
			return;
		}
		
		auto dummyParser = split(opt, ":");
		if( dummyParser.size() != 1){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto atomI = Lat.getAtom();
		if( atomI.atomName != optList[0])	return;
		
		auto vectorB = parseSiteString(atomI, varList[0]);
		r_var multiplyB = 1.0;
		if( varList.size() == 2) multiplyB = parseSiteString(atomI, varList[1])[0].real();
		
		vectorB = vectorB * multiplyB;
		
		string orbital = optList[1];
		r_var Sx =-vectorB[0].real();
		r_var Sy =-vectorB[1].real();
		r_var Sz =-vectorB[2].real();
		switch( Lat.HSpace() ){
			case NORMAL: {
				unsigned indexNu = atomI.index(orbital+"u");
				unsigned indexNd = atomI.index(orbital+"d");
				
				hamElementList.push_back(MatrixElement(indexNu, indexNd, Sx-Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexNd, indexNu, Sx+Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexNu, indexNu, Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexNd, indexNd,-Sz			, vec(0,0,0)));
				break;
			}
			case NAMBU: {
				unsigned indexAu = atomI.index(orbital+"Au");
				unsigned indexBd = atomI.index(orbital+"Bd");
				hamElementList.push_back(MatrixElement(indexAu, indexAu, Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBd, Sz			, vec(0,0,0)));
				if( Sx > 0.000001 or Sy > 0.000001){
					ErrorMessage("Error, spin-flip operation:\n"+
							 opt+
							 ".\n Cannot be applied for a spin-depende Nambu space Hamiltonian.");
				}
				break;
			}
			case EXNAMBU: {
				unsigned indexAu = atomI.index(orbital+"Au");
				unsigned indexAd = atomI.index(orbital+"Ad");
				unsigned indexBu = atomI.index(orbital+"Bu");
				unsigned indexBd = atomI.index(orbital+"Bd");
				
				hamElementList.push_back(MatrixElement(indexAu, indexAd, Sx-Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexAd, indexAu, Sx+Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexAu, indexAu, Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexAd, indexAd,-Sz			, vec(0,0,0)));
				
				hamElementList.push_back(MatrixElement(indexBu, indexBd,-Sx+Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBu,-Sx-Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBu, indexBu,-Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBd, Sz			, vec(0,0,0)));
				break;
			}
		}
	}
	
	void		addChemicalPotential(r_var mu)				{
		
		// Apply the Chemical potential
		//r_var mu = Lat.parameter.VAR("Mu",0).real();
		while( Lat.iterate() ){
			auto atomI = Lat.getAtom();
			for( auto & index: atomI.allIndexList()){
				auto labelParser = split(index.first, ".");
				
				if( labelParser.size() == 2){
					string sub_space = labelParser[1];
					if(sub_space[0] == 'n' or sub_space[0] == 'u' or sub_space[0] == 'd'){
						hamElementList.push_back(MatrixElement(index.second, index.second,-mu, vec(0,0,0)));
					}
					if(sub_space[0] == 'A'){
						hamElementList.push_back(MatrixElement(index.second, index.second,-mu, vec(0,0,0)));
					}
					if(sub_space[0] == 'B'){
						hamElementList.push_back(MatrixElement(index.second, index.second, mu, vec(0,0,0)));
					}
				}
			}
		}
	}
	
	/*------------------------------------------------
	 Using these methods to construct and diagonalize (in k-space) Hamiltonian.
	 -------------------------------------------------*/
	void		constructTBMHam	(bool withMu = true)		{
		
		initHam();
		
		//if( hamElementList.size() > 0 ) return;
		
		maxE = -10000;
		minE =  10000;
		
		//unsigned iter = 0;
		//cout<<"Entering debugging view"<<endl;
		//
		//auto  pair = Lat.getPair(0, vec(1,0,0)); //***
		//while(true){
		//	auto ddd = split("AAA,BBB,CCC,EEE", ",");
		//	hamElementList.push_back(MatrixElement(0,0,0,vec(0,0,0)));
		//	iter++;
		//	cout<<iter<<" "<<hamElementList.size()<<endl;
		//	if( iter % 2000 == 0){
		//		hamElementList.clear();
		//	}
		//	while( Lat.iterate() ){
		//		auto atom = Lat.getAtom();
		//		x_mat	sVec = parseSiteString(atom, "0");
		//		Lat.vec("+1+0+0#"); //---
		//		Lat.vec("+1+0+0"); //---
		//		auto  pair = Lat.getPair(vec(1,0,0)); //***
		//		x_mat tmpvec = parseBondString(pair, "tdp"); //---
		//	}
		//}
		//cout<<"Leaving debugging view"<<endl;
		
		/* Construct the hamElementList from Lat.hamParser.hamOperationList from the 'xxx.lat.tbm'. */
		hamElementList.clear();
		while( Lat.iterate() ) {
			
			for( auto & iter : Lat.hamParser.getOperationListMap("hundSpin"))	{	addHundSpin		(iter.get<0>(), iter.get<1>(), iter.get<2>());}
			for( auto & iter : Lat.hamParser.getOperationListMap("orbital")	)	{	addOrbitalEng	(iter.get<0>(), iter.get<1>(), iter.get<2>());}
			for( auto & iter : Lat.hamParser.getOperationListMap("site")	)	{	addSiteCouple	(iter.get<0>(), iter.get<1>(), iter.get<2>());}
			for( auto & iter : Lat.hamParser.getOperationListMap("siteHc")	)	{	addSiteCoupleHc	(iter.get<0>(), iter.get<1>(), iter.get<2>());}
			for( auto & iter : Lat.hamParser.getOperationListMap("hopping")	)	{	addHoppingInt	(iter.get<0>(), iter.get<1>(), iter.get<2>());}
			for( auto & iter : Lat.hamParser.getOperationListMap("hoppingHc"))	{	addHoppingIntHc	(iter.get<0>(), iter.get<1>(), iter.get<2>());}
			for( auto & iter : Lat.hamParser.getOperationListMap("bond")	)	{	addBondCouple	(iter.get<0>(), iter.get<1>(), iter.get<2>());}
			for( auto & iter : Lat.hamParser.getOperationListMap("bondHc")	)	{	addBondCoupleHc	(iter.get<0>(), iter.get<1>(), iter.get<2>());}
			for( auto & iter : Lat.hamParser.getOperationListMap("screenCoulomb")){	addScreenCoulomb(iter.get<0>(), iter.get<1>(), iter.get<2>());}
			for( auto & iter : Lat.hamParser.getOperationListMap("fieldB")	)	{	addFieldB		(iter.get<0>(), iter.get<1>(), iter.get<2>());}
			
			// If space==normal The following part will be ignored.
			if( Lat.HSpace() == NORMAL) continue;
			
			for( auto & iter : Lat.hamParser.getOperationListMap("pairingS"))	{	addPairingS	(iter.get<0>(), iter.get<1>(), iter.get<2>()); }
			for( auto & iter : Lat.hamParser.getOperationListMap("pairingU"))	{	addPairingU	(iter.get<0>(), iter.get<1>(), iter.get<2>()); }
			for( auto & iter : Lat.hamParser.getOperationListMap("pairingD"))	{	addPairingD	(iter.get<0>(), iter.get<1>(), iter.get<2>()); }
			for( auto & iter : Lat.hamParser.getOperationListMap("pairingT"))	{	addPairingT	(iter.get<0>(), iter.get<1>(), iter.get<2>()); }
		}
		if( withMu ) addChemicalPotential(Lat.parameter.VAR("Mu").real());
	}

	EigVec  	HamEvd			(r_mat kp)					{
		
		Ham.zerolize();
		
		for( auto & elem: hamElementList){
			
			auto bvec = elem.bondVector;
			x_var phase = kp[0]*bvec[0]+ kp[1]*bvec[1]+ kp[2]*bvec[2];
			auto	exp_IJ=exp(-Im*phase);
			
			Ham(elem.I, elem.J) += elem.val * exp_IJ;
		}
		
		EigVec					tmpEigValVec;	// A temporary storage of eigen value vector for a specific K-point.
		tmpEigValVec.kPoint = kp;
		
		tmpEigValVec.message = Ham.evd(tmpEigValVec.eigenValue, tmpEigValVec.eigenVector);
		
		if ( maxE < tmpEigValVec.eigenValue[Lat.indexSize()-1]	){ maxE = tmpEigValVec.eigenValue[Lat.indexSize()-1]; }
		if ( minE > tmpEigValVec.eigenValue[0]					){ minE = tmpEigValVec.eigenValue[0]; }
		
		return tmpEigValVec;
	}
	
	x_var		getDensityMatrix		(unsigned index_i, unsigned index_j, double Mu=0){
		
		if( Lat.HSpace() != NORMAL) Mu = 0;
		
		auto Temperature = Lat.parameter.VAR("Temperature", 0.00001);
		if ( KEigenValVec.size()>0  and index_i < Lat.indexSize() and index_j < Lat.indexSize() ) {
			x_var den_ij=0;
			for (unsigned k=0; k<KEigenValVec.size(); k++) {
				auto & kEigVal = KEigenValVec[k].eigenValue;
				auto & kEigVec = KEigenValVec[k].eigenVector;
    
				for( int i=0 ; i<kEigVal.size() ; i++ ){
					double EigF = 1.0/(1+exp( (kEigVal[i]-Mu)/Temperature.real()));
					x_var vI = kEigVec(index_i,i);
					x_var vJ = kEigVec(index_j,i);
					den_ij += conj(vI)*vJ*EigF;
				}
			}
			return den_ij/KEigenValVec.size();
		}
		
		return 0;
	}
	
	void		calculate4DensityOrder	()							{
		if( Lat.parameter.STR("spin")!="on" ){
			ErrorMessage("Error, the spin space suld be turned \'on\' for the 4-density calculation");
		}
		if( Lat.HSpace()==NAMBU ){
			ErrorMessage("Error, NAMBU space is not applicable forthe 4-density calculation");
		}
		
		//OrderParameter	newOrder( order );
			
		while( Lat.iterate() ){
			auto atomI = Lat.getAtom();
			
			map<string, x_mat> fourDensity;
			x_mat	totalDen(1,1);
			
			auto & indexLabel= atomI.allIndexList() ;
			for( unsigned i=0 ; i<indexLabel.size() ; i++)
			for( unsigned j=0 ; j<indexLabel.size() ; j++) {
				auto & indexLabel_I = indexLabel[i];
				auto & indexLabel_J = indexLabel[j];
				
				auto index_I = indexLabel_I.second;
				auto index_J = indexLabel_J.second;
				auto parser_I = split(indexLabel_I.first, ".");
				auto parser_J = split(indexLabel_J.first, ".");
				
				bool theSame_AB_space = true;
				if( Lat.HSpace() == EXNAMBU and parser_I[1][0] != parser_J[1][0]) theSame_AB_space = false;
				
				
				if( parser_I[0] == parser_J[0] and theSame_AB_space){
					
					if( fourDensity.find(parser_I[0]) == fourDensity.end() ) {
						fourDensity[parser_I[0]] = x_mat(1,4);
					}
					auto & r4den = fourDensity[parser_I[0]];
					
					if( parser_I[1] == "u" and parser_J[1] == "u"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += tmpDen;
						r4den[3] += tmpDen;
						totalDen[0] += tmpDen;
					}
					if( parser_I[1] == "d" and parser_J[1] == "d"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += tmpDen;
						r4den[3] -= tmpDen;
						totalDen[0] += tmpDen;
					}
					if( parser_I[1] == "u" and parser_J[1] == "d"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[1] += tmpDen.real();
						r4den[2] -= tmpDen.imag();
					}
					if( parser_I[1] == "d" and parser_J[1] == "u"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[1] += tmpDen.real();
						r4den[2] += tmpDen.imag();
					}
					
					if( parser_I[1] == "Au" and parser_J[1] == "Au"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += tmpDen*0.5;
						r4den[3] += tmpDen*0.5;
						totalDen[0] += tmpDen*0.5;
					}
					if( parser_I[1] == "Ad" and parser_J[1] == "Ad"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += tmpDen*0.5;
						r4den[3] -= tmpDen*0.5;
						totalDen[0] += tmpDen*0.5;
					}                     
					if( parser_I[1] == "Au" and parser_J[1] == "Ad"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[1] += tmpDen.real()*0.5;
						r4den[2] -= tmpDen.imag()*0.5;
					}                            
					if( parser_I[1] == "Ad" and parser_J[1] == "Au"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[1] += tmpDen.real()*0.5;
						r4den[2] += tmpDen.imag()*0.5;
					}                            
					
					if( parser_I[1] == "Bu" and parser_J[1] == "Bu"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += 1-tmpDen*0.5;
						r4den[3] += 1-tmpDen*0.5;
						totalDen[0] += 1-tmpDen*0.5;
					}                       
					if( parser_I[1] == "Bd" and parser_J[1] == "Bd"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += 1-tmpDen*0.5;
						r4den[3] -= 1-tmpDen*0.5;
						totalDen[0] += 1-tmpDen*0.5;
					}                       
					if( parser_I[1] == "Bu" and parser_J[1] == "Bd"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[1] += -tmpDen.real()*0.5;
						r4den[2] -= -tmpDen.imag()*0.5;
					}                             
					if( parser_I[1] == "Bd" and parser_J[1] == "Bu"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[1] += -tmpDen.real()*0.5;
						r4den[2] += -tmpDen.imag()*0.5;
					}
					for( unsigned ii=0 ; ii<r4den.size() ; ii++)
						if( abs(r4den[ii]) < 0.000001 ) r4den[ii] = 0;
				}
			}
			
			//order(atomI.atomName+" den") = totalDen;
			
			for( auto & iter: fourDensity){
				order(atomI.atomName+" "+iter.first+":4den") = iter.second;
			}
		}
		order.save();
	}
	void		calculateEnergy			()							{

		double totalE=0;
		auto Temperature = Lat.parameter.VAR("Temperature", 0.00001).real();
		for ( auto & iter: KEigenValVec ){
			auto kEig = iter.eigenValue;
			for (unsigned i=0; i<kEig.size(); i++) {
				double EigF = 1.0/(1+exp( kEig[i]/Temperature));
				totalE+=EigF*kEig[i];
			}
		}
		
		totalE = totalE / KEigenValVec.size();
		
		if( KEigenValVec.size() > 0 ){
			energyMap["1.Q Eng"] = totalE;
		}
		else {
			energyMap["1.Q Eng"] = 0;
		}
	}
	
	x_mat		parseSiteString(Atom & at, string svar)			{
		// @:cspin .. Jh .. 1.0 .. [ 1, 2, 3] ...
		removeSpace(svar);
		
		x_mat xvar;
		
		if( svar[0] == '@'){ // Order parameter
			auto found = order.findOrder(at.atomName, svar);
			if( found.first ) { xvar = found.second; }
		}
		else if( svar[0] == '[' and svar[svar.size()-1] == ']'){ // Vector type
			replaceAll(svar, "[", "");
			replaceAll(svar, "]", "");
			xvar = StrToXVec(svar);
		}
		else if( IsFloatStr(svar)){ // a number
			x_mat tmp(1,1);
			tmp[0] = StrToComplex(svar);
			xvar = tmp;
		}
		else { // Take defined input variable.
			x_mat tmp(1,1);
			tmp[0] = Lat.parameter.VAR(svar, 1);
			xvar = tmp;
		}
		
		return xvar;
	}
	
	x_mat		parseBondString(AtomPair & ap, string svar)		{
		// @:cspin .. Jh .. 1.0 .. [ 1, 2, 3] ...
		removeSpace(svar);
		
		x_mat xvar;
		
		
		replaceAll(svar, "@", vecToStr(ap.bondIJ()) );
		if( svar[0] == '@'){ // Order parameter
			auto found = order.findOrder(ap.atomI.atomName, svar);
			if( found.first ) { xvar = found.second; }
		}
		else if( svar[0] == '[' and svar[svar.size()-1] == ']'){ // Vector type
			replaceAll(svar, "[", "");
			replaceAll(svar, "]", "");
			xvar = StrToXVec(svar);
		}
		else if( svar[0] == '(' and svar[svar.size()-1] == ')'){ // Complex type
			x_mat tmp(1,1);
			tmp[0] = StrToComplex(svar);
			xvar = tmp;
		}
		else if( IsFloatStr(svar)){ // a number
			x_mat tmp(1,1);
			tmp[0] = StrToComplex(svar);
			xvar = tmp;
		}
		else { // Take defined input variable.
			x_mat tmp(1,1);
			tmp[0] = Lat.parameter.VAR(svar, 1);
			xvar = tmp;
		}
		
		return xvar;
	}
	
	vector<x_mat>	siteVec;
	vector<x_mat> &	parseSiteVecString(Atom & at, string svec)		{ // @:cspin * Jh * ...
		
		siteVec.clear();
		removeSpace(svec);
		auto varList = split(svec, "*");
		
		//vector<x_mat> xvec;
		for( auto & elem: varList) siteVec.push_back( parseSiteString(at,elem));
		
		return siteVec;
	}
	
	vector<x_mat>	bondVec;
	vector<x_mat> &	parseBondVecString(AtomPair & ap, string svec)	{ // +1+0+0:1:1:pairS * Jh * ...
		
		bondVec.clear();
		removeSpace(svec);
		auto varList = split(svec, "*");
		
		//vector<x_mat> xvec;
		for( auto & elem: varList) bondVec.push_back(parseBondString(ap, elem));
		
		return bondVec;
	}

};



