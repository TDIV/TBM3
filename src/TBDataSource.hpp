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
//  TBDataSource.hpp
//  TBM^3
//

class	MatrixElement{
public:
	int		I,J;
	r_mat	bondVector;
	x_var	val;
	
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
	string	status;
	r_mat	kPoint;
	r_mat	eigenValue;
	x_mat	eigenVector;
};

/*
 This structure will be passed through all future models.
 A. It can be used to construct all matrix element, for example:
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
	
	Lattice &				Lat;
	OrderParameter			order;
	OrderParameter			order_old;
	
	vector<MatrixElement>	hamElementList;	// A list to store all the matrix elements
	x_mat					Ham;			// The body of the Hamiltonian (matrix).
	vector<EigVec>			KEigenValVec;	// Calculated K-space dependent eigen value and vectore.
	EigVec					tmpEigValVec;	// A temporary storage of eigen value vector for a specific K-point.
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
	void addHundSpin	(string opt, string svar)	{
		/* Site operation: "Fe 1".	(operate for only a single orbital) */
		
		replaceAll(opt, "\t", " ");
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage("Error, 'HundSpin' operation:\n"+opt+".\n Cannot be applied for a spin-independent Hamiltonian.");
		}
		
		auto parser = split(opt, " ");
		if( parser.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto atomI = Lat.getAtom();
		if( atomI.atomName != parser[0]){ return; } // Name does not match.
		if(!atomI.hasOrbital(parser[1])){ return; } // No such orbital.
		
		auto svarParser = split(svar, "*");
		auto orderKey = svarParser[0];
		r_var Jh = 1.0;
		x_mat xvec = parseSiteString(atomI, svarParser[0]);
		if( svarParser.size() >= 2){
			for( unsigned i=1 ; i<svarParser.size() ; i++){
				Jh = Jh * parseSiteString(atomI, svarParser[i])[0].real();
			}
		}
		
		normalizeXVec(xvec);
		x_mat sVec = xvec * Jh;
	
		if( sVec.size() != 3 ){
			ErrorMessage("Error, operation:\n"+opt+".\n Has been applied with a wrong spin-variable with size:"+IntToStr(sVec.size()));
		}
		
		//normalizeXVec( sVec );
		
		string orbital = parser[1];
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
	void addOrbitalEng	(string opt, string svar)	{
		/*
		 Site operation: "Fe 1".	(operate for only a single orbital)
		 */
		replaceAll(opt, "\t", " ");
		auto parser = split(opt, " ");
		if( parser.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto atomI = Lat.getAtom();
		if( atomI.atomName != parser[0]){ return; } // Name does not match.
		if(!atomI.hasOrbital(parser[1])){ return; } // No such orbital.
		
		auto xvec = parseSiteVecString(atomI, svar);
		x_mat sVec(1,3);
		if( xvec.size() >= 1 ){ sVec = xvec[0]; }
		for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
		r_var val = sVec[0].real();
		
		auto subIndex = atomI.orbitalIndexList(parser[1]);
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
	void addSiteCouple	(string opt, string svar)	{
		/*
		 Site operation: "Fe 1u:1u".	(operate for only a single orbital)
		 Site operation: "Fe 1u:2d".	(operate between orbitals and spin)
		 */
		replaceAll(opt, "\t", " ");
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage("Error, operation:\n"+opt+".\n Cannot be applied for a spin-independent Hamiltonian.");
		}
		
		auto parser = split(opt, " ");
		if( parser.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(parser[0], ":");
		auto secondSec = split(parser[1], ":");
		if( firstSec.size() != 1 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		if( secondSec[0][secondSec[0].size()-1] != secondSec[1][secondSec[1].size()-1] and
		   Lat.HSpace() == NAMBU) {
			ErrorMessage("Error, operation:\n"+opt+".\n Spin-flip term cannot be applied under 'nambu' space.\n Please set the 'space' to 'normal' or 'exnambu'.");
		}
		
		auto atomI = Lat.getAtom();
		if( parser[0] != atomI.atomName )	return;
		
		auto xvec = parseSiteVecString(atomI, svar);
		x_mat sVec(1,3);
		if( xvec.size() >= 1 ){ sVec = xvec[0]; }
		for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
		r_var val = sVec[0].real();
		
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
			//auto & label_j = spinIndexList_j[ii].first;
			auto & index_i = spinIndexList_i[ii].second;
			auto & index_j = spinIndexList_j[ii].second;
			
			if( label_i == "N" or label_i == "A")
				hamElementList.push_back(MatrixElement(index_i, index_j, val, vec(0,0,0)));
			if( label_i == "B")
				hamElementList.push_back(MatrixElement(index_i, index_j,-val, vec(0,0,0)));
		}
	}
	void addSiteCoupleHc(string opt, string svar)	{
		/*
		 Site operation: "Fe 1u:1u".	(operate for only a single orbital)
		 Site operation: "Fe 1u:2d".	(operate between orbitals and spin)
		 */
		replaceAll(opt, "\t", " ");
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage("Error, operation:\n"+opt+".\n Cannot apply operation with a spin-independent Hamiltonian.");
		}
		
		auto parser = split(opt, " ");
		if( parser.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(parser[0], ":");
		auto secondSec = split(parser[1], ":");
		if( firstSec.size() != 1 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		if( secondSec[0][secondSec[0].size()-1] != secondSec[1][secondSec[1].size()-1] and
		   Lat.HSpace() == NAMBU) {
			ErrorMessage("Error, operation:\n"+opt+".\n Spin-flip term cannot be applied under 'nambu' space.\n Please set the 'space' to 'normal' or 'exnambu'.");
		}
		
		auto atomI = Lat.getAtom();
		if( parser[0] != atomI.atomName )	return;
		
		auto xvec = parseSiteVecString(atomI, svar);
		x_mat sVec(1,3);
		if( xvec.size() >= 1 ){ sVec = xvec[0]; }
		for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
		x_var val = sVec[0];
		
		auto spinIndexList_i = atomI.spinIndexList(secondSec[0]);
		auto spinIndexList_j = atomI.spinIndexList(secondSec[1]);
		
		if(spinIndexList_i.size() <= 0)
			ErrorMessage("Error, operation: \n'"+opt+
						 "'\n atom:'"+atomI.atomName+"' don't have any orbital to operate.");
		if(spinIndexList_j.size() <= 0)
			ErrorMessage("Error, operation: \n'"+opt+
						 "'\n atom:'"+atomI.atomName+"' don't have any orbital to operate.");
		
		if(spinIndexList_i.size() != spinIndexList_j.size()){
			ErrorMessage("Error, operation:\n"+opt+".\n The 'i' and 'j' spin-operation-list does not match eatch other.");
		}
		
		for( unsigned ii=0 ; ii<spinIndexList_i.size() ; ii++){
			auto & label_i = spinIndexList_i[ii].first;
			//auto & label_j = spinIndexList_j[ii].first;
			auto & index_i = spinIndexList_i[ii].second;
			auto & index_j = spinIndexList_j[ii].second;
			
			if( label_i == "N" or label_i == "A"){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(index_j, index_i,conj(val), vec(0,0,0)));
			}
			if( label_i == "B"){
				hamElementList.push_back(MatrixElement(index_i, index_j,-val, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(index_j, index_i,-conj(val), vec(0,0,0)));
			}
		}
	}
	
	void addHoppingInt	(string opt, string svar)	{
		/*
		 Bond operation: "Fe:O:+1+0+0 1:2".	(operate for only a single orbital)
		 */
		replaceAll(opt, "\t", " ");
		auto parser = split(opt, " ");
		if( parser.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(parser[0], ":");
		auto secondSec = split(parser[1], ":");
		if( firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto pair = Lat.getPair(firstSec[2]);
		if( pair.atomI.atomName != firstSec[0])	return;
		if( pair.atomJ.atomName != firstSec[1])	return;
		if( !pair.withinRange() )				return;
		
		if(!pair.atomI.hasOrbital(secondSec[0]))return;
		if(!pair.atomJ.hasOrbital(secondSec[1]))return;
		
		auto xvec = parseBondVecString(pair, svar);
		x_mat sVec(1,3);
		if( xvec.size() >= 1 ){ sVec = xvec[0]; }
		for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
		x_var val = sVec[0];
		
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
	void addHoppingIntHc(string opt, string svar)	{
		/*
		 Pair operation: "Fe:O:+1+0+0 1:2".	(operate for only a single orbital)
		 */
		replaceAll(opt, "\t", " ");
		auto parser = split(opt, " ");
		if( parser.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(parser[0], ":");
		auto secondSec = split(parser[1], ":");
		if( firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto pair = Lat.getPair(firstSec[2]);
		if( pair.atomI.atomName != firstSec[0])	return;
		if( pair.atomJ.atomName != firstSec[1])	return;
		if( !pair.withinRange() )				return;
		
		if(!pair.atomI.hasOrbital(secondSec[0]))return;
		if(!pair.atomJ.hasOrbital(secondSec[1]))return;
		
		auto xvec = parseBondVecString(pair, svar);
		x_mat sVec(1,3);
		if( xvec.size() >= 1 ){ sVec = xvec[0]; }
		for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
		x_var val = sVec[0];
		
		auto subIndexI = pair.atomI.orbitalIndexList(secondSec[0]);
		auto subIndexJ = pair.atomJ.orbitalIndexList(secondSec[1]);
		for( unsigned ii = 0; ii<subIndexI.size() ; ii++){
			auto & iterI = subIndexI[ii];
			auto & iterJ = subIndexJ[ii];
			char & tmpChar = iterI.first[0];
			auto bondIJ = pair.bondIJ();
			
			if(tmpChar == 'n' or tmpChar == 'u' or tmpChar == 'd'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, bondIJ));
				hamElementList.push_back(MatrixElement(iterJ.second, iterI.second, conj(val),0-bondIJ));
			}
			if(tmpChar == 'A'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, bondIJ));
				hamElementList.push_back(MatrixElement(iterJ.second, iterI.second, conj(val),0-bondIJ));
			}
			if(tmpChar == 'B'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,-val, bondIJ));
				hamElementList.push_back(MatrixElement(iterJ.second, iterI.second,conj(-val),0-bondIJ));
			}
		}
	}

	void addBondCouple	(string opt, string svar)	{
		/*
		 Bond operation: "Fe:Fe:+1+0+0 1u:1u".	(operate for only a single orbital)
		 Bond operation: "Fe:O:+1+0+0# 1u:2d".	(operate between orbitals and spin)
		 */
		replaceAll(opt, "\t", " ");
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage("Error, operation:\n"+opt+".\n Cannot apply operation with a spin-independent Hamiltonian.");
		}
		
		auto parser = split(opt, " ");
		if( parser.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(parser[0], ":");
		auto secondSec = split(parser[1], ":");
		if( firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		if( secondSec[0][secondSec[0].size()-1] != secondSec[1][secondSec[1].size()-1] and
		   Lat.HSpace() == NAMBU) {
			ErrorMessage("Error, operation:\n"+opt+".\n Spin-flip term cannot be applied under 'nambu' space.\n Please set the 'space' to 'normal' or 'exnambu'.");
		}
		
		auto pair = Lat.getPair(firstSec[2]);
		if( pair.atomI.atomName != firstSec[0])	return;
		if( pair.atomJ.atomName != firstSec[1])	return;
		if( !pair.withinRange() )				return;
		
		auto xvec = parseBondVecString(pair, svar);
		x_mat sVec(1,3);
		if( xvec.size() >= 1 ){ sVec = xvec[0]; }
		for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
		r_var val = sVec[0].real();
		
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
			auto & label_i = spinIndexList_i[ii].first;
			//auto & label_j = spinIndexList_j[ii].first;
			auto & index_i = spinIndexList_i[ii].second;
			auto & index_j = spinIndexList_j[ii].second;
			
			if( label_i == "N" or label_i == "A")
				hamElementList.push_back(MatrixElement(index_i, index_j, val, pair.bondIJ()));
			if( label_i == "B")
				hamElementList.push_back(MatrixElement(index_i, index_j,-val, pair.bondIJ()));
		}
	}
	void addBondCoupleHc(string opt, string svar)	{
		/*
		 Bond operation: "Fe:Fe:+1+0+0 1u:1u".	(operate for only a single orbital)
		 Bond operation: "Fe:O:+1+0+0# 1u:2d".	(operate between orbitals and spin)
		 */
		replaceAll(opt, "\t", " ");
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage("Error, operation:\n"+opt+".\n Cannot apply operation with a spin-independent Hamiltonian.");
		}
		
		auto parser = split(opt, " ");
		if( parser.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(parser[0], ":");
		auto secondSec = split(parser[1], ":");
		if( firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		if( secondSec[0][secondSec[0].size()-1] != secondSec[1][secondSec[1].size()-1] and
		   Lat.HSpace() == NAMBU) {
			ErrorMessage("Error, operation:\n"+opt+".\n Spin-flip term cannot be applied under 'nambu' space.\n Please set the 'space' to 'normal' or 'exnambu'.");
		}
		
		auto pair = Lat.getPair(firstSec[2]);
		if( pair.atomI.atomName != firstSec[0])	return;
		if( pair.atomJ.atomName != firstSec[1])	return;
		if( !pair.withinRange() )				return;
		
		auto xvec = parseBondVecString(pair, svar);
		x_mat sVec(1,3);
		if( xvec.size() >= 1 ){ sVec = xvec[0]; }
		for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
		x_var val = sVec[0];
		
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
			auto & label_i = spinIndexList_i[ii].first;
			//auto & label_j = spinIndexList_j[ii].first;
			auto & index_i = spinIndexList_i[ii].second;
			auto & index_j = spinIndexList_j[ii].second;
			auto bondIJ = pair.bondIJ();
			if( label_i == "N" or label_i == "A"){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, bondIJ));
				hamElementList.push_back(MatrixElement(index_j, index_i,conj(val), 0-bondIJ));
			}
			if( label_i == "B"){
				hamElementList.push_back(MatrixElement(index_i, index_j,-val, bondIJ));
				hamElementList.push_back(MatrixElement(index_j, index_i,conj(-val), bondIJ));
			}
		}
	}
	
	void addPairingS	(string opt, string svar)	{
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
		
		
		replaceAll(opt, "\t", " ");
		auto parser = split(opt, " ");
		if( parser.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(parser[0], ":");
		auto secondSec = split(parser[1], ":");
		if( firstSec.size() != 1 and firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		
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
		if( firstSec.size() == 2 ){
			auto xvec = parseSiteVecString(atomI, svar);
			x_mat sVec(1,3);
			if( xvec.size() >= 1 ){ sVec = xvec[0]; }
			for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
			val = sVec[0].real();
		}
		if( firstSec.size() == 3 ){
			auto xvec = parseBondVecString(pair, svar);
			x_mat sVec(1,3);
			if( xvec.size() >= 1 ){ sVec = xvec[0]; }
			for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
			val = sVec[0].real();
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
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val), 0-pair.bondIJ()));
				continue;
			}
		}
	}
	void addPairingU	(string opt, string svar)	{
		/*
		 Bond operation: "Fe:O:+1+0+0 1:2".	(operate for only a single orbital)
		 */
		if( Lat.HSpace() == NORMAL or Lat.parameter.STR("spin") == "off"){
			 ErrorMessage("Error, pairing operation:\n"+opt+".\n Cannot be applited with normal space.");
		}
		if( Lat.HSpace() == NAMBU){
			 ErrorMessage("Error, Up-Up pairing operation:\n"+opt+".\n Cannot be applited with spinful-nambu space.");
		}
		
		replaceAll(opt, "\t", " ");
		auto parser = split(opt, " ");
		if( parser.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(parser[0], ":");
		auto secondSec = split(parser[1], ":");
		if( firstSec.size() != 1 and firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
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
		if( firstSec.size() == 2 ){
			auto xvec = parseSiteVecString(atomI, svar);
			x_mat sVec(1,3);
			if( xvec.size() >= 1 ){ sVec = xvec[0]; }
			for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
			val = sVec[0].real();
		}
		if( firstSec.size() == 3 ){
			auto xvec = parseBondVecString(pair, svar);
			x_mat sVec(1,3);
			if( xvec.size() >= 1 ){ sVec = xvec[0]; }
			for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
			val = sVec[0].real();
		}
		
		auto subIndexI = pair.atomI.orbitalIndexList(secondSec[0]);
		auto subIndexJ = pair.atomJ.orbitalIndexList(secondSec[1]);
		for( unsigned ii = 0; ii<subIndexI.size() ; ii++)
		for( unsigned jj = 0; jj<subIndexI.size() ; jj++){
			auto & iterI = subIndexI[ii];
			auto & iterJ = subIndexJ[jj];
			
			if(	(iterI.first == "Au" and iterJ.first == "Bu")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));		continue;}
			if(	(iterI.first == "Bu" and iterJ.first == "Au")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val), 0-pair.bondIJ()));	continue;}
		}
	}
	void addPairingD	(string opt, string svar)	{
		/*
		 Bond operation: "Fe:O:+1+0+0 1:2".	(operate for only a single orbital)
		 */
		if( Lat.HSpace() == NORMAL or Lat.parameter.STR("spin") == "off"){
			 ErrorMessage("Error, pairing operation:\n"+opt+".\n Cannot be applited with normal space.");
		}
		if( Lat.HSpace() == NAMBU){
			 ErrorMessage("Error, Dn-Dn pairing operation:\n"+opt+".\n Cannot be applited with spinful-nambu space.");
		}
		
		replaceAll(opt, "\t", " ");
		auto parser = split(opt, " ");
		if( parser.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(parser[0], ":");
		auto secondSec = split(parser[1], ":");
		if( firstSec.size() != 1 and firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
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
		if( firstSec.size() == 2 ){
			auto xvec = parseSiteVecString(atomI, svar);
			x_mat sVec(1,3);
			if( xvec.size() >= 1 ){ sVec = xvec[0]; }
			for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
			val = sVec[0].real();
		}
		if( firstSec.size() == 3 ){
			auto xvec = parseBondVecString(pair, svar);
			x_mat sVec(1,3);
			if( xvec.size() >= 1 ){ sVec = xvec[0]; }
			for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
			val = sVec[0].real();
		}
		
		auto subIndexI = pair.atomI.orbitalIndexList(secondSec[0]);
		auto subIndexJ = pair.atomJ.orbitalIndexList(secondSec[1]);
		for( unsigned ii = 0; ii<subIndexI.size() ; ii++)
		for( unsigned jj = 0; jj<subIndexI.size() ; jj++){
			auto & iterI = subIndexI[ii];
			auto & iterJ = subIndexJ[jj];
			
			if(	(iterI.first == "Ad" and iterJ.first == "Bd")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));		continue;}
			if(	(iterI.first == "Bd" and iterJ.first == "Ad")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val), 0-pair.bondIJ()));	continue;}
		}
	}
	void addPairingT	(string opt, string svar)	{
		
		/*
		 Bond operation: "Fe:O:+1+0+0 1:2".	(operate for only a single orbital)
		 */
		if( Lat.HSpace() == NORMAL)
			{ ErrorMessage("Error, triplet pairing operation:\n"+opt+".\n Cannot be applited with normal space."); }
		if( Lat.HSpace() == NAMBU and Lat.parameter.STR("spin") == "on")
			{ ErrorMessage("Error, triplet pairing operation:\n"+opt+".\n Cannot be applited with spinful-Nambu space."); }
		
		replaceAll(opt, "\t", " ");
		auto parser = split(opt, " ");
		if( parser.size() != 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		auto firstSec = split(parser[0], ":");
		auto secondSec = split(parser[1], ":");
		if( firstSec.size() != 1 and firstSec.size() != 3 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( secondSec.size() != 2 )	{ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
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
		if( firstSec.size() == 2 ){
			auto xvec = parseSiteVecString(atomI, svar);
			x_mat sVec(1,3);
			if( xvec.size() >= 1 ){ sVec = xvec[0]; }
			for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
			val = sVec[0].real();
		}
		if( firstSec.size() == 3 ){
			auto xvec = parseBondVecString(pair, svar);
			x_mat sVec(1,3);
			if( xvec.size() >= 1 ){ sVec = xvec[0]; }
			for( unsigned i=1 ; i<xvec.size() ; i++){ sVec = sVec * xvec[i][0].real(); }
			val = sVec[0].real();
		}
		
		auto subIndexI = pair.atomI.orbitalIndexList(secondSec[0]);
		auto subIndexJ = pair.atomJ.orbitalIndexList(secondSec[1]);
		for( unsigned ii = 0; ii<subIndexI.size() ; ii++)
		for( unsigned jj = 0; jj<subIndexI.size() ; jj++){
			auto & iterI = subIndexI[ii];
			auto & iterJ = subIndexJ[jj];
			
			if(		(iterI.first == "A"  and iterJ.first == "B" )){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,	val,	pair.bondIJ())); continue;}
			if(		(iterI.first == "B"  and iterJ.first == "A" )){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val),0-pair.bondIJ())); continue;}
			if(		(iterI.first == "Au" and iterJ.first == "Bu")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,	val,	pair.bondIJ())); continue;}
			if(		(iterI.first == "Bu" and iterJ.first == "Au")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val),0-pair.bondIJ())); continue;}
			if(		(iterI.first == "Ad" and iterJ.first == "Bd")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,	val,	pair.bondIJ())); continue;}
			if(		(iterI.first == "Bd" and iterJ.first == "Ad")){ hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val),0-pair.bondIJ())); continue;}
		}
	}
	
	void addScreenCoulomb(string opt, string svar)	{
		replaceAll(opt, "\t", " ");
		
		auto parser = split(opt, " ");
		if( parser.size() > 2)	{ ErrorMessage("Error, not a valid operation:\n screenCoulomb > "+opt); }
		
		double radius = Lat.parameter.VAR("bondRadius", 0).real();
		
		if( parser.size() == 2){
			if( parser[1][0] != '~'){ ErrorMessage("Error, not a valid operation:\n screenCoulomb > "+opt); }
			
			parser[1][0] = ' ';
			double tmpRadius = StrToDouble(parser[1]);
			if( tmpRadius > radius){
				 ErrorMessage("Error, in operation:\n screenCoulomb > "+opt+" > "+svar
							  +"\n the given radius,~"+DoubleToStr(tmpRadius)+", "
							  +"is oarger the 'bondRadius(="+DoubleToStr(radius)+")'.");
			}
			radius = tmpRadius;
		}
		
		auto atomI = Lat.getAtom();
		if( atomI.atomName != parser[0]){ return; } // Name does not match.
		
		removeSpace(svar);
		auto svarParser = split(svar, "*");
		auto orderKey = svarParser[0];
		r_var alpha = 0.0;
		if( svarParser.size() == 2) alpha = parseSiteString(atomI, svarParser[1])[0].real();

		auto rboxAtoms = Lat.getRBox(radius);
		//cout<<alpha<<" "<<atomI.atomName<<" "<<rboxAtoms.second.size()<<endl;
		
		double sumScreenCoulomb = 0;
		for( auto & atomJ: rboxAtoms.second ){
			auto orderJ = order.findOrder(atomJ, orderKey);
			
			double Charge = -Lat.coreCharge.getCharge(atomJ.atomName);
			if( orderJ.first ){ Charge += orderJ.second[0].real(); }
			
			r_mat vecIJ= atomJ.pos - atomI.pos;
			double distIJ = sqrt(cdot(vecIJ,vecIJ));
			
			if( distIJ > 0)
			sumScreenCoulomb += alpha * Charge * exp(-distIJ/radius) / distIJ;
		}
		
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
	
	void addFieldB		(string opt, string svar)	{ // svar === " @:cspin * Jse "
		replaceAll(opt, "\t", " ");
		auto  parser = split(opt, " ");
		if( parser.size() > 2){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		if( parser.size() == 1){ return; }
		
		removeSpace(parser[0]);
		auto dummyParser = split(opt, ":");
		if( dummyParser.size() != 1){ ErrorMessage("Error, not a valid operation:\n"+opt); }
		
		
		auto atomI = Lat.getAtom();
		if( atomI.atomName != parser[0])	return;
		
		
		removeSpace(svar);
		auto svarParser = split(svar, "*");
		auto vectorB = parseSiteString(atomI, svarParser[0]);
		r_var multiplyB = 1.0;
		if( svarParser.size() == 2) multiplyB = parseSiteString(atomI, svarParser[1])[0].real();
		
		vectorB = vectorB * multiplyB;
		
		string orbital = parser[1];
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
	
	/*------------------------------------------------
	 Using these methods to construct and diagonalize (in k-space) Hamiltonian.
	 -------------------------------------------------*/
	void constructTBMHam	()					{
		clear();
		initHam();

		/* Construct the hamElementList from Lat.hamParser.hamOperationList from the 'xxx.lat.tbm'. */
		while( Lat.iterate() ) {
			
			auto atomI = Lat.getAtom();
			
			for( auto & iter : Lat.hamParser.getOperationList("hundSpin")	)	{	addHundSpin(iter.first, iter.second		);	}
			for( auto & iter : Lat.hamParser.getOperationList("orbital")	)	{	addOrbitalEng(iter.first, iter.second	);	}
			for( auto & iter : Lat.hamParser.getOperationList("site")		)	{	addSiteCouple(iter.first, iter.second	);	}
			for( auto & iter : Lat.hamParser.getOperationList("siteHc")		)	{	addSiteCoupleHc(iter.first, iter.second	);	}
			for( auto & iter : Lat.hamParser.getOperationList("hopping")	)	{	addHoppingInt(iter.first, iter.second	);	}
			for( auto & iter : Lat.hamParser.getOperationList("hoppingHc")	)	{	addHoppingIntHc(iter.first, iter.second	);	}
			for( auto & iter : Lat.hamParser.getOperationList("bond")		)	{	addBondCouple(iter.first, iter.second	);	}
			for( auto & iter : Lat.hamParser.getOperationList("bondHc")		)	{	addBondCoupleHc(iter.first, iter.second	);	}
			for( auto & iter : Lat.hamParser.getOperationList("screenCoulomb"))	{	addScreenCoulomb(iter.first, iter.second);	}
			for( auto & iter : Lat.hamParser.getOperationList("fieldB")	)		{	addFieldB(iter.first, iter.second		);	}
			
			
			if( Lat.HSpace() != NORMAL){
				for( auto & iter : Lat.hamParser.getOperationList("pairingS") ){ addPairingS	(iter.first, iter.second	); }
				for( auto & iter : Lat.hamParser.getOperationList("pairingU") ){ addPairingU	(iter.first, iter.second	); }
				for( auto & iter : Lat.hamParser.getOperationList("pairingD") ){ addPairingD	(iter.first, iter.second	); }
				for( auto & iter : Lat.hamParser.getOperationList("pairingT") ){ addPairingT	(iter.first, iter.second	); }
			}
		}
	}
	void addChemicalPotential(r_var mu)			{
		
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
	
	EigVec &	HamEvd			(r_mat kp)		{
		
		Ham.zerolize();
		for( auto & elem: hamElementList){
			
			auto & bvec = elem.bondVector;
			x_var phase = kp[0]*bvec[0]+ kp[1]*bvec[1]+ kp[2]*bvec[2];
			auto	exp_IJ=exp(-Im*phase);
			
			Ham(elem.I, elem.J) += elem.val * exp_IJ;
		}
		
		tmpEigValVec.kPoint = kp;
		tmpEigValVec.status = Ham.evd(tmpEigValVec.eigenValue, tmpEigValVec.eigenVector);
		
		//if( tmpEigValVec.status != "Success" ){
		//	cout<<"at:"<<kp<<" "<<status<<endl;
		//}
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
	
	void		calculate4DensityOrder	()		{
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
					}
					if( parser_I[1] == "d" and parser_J[1] == "d"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += tmpDen;
						r4den[3] -= tmpDen;
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
					}
					if( parser_I[1] == "Ad" and parser_J[1] == "Ad"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += tmpDen*0.5;
						r4den[3] -= tmpDen*0.5;
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
					}                       
					if( parser_I[1] == "Bd" and parser_J[1] == "Bd"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += 1-tmpDen*0.5;
						r4den[3] -= 1-tmpDen*0.5;
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
			
			for( auto & iter: fourDensity){
				order(atomI.atomName+" "+iter.first+":4den") = iter.second;
			}
		}
		order.save();
	}
	void		calculateEnergy			()		{

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
	
	x_mat			parseSiteString(Atom & at, string svar)			{ // @:cspin .. Jh .. 1.0 .. [ 1, 2, 3] ...
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
	vector<x_mat>	parseSiteVecString(Atom & at, string svec)		{ // @:cspin * Jh * ...
		
		removeSpace(svec);
		auto varList = split(svec, "*");
		
		vector<x_mat> xvec;
		for( auto & elem: varList) xvec.push_back( parseSiteString(at,elem));
		
		return xvec;
	}
	x_mat			parseBondString(AtomPair & ap, string svar)		{ // @:cspin .. Jh .. 1.0 .. [ 1, 2, 3] ...
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
	vector<x_mat>	parseBondVecString(AtomPair & ap, string svec)	{ // +1+0+0:1:1:pairS * Jh * ...
		
		removeSpace(svec);
		auto varList = split(svec, "*");
		
		vector<x_mat> xvec;
		for( auto & elem: varList) xvec.push_back(parseBondString(ap, elem));
		
		return xvec;
	}

private:
	void clear			()						{
		Ham.zerolize();
		KEigenValVec.clear();
		hamElementList.clear();
		maxE = -10000;
		minE =  10000;
	}
};







































