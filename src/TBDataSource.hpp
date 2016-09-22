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

struct	EigVec		{
	string	message;
	r_mat	kPoint;
	r_mat	eigenValue;
	x_mat	eigenVector;
};

/*
 This structure will be passed through all future models.
 It can be used to construct all matrix element, for example:
 	orbital		> Fe 1		> 1		%Create the on-site energy for the first orbital of Fe.
 	hundSpin	> Fe 1		> @:cspin * Jh %Create a classical spin (x-direction) coupled to orbital 1.
 
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
	TBDataSource(Lattice & _lat, TBMParser & _tbm):Lat(_lat), order(_lat), order_old(_lat), tbm(_tbm)	{	}
	~TBDataSource(){
		KEigenValVec.clear();
	}
	
	Lattice &				Lat;
	OrderParameter			order;
	OrderParameter			order_old;
	TBMParser	&			tbm;
	
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

	void		addHundSpin		(const PreprocessorInfo & preInfo)	{
		/*****************************************
		* hund > Fe 1 > @:cspin * Jh
		*****************************************/
		auto & optList = preInfo.optList;
		auto & varList = preInfo.varList;
		
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage("Error, 'HundSpin' operation:\n"+preInfo.line+"\nCannot be applied for a spin-independent Hamiltonian.");
		}
		
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
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+
						 " Has been applied with a wrong spin-variable with size: " +IntToStr(sVec.size()));
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
					ErrorMessage(preInfo.filename,
								 preInfo.lineNumber,
								 " \""+preInfo.line+"\"\n"+
								 " Cannot be applied for a spin-depende Nambu space Hamiltonian.");
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
	void		addOrbitalEng	(const PreprocessorInfo & preInfo)	{
		/*****************************************
		 * orbital > Fe 1 > 1
		 *****************************************/
		auto & optList = preInfo.optList;
		auto & varList = preInfo.varList;
		
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

	void		addSiteCouple	(const PreprocessorInfo & preInfo, bool needHc = false)	{
		/*****************************************
		 * site > Fe 1u:2d		> var
		 * site > Fe px.u:py.d	> var
		 *****************************************/
		auto & optList = preInfo.optList;
		auto & varList = preInfo.varList;
		
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+ " Cannot be applied for a spin-independent Hamiltonian.");
		}
		
		if( optList[1][optList[1].size()-1] != optList[2][optList[2].size()-1] and
		   Lat.HSpace() == NAMBU) {
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+ " Spin-flip term cannot be applied under 'nambu' space.");
		}
		
		auto atomI = Lat.getAtom();
		if( optList[0] != atomI.atomName )	return;
		

		auto spinIndexList_i = atomI.spinIndexList(optList[1]);
		auto spinIndexList_j = atomI.spinIndexList(optList[2]);
		
		if(spinIndexList_i.size() <= 0 or spinIndexList_j.size() <= 0){
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+ "atom \'"+atomI.atomName+"\' don't have any orbital to operate..");
		}
		
		x_var val = parseSiteString(atomI, varList[0])[0];
		for( unsigned i=1 ; i<varList.size() ; i++){
			val = val * parseSiteString(atomI, varList[i])[0];
		}
		
		for( unsigned ii=0 ; ii<spinIndexList_i.size() ; ii++){
			auto & tmpChar = spinIndexList_i[ii].first[0];
			
			auto & index_i = spinIndexList_i[ii].second;
			auto & index_j = spinIndexList_j[ii].second;
			
			if(tmpChar == 'n' or tmpChar == 'u' or tmpChar == 'd'){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, vec(0,0,0)));
				if(needHc)
					hamElementList.push_back(MatrixElement(index_j, index_i,conj(val), vec(0,0,0)));
			}
			if(tmpChar == 'A'){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, vec(0,0,0)));
				if(needHc)
					hamElementList.push_back(MatrixElement(index_j, index_i,conj(val), vec(0,0,0)));
			}
			if(tmpChar == 'B'){
				hamElementList.push_back(MatrixElement(index_i, index_j,-val, vec(0,0,0)));
				if(needHc)
					hamElementList.push_back(MatrixElement(index_j, index_i,-conj(val), vec(0,0,0)));
			}
		}
	}
	void		addHoppingInt	(const PreprocessorInfo & preInfo, bool needHc = false)	{
		/*****************************************
		 * hopping > Fe:O:+1+0+0# 1:2 >  0.5
		 * hopping > Fe:O:+1+0+0# 1:px >  0.5
		 *****************************************/
		
		auto & optList = preInfo.optList;
		auto & varList = preInfo.varList;
		
		auto pair = Lat.getPair(optList[2]);
		if( pair.atomI.atomName != optList[0])	return;
		if( pair.atomJ.atomName != optList[1])	return;
		if( !pair.withinRange() )				return;
		
		if(!pair.atomI.hasOrbital(optList[3]))return;
		if(!pair.atomJ.hasOrbital(optList[4]))return;
		
		x_mat sVec = parseBondString(pair, varList[0]);
		x_var val = sVec[0];
		for( unsigned i=1 ; i<varList.size() ; i++){
			val = val * parseBondString(pair, varList[i])[0];
		}
		
		auto subIndexI = pair.atomI.orbitalIndexList(optList[3]);
		auto subIndexJ = pair.atomJ.orbitalIndexList(optList[4]);
		for( unsigned ii = 0; ii<subIndexI.size() ; ii++){
			auto & iterI = subIndexI[ii];
			auto & iterJ = subIndexJ[ii];
			char & tmpChar = iterI.first[0];
			if(tmpChar == 'n' or tmpChar == 'u' or tmpChar == 'd'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));
				if(needHc)
					hamElementList.push_back(MatrixElement(iterJ.second, iterI.second, conj(val),	-pair.bondIJ()));
			}
			if(tmpChar == 'A'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));
				if(needHc)
					hamElementList.push_back(MatrixElement(iterJ.second, iterI.second, conj(val),	-pair.bondIJ()));
			}
			if(tmpChar == 'B'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,-val, pair.bondIJ()));
				if(needHc)
					hamElementList.push_back(MatrixElement(iterJ.second, iterI.second,conj(-val),	-pair.bondIJ()));
			}
		}
	}
	void		addBondCouple	(const PreprocessorInfo & preInfo, bool needHc = false)	{
		/*****************************************
		 * bond > Fe:O:+1+0+0# 1:1 >  0.5
		 * bond > Fe:O:+1+0+0# 1:px >  0.5
		 *****************************************/
		auto & optList = preInfo.optList;
		auto & varList = preInfo.varList;
		
		auto pair = Lat.getPair(optList[2]);
		if( pair.atomI.atomName != optList[0])	return;
		if( pair.atomJ.atomName != optList[1])	return;
		if( !pair.withinRange() )				return;
		
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+ " Cannot be applied for a spin-independent Hamiltonian.");
		}
		
		if( optList[3][optList[3].size()-1] != optList[4][optList[4].size()-1] and
		   Lat.HSpace() == NAMBU) {
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+ " Spin-flip term cannot be applied under 'nambu' space.");
		}
		
		x_mat sVec = parseBondString(pair, varList[0]);
		x_var val = sVec[0];
		for( unsigned i=1 ; i<varList.size() ; i++){
			val = val * parseBondString(pair, varList[i])[0];
		}
		
		auto spinIndexList_i = pair.atomI.spinIndexList(optList[3]);
		auto spinIndexList_j = pair.atomJ.spinIndexList(optList[4]);
		
		if(spinIndexList_i.size() <= 0){
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+ "atom \'"+pair.atomI.atomName+"\' don't have such orbital to operate.");
		}
		if(spinIndexList_j.size() <= 0){
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+ "atom \'"+pair.atomJ.atomName+"\' don't have such orbital to operate.");
		}
		
		for( unsigned ii=0 ; ii<spinIndexList_i.size() ; ii++){
			auto & tmpChar = spinIndexList_i[ii].first[0];
			
			auto & index_i = spinIndexList_i[ii].second;
			auto & index_j = spinIndexList_j[ii].second;
			
			if(tmpChar == 'n' or tmpChar == 'u' or tmpChar == 'd'){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, pair.bondIJ()));
				if( needHc )
					hamElementList.push_back(MatrixElement(index_j, index_i, conj(val), -pair.bondIJ()));
			}
			if(tmpChar == 'A'){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, pair.bondIJ()));
				if( needHc )
					hamElementList.push_back(MatrixElement(index_j, index_i, conj(val), -pair.bondIJ()));
			}
			if(tmpChar == 'B'){
				hamElementList.push_back(MatrixElement(index_i, index_j,-val, pair.bondIJ()));
				if( needHc )
					hamElementList.push_back(MatrixElement(index_j, index_i, conj(-val),-pair.bondIJ()));
			}
		}
	}
	
	void		addScreenCoulomb(const PreprocessorInfo & preInfo)	{
		/*****************************************
		 * screenCoulomb > Fe ~2 >  @:den * alpha
		 *****************************************/
		auto & optList = preInfo.optList;
		auto & varList = preInfo.varList;
		
		double radius = Lat.parameter.VAR("bondRadius", 0).real();
		
		string coulombRadiusStr = optList[1];
		coulombRadiusStr[0] = ' ';
		double tmpRadius = StrToDouble(coulombRadiusStr);
		if( tmpRadius > radius){
			
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+
						 " the given radius "+ DoubleToStr(tmpRadius) +
						 "is larger the 'bondRadius(="+DoubleToStr(radius)+")'.");
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
		
		auto orderI = order.findOrder(atomI, orderKey);
		auto Z_I = tbm.coreCharge.getCharge(atomI.atomName);
		double ChargeI = orderI.second[0].real() - Z_I;
		
		double sumScreenCoulomb = 0;
		for( auto & atomJ: rboxAtoms.second ){
			auto orderJ = order.findOrder(atomJ, orderKey);
			auto Z_J = tbm.coreCharge.getCharge(atomJ.atomName);
			
			double ChargeJ = -Z_J;
			if( orderJ.first ){ ChargeJ += orderJ.second[0].real(); }
			
			r_mat vecIJ(1,3);
			for(unsigned i=0 ; i<vecIJ.size() ; i++)
				vecIJ[i] = atomJ.pos[i] - atomI.pos[i];
			
			double distIJ = sqrt(cdot(vecIJ,vecIJ));
			
			if( distIJ > 0){
				double factor = alpha * exp(-distIJ/radius) / distIJ;
				sumScreenCoulomb		+= factor * ChargeJ;
				energyMap["2.Coul Eng"] -= factor * orderI.second[0].real() * orderJ.second[0].real();
				energyMap["2.Coul Eng"] += factor * Z_I * Z_J;
			}
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
	
	void		addFieldB		(const PreprocessorInfo & preInfo)	{ // svar === " @:cspin * Jse "
		// *******************************
		// fieldB	> Fe   > [0,0,1] * B
		// fieldB	> Fe 1 > [0,0,1] * B
		// fieldB	> Fe 2 > [0,0,1] * B
		// *******************************
		auto & optList = preInfo.optList;
		auto & varList = preInfo.varList;
		
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
					ErrorMessage(preInfo.filename, preInfo.lineNumber,
								 " \""+preInfo.line+"\"\n"+
								 " Cannot be applied for a spin-depende Nambu space Hamiltonian.");
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

	void		addPairing		(const PreprocessorInfo & preInfo, string pairingType)	{
		/*****************************************
		 * pairing(X)	> Fe 1:1			> 1		% On-site Singlet pairing
		 * pairing(X)	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site Singlet pairing
		 *****************************************/
		auto & optList = preInfo.optList;
		auto & varList = preInfo.varList;
		
		if(pairingType == "S")
		if( Lat.parameter.STR("spin") == "off" and Lat.HSpace() == NAMBU){
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+
						 " Singlet pairing cannot be applited with spinless-nambu space.");
		}
		if(pairingType == "U")
		if( Lat.HSpace() != EXNAMBU){
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+
						 " Up-triplet pairing cannot be applited without exnambu space.");
		}
		if(pairingType == "D")
		if( Lat.HSpace() == EXNAMBU){
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+
						 " Dn-triplet pairing cannot be applited without exnambu space.");
		}
		if(pairingType == "T")
		if( Lat.HSpace() == EXNAMBU){
			ErrorMessage(preInfo.filename, preInfo.lineNumber,
						 " \""+preInfo.line+"\"\n"+
						 " Triplet pairing cannot be applited without exnambu space.");
		}
		
		
		deque<string> firstSec;
		deque<string> secondSec;
		if( optList.size() == 3){
			firstSec.push_back(optList[0]);
			secondSec.push_back(optList[1]);
			secondSec.push_back(optList[2]);
		}
		else if( optList.size() == 5){
			firstSec.push_back(optList[0]);
			firstSec.push_back(optList[1]);
			firstSec.push_back(optList[2]);
			secondSec.push_back(optList[3]);
			secondSec.push_back(optList[4]);
		}
		
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
			
			if( pairingType == "S"){
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
			if( pairingType == "U" or pairingType == "T"){
				if(		(iterI.first == "Au" and iterJ.first == "Bu")	){
					hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));
					continue;
				}
				if(		(iterI.first == "Bu" and iterJ.first == "Au")	){
					hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val), -pair.bondIJ()));
					continue;
				}
			}
			if( pairingType == "D" or pairingType == "T"){
				if(		(iterI.first == "Ad" and iterJ.first == "Bd")	){
					hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));
					continue;
				}
				if(		(iterI.first == "Bd" and iterJ.first == "Ad")	){
					hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val), -pair.bondIJ()));
					continue;
				}
			}
			if( pairingType == "T" ){
				if(		(iterI.first == "A"  and iterJ.first == "B" )){
					hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,	val,	pair.bondIJ()));
					continue;
				}
				if(		(iterI.first == "B"  and iterJ.first == "A" )){
					hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(val),-pair.bondIJ()));
					continue;
				}
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
	void			constructTBMHam	(bool withMu = true)		{
		
		initHam();
		
		//if( hamElementList.size() > 0 ) return;
		
		maxE = -10000;
		minE =  10000;
		
		energyMap["2.Coul Eng"] = 0;
		
		/*
		 // Performance testing loop
		 unsigned iter = 0;
		 cout<<"Entering debugging view"<<endl;
		 
		 auto  pair = Lat.getPair(0, vec(1,0,0)); //***
		 while(true){
		 	auto ddd = split("AAA,BBB,CCC,EEE", ",");
		 	hamElementList.push_back(MatrixElement(0,0,0,vec(0,0,0)));
		 	iter++;
		 	cout<<iter<<" "<<hamElementList.size()<<endl;
		 	if( iter % 2000 == 0){
		 		hamElementList.clear();
		 	}
		 	while( Lat.iterate() ){
		 		auto atom = Lat.getAtom();
		 		x_mat	sVec = parseSiteString(atom, "0");
		 		Lat.vec("+1+0+0#"); //---
		 		Lat.vec("+1+0+0"); //---
		 		auto  pair = Lat.getPair(vec(1,0,0)); //***
		 		x_mat tmpvec = parseBondString(pair, "tdp"); //---
		 	}
		 }
		 cout<<"Leaving debugging view"<<endl;
		*/
		
		/* Construct the hamElementList from Lat.hamParser.hamOperationList from the 'xxx.lat.tbm'. */
		hamElementList.clear();
		while( Lat.iterate() ) {
			for( auto & iter : tbm.hamPreprocessor.list_OrbitalEng)		{	addOrbitalEng	(iter);}
			for( auto & iter : tbm.hamPreprocessor.list_HundSpin)		{	addHundSpin		(iter);}
			
			for( auto & iter : tbm.hamPreprocessor.list_SiteCouple)		{	addSiteCouple	(iter);}
			for( auto & iter : tbm.hamPreprocessor.list_SiteCoupleHc)	{	addSiteCouple	(iter,true);}
			for( auto & iter : tbm.hamPreprocessor.list_HoppingInt)		{	addHoppingInt	(iter);}
			for( auto & iter : tbm.hamPreprocessor.list_HoppingIntHc)	{	addHoppingInt	(iter,true);}
			for( auto & iter : tbm.hamPreprocessor.list_BondCouple)		{	addBondCouple	(iter);}
			for( auto & iter : tbm.hamPreprocessor.list_BondCoupleHc)	{	addBondCouple	(iter,true);}
			for( auto & iter : tbm.hamPreprocessor.list_ScreenCoulomb)	{	addScreenCoulomb(iter);}
			for( auto & iter : tbm.hamPreprocessor.list_QuantumFieldB)	{	addFieldB		(iter);}
			
			// If space==normal The following part will be ignored.
			if( Lat.HSpace() == NORMAL) continue;
			
			for( auto & iter : tbm.hamPreprocessor.list_PairingS)		{	addPairing		(iter,"S"); }
			for( auto & iter : tbm.hamPreprocessor.list_PairingU)		{	addPairing		(iter,"U"); }
			for( auto & iter : tbm.hamPreprocessor.list_PairingD)		{	addPairing		(iter,"D"); }
			for( auto & iter : tbm.hamPreprocessor.list_PairingT)		{	addPairing		(iter,"T"); }
			
		}
		if( withMu ) addChemicalPotential(Lat.parameter.VAR("Mu").real());
	}

	EigVec	&		HamEvd			(r_mat kp)					{
		
		Ham.zerolize();
		
		//cout<<"KKK"<<kp<<endl;
		for( auto & elem: hamElementList){
			
			auto bvec = elem.bondVector;
			x_var phase = kp[0]*bvec[0]+ kp[1]*bvec[1]+ kp[2]*bvec[2];
			auto	exp_IJ=exp(-Im*phase);
			
			Ham(elem.I, elem.J) += elem.val * exp_IJ;
		}
		//cout<<Ham<<endl<<endl;
		//for( unsigned i=0 ; i<Ham.cols() ; i++)
		//for( unsigned j=0 ; j<Ham.cols() ; j++){
		//	if( abs(Ham(i,j)) > 0.0001 )
		//	cout<<fformat(i)<<fformat(j)<<" "<<Ham(i,j)<<endl;
		//}
		
		//EigVec					tmpEigValVec;	// A temporary storage of eigen value vector for a specific K-point.
		tmpEigValVec.kPoint = kp;
		
		tmpEigValVec.message = Ham.evd(tmpEigValVec.eigenValue, tmpEigValVec.eigenVector);
		
		if ( maxE < tmpEigValVec.eigenValue[Lat.indexSize()-1]	){ maxE = tmpEigValVec.eigenValue[Lat.indexSize()-1]; }
		if ( minE > tmpEigValVec.eigenValue[0]					){ minE = tmpEigValVec.eigenValue[0]; }
		
		return tmpEigValVec;
	}
	
	x_var			getDensityMatrix		(unsigned index_i, unsigned index_j, double Mu=0){
		
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
	
	void			calculate4DensityOrder	()						{
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
	void			calculateEnergy			()						{

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
	
	x_mat			parseSiteString(Atom & at, string svar)			{
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
	
	x_mat			parseBondString(AtomPair & ap, string svar)		{
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



