/*-----------------------------------------------------------|
| Copyright (C) 2016 Yuan-Yen Tai, Hongchul Choi,            |
|                    Jian-Xin Zhu                            |
|                                                            |
| This file is distributed under the terms of the BSD        |
| Berkeley Software Distribution. See the file `LICENSE' in  |
| the root directory of the present distribution.            |
|                                                            |
|-----------------------------------------------------------*/
//
//  TBDataSource.hpp
//  TBM^3
//
//  Created by Yuan-Yen Tai on 7/19/16.
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

struct PairingElement{
public:
	AtomPair	atomPair;
	string		orderStr, orderKey;
	double		V;
	
	typedef map<string,unsigned > PairIndexMap;
	PairIndexMap pIndexMap;
	
	PairingElement(){}
	PairingElement(AtomPair ap, double _v, string _orderStr, string _orderKey, PairIndexMap _indexMap){
		V			= _v;
		atomPair	= ap;
		orderStr	= _orderStr;
		orderKey	= _orderKey;
		pIndexMap	= _indexMap;
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
	pairingUU	> Fe 1:1			> 1		% On-site Triplet pairing for 'up' spin
	pairingUU	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site Triplet pairing for 'up' spin
	pairingDD	> Fe 1:1			> 1		% On-site Triplet pairing for 'dn' spin
	pairingDD	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site Triplet pairing for 'dn' spin
	pairingUD	> Fe 1:1			> 1		% On-site Triplet pairing
	pairingUD	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site Triplet pairing
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

	vector<PairingElement>	pairElementList;	// A list to store all the matrix elements

	map<string, double>		energyMap;

	void		initHam		()		{
		Ham = x_mat(Lat.indexSize(), Lat.indexSize());
		initOrder();
	}
	void		initOrder	()		{
		order.clear();
		order.load();
	}

	/*------------------------------------------------
	 Using these methods to construct the "hamElementList".
	 -------------------------------------------------*/
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
				
				hamElementList.push_back(MatrixElement(indexBu, indexBd,-Sx-Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBu,-Sx+Im*Sy	, vec(0,0,0)));
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
		
		// Stretch the spin character of I and J part.
		char spinCharI = optList[1][optList[1].size()-1];
		char spinCharJ = optList[2][optList[2].size()-1];
		
		// Get a list of all the operable spin-dependent indexList.
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
			auto & tmpCharI = spinIndexList_i[ii].first[0];
			
			auto & index_i = spinIndexList_i[ii].second;
			auto & index_j = spinIndexList_j[ii].second;
			
			if(tmpCharI == 'n' or tmpCharI == 'u' or tmpCharI == 'd'){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, vec(0,0,0)));
				if(needHc)
					hamElementList.push_back(MatrixElement(index_j, index_i,conj(val), vec(0,0,0)));
			}
			if(tmpCharI == 'A'){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, vec(0,0,0)));
				if(needHc)
					hamElementList.push_back(MatrixElement(index_j, index_i,conj(val), vec(0,0,0)));
			}
			if(tmpCharI == 'B'){
				if( spinCharI == spinCharJ ){
					hamElementList.push_back(MatrixElement(index_i, index_j,conj(-val), vec(0,0,0)));
					if(needHc)
						hamElementList.push_back(MatrixElement(index_j, index_i,-val, vec(0,0,0)));
				}
				else{
					hamElementList.push_back(MatrixElement(index_i, index_j,conj(val), vec(0,0,0)));
					if(needHc)
						hamElementList.push_back(MatrixElement(index_j, index_i,val, vec(0,0,0)));
				}
			}
		}
	}
	void		addHoppingInt	(const PreprocessorInfo & preInfo, bool needHc = false)	{
		/*****************************************
		 * hopping > Fe:O:+1+0+0# 1:2 >  0.5
		 * hopping > Fe:O:+1+0+0# 1:px >  0.5
		 *****************************************/
		// Note that, there is no spin-flip term in the regular hopping term.
		
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
			char & tmpCharI = iterI.first[0];
			if(tmpCharI == 'n' or tmpCharI == 'u' or tmpCharI == 'd'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));
				if(needHc)
					hamElementList.push_back(MatrixElement(iterJ.second, iterI.second, conj(val),	-pair.bondIJ()));
			}
			if(tmpCharI == 'A'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second, val, pair.bondIJ()));
				if(needHc)
					hamElementList.push_back(MatrixElement(iterJ.second, iterI.second, conj(val),	-pair.bondIJ()));
			}
			if(tmpCharI == 'B'){
				hamElementList.push_back(MatrixElement(iterI.second, iterJ.second,conj(-val), pair.bondIJ()));
				if(needHc)
					hamElementList.push_back(MatrixElement(iterJ.second, iterI.second,-val,	-pair.bondIJ()));
			}
		}
	}
	void		addBondCouple	(const PreprocessorInfo & preInfo, bool needHc = false)	{
		/*****************************************
		 * bond > Fe:O:+1+0+0# 1u:1d >  0.5
		 * bond > Fe:O:+1+0+0# 1u:px.d >  0.5
		 *****************************************/
		// Note that, the spin-flip term is applicable in the bond term.
		
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
		
		// Stretch the spin character of I and J part.
		char spinCharI = optList[3][optList[3].size()-1];
		char spinCharJ = optList[4][optList[4].size()-1];
		
		// Get a list of all the operable spin-dependent indexList.
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
			auto & tmpCharI = spinIndexList_i[ii].first[0];
			
			auto & index_i = spinIndexList_i[ii].second;
			auto & index_j = spinIndexList_j[ii].second;
			
			if(tmpCharI == 'n' or tmpCharI == 'u' or tmpCharI == 'd'){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, pair.bondIJ()));
				if( needHc )
					hamElementList.push_back(MatrixElement(index_j, index_i, conj(val), -pair.bondIJ()));
			}
			if(tmpCharI == 'A'){
				hamElementList.push_back(MatrixElement(index_i, index_j, val, pair.bondIJ()));
				if( needHc )
					hamElementList.push_back(MatrixElement(index_j, index_i, conj(val), -pair.bondIJ()));
			}
			if(tmpCharI == 'B'){
				if( spinCharI == spinCharJ ){
					hamElementList.push_back(MatrixElement(index_i, index_j,conj(-val), pair.bondIJ()));
					if( needHc )
						hamElementList.push_back(MatrixElement(index_j, index_i, -val,-pair.bondIJ()));
				}
				else{
					hamElementList.push_back(MatrixElement(index_i, index_j,conj(val), pair.bondIJ()));
					if( needHc )
						hamElementList.push_back(MatrixElement(index_j, index_i, val,-pair.bondIJ()));
				}
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
				
				hamElementList.push_back(MatrixElement(indexBu, indexBd,-Sx-Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBu,-Sx+Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBu, indexBu,-Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBd, Sz			, vec(0,0,0)));
				break;
			}
		}
	}
	void		addPairing		(const PreprocessorInfo & preInfo, string pairingType)	{
		/*****************************************
		pairingX	> Fe 1:1			> 1		% On-site  intra-orbital X-pairing
		pairingX	> Fe:Fe:+1+0+0 1:1	> 1		% Off-site intra-orbital X-pairing
		pairingX	> Fe 1:2			> 1		% On-site  inter-orbital X-pairing
		pairingX	> Fe:Fe:+1+0+0 1:2	> 1		% Off-site inter-orbital X-pairing
		X = S  --> Spin Singlet pairing
		X = UU --> Up-Up spin	triplet pairing
		X = DD --> Dn-Dn Spin	triplet pairing
		X = UD --> UP-Dn Spin	triplet pairing
		 *****************************************/
		auto & optList = preInfo.optList;
		auto & varList = preInfo.varList;
		
		string orderStr	= "";
		if(pairingType == "S"){
			if( Lat.parameter.STR("spin") == "off" and Lat.HSpace() == NAMBU){
				ErrorMessage(preInfo.filename, preInfo.lineNumber,
							 " \""+preInfo.line+"\"\n"+
							 " Singlet pairing cannot be applited with spinless-nambu space.");
			}
			orderStr = "spair";
		}
		if(pairingType == "UD"){
			if( Lat.parameter.STR("spin") == "off" and Lat.HSpace() == NAMBU){
				ErrorMessage(preInfo.filename, preInfo.lineNumber,
							 " \""+preInfo.line+"\"\n"+
							 " UD-Triplet pairing cannot be applited with spinless-nambu space.");
			}
			orderStr = "udpair";
		}
		if(pairingType == "UU"){
			if( Lat.HSpace() != EXNAMBU){
				ErrorMessage(preInfo.filename, preInfo.lineNumber,
							 " \""+preInfo.line+"\"\n"+
							 " UU-triplet pairing cannot be applited without exnambu space.");
			}
			orderStr = "uupair";
		}
		if(pairingType == "DD"){
			if( Lat.HSpace() != EXNAMBU){
				ErrorMessage(preInfo.filename, preInfo.lineNumber,
							 " \""+preInfo.line+"\"\n"+
							 " DD-triplet pairing cannot be applited without exnambu space.");
			}
			orderStr = "ddpair";
		}


		string orderKey = "";
		deque<string> firstSec;
		deque<string> secondSec;
		if( optList.size() == 3){
			firstSec.push_back(optList[0]);
			secondSec.push_back(optList[1]);
			secondSec.push_back(optList[2]);
			orderKey = "@:"+optList[1]+":"+optList[2]+":"+orderStr;
		}
		else if( optList.size() == 5){
			firstSec.push_back(optList[0]);
			firstSec.push_back(optList[1]);
			firstSec.push_back(optList[2]);
			secondSec.push_back(optList[3]);
			secondSec.push_back(optList[4]);
			r_mat orderVec = Lat.vec(optList[2]);
			orderKey = vecToStr(orderVec)+":"+optList[3]+":"+optList[4]+":"+orderStr;
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
		
		string orb1 = secondSec[0];
		string orb2 = secondSec[1];
		
		// The pairPoetnetial will sent into pairElementList for further calculation.
		r_var pairPotential = parseSiteString(atomI, varList[0])[0].real();
		for( unsigned i=1 ; i<varList.size() ; i++){
			pairPotential = pairPotential * parseSiteString(atomI, varList[i])[0].real();
		}
		pairPotential *= -1;
		
		x_var val = 0;
		auto orderI = order.findOrder(atomI, orderKey);
		if( orderI.first ){
			val = orderI.second[0];
		}
		else{
			// If the order does not exist inside *.lat.ord,
			// assign a pair term to it and give a selection of random phase,
			//
			// if( pair_phase == "random" ) --> the random phase will be assigned.
			// else --> no phase.
			//
			// Note that, in the future we could design more sophisticate phase assignment.
			val = tbm.parameter.VAR("init_pair_magnitude", 1.0);
			
			string pair_phase_flag = tbm.parameter.STR("init_pair_phase", "none");
			if( pair_phase_flag == "random"){
				double phaseTheta = (double)rand()/RAND_MAX;
				x_var phase = exp(Im * 2 * pi * phaseTheta);
				val = val * phase * pairPotential;
			}
		}
		
		energyMap["5.Pair Eng"] += (conj(val)*val).real() / pairPotential;
		
		switch( Lat.HSpace() ){
			case NAMBU: {
				if( pairingType == "S" ){
					
					unsigned indexI1Au = pair.atomI.index(orb1+"Au");
					unsigned indexJ2Bd = pair.atomJ.index(orb2+"Bd");
					hamElementList.push_back(MatrixElement(indexI1Au, indexJ2Bd, val		, pair.bondIJ()));
					
					unsigned indexI1Bd = pair.atomI.index(orb1+"Bd"); //Conjugate part
					unsigned indexJ2Au = pair.atomJ.index(orb2+"Au"); //Conjugate part
					hamElementList.push_back(MatrixElement(indexI1Bd, indexJ2Au, conj(val)	, pair.bondIJ()));
					
					PairingElement::PairIndexMap pIndexMap;
					pIndexMap["iAu"]=indexI1Au;
					pIndexMap["jBd"]=indexJ2Bd;
					pIndexMap["iBd"]=indexI1Bd;
					pIndexMap["jAu"]=indexJ2Au;
					pairElementList.push_back(PairingElement(pair, pairPotential, orderStr, orderKey, pIndexMap));
					break;
				}
			}
			case EXNAMBU: {
				if( pairingType == "S"){
					
					unsigned indexI1Au = pair.atomI.index(orb1+"Au");
					unsigned indexJ2Bd = pair.atomJ.index(orb2+"Bd");
					hamElementList.push_back(MatrixElement(indexI1Au, indexJ2Bd, val		, pair.bondIJ()));
					
					unsigned indexI1Ad = pair.atomI.index(orb1+"Ad");
					unsigned indexJ2Bu = pair.atomJ.index(orb2+"Bu");
					hamElementList.push_back(MatrixElement(indexI1Ad, indexJ2Bu, val		, pair.bondIJ()));
					
					unsigned indexI1Bu = pair.atomI.index(orb1+"Bu"); //Conjugate part
					unsigned indexJ2Ad = pair.atomJ.index(orb2+"Ad"); //Conjugate part
					hamElementList.push_back(MatrixElement(indexI1Bu, indexJ2Ad, conj(val)	, pair.bondIJ()));
					
					unsigned indexI1Bd = pair.atomI.index(orb1+"Bd"); //Conjugate part
					unsigned indexJ2Au = pair.atomJ.index(orb2+"Au"); //Conjugate part
					hamElementList.push_back(MatrixElement(indexI1Bd, indexJ2Au, conj(val)	, pair.bondIJ()));
					
					PairingElement::PairIndexMap pIndexMap;
					pIndexMap["iAu"]=indexI1Au;
					pIndexMap["iAd"]=indexI1Ad;
					pIndexMap["iBu"]=indexI1Bu;
					pIndexMap["iBd"]=indexI1Bd;
					pIndexMap["jAu"]=indexJ2Au;
					pIndexMap["jAd"]=indexJ2Ad;
					pIndexMap["jBu"]=indexJ2Bu;
					pIndexMap["jBd"]=indexJ2Bd;
					pairElementList.push_back(PairingElement(pair, pairPotential, orderStr, orderKey, pIndexMap));
					break;
				}
				if( pairingType == "UD"){
					
					unsigned indexI1Au = pair.atomI.index(orb1+"Au");
					unsigned indexJ2Bd = pair.atomJ.index(orb2+"Bd");
					hamElementList.push_back(MatrixElement(indexI1Au, indexJ2Bd, val		, pair.bondIJ()));
					
					unsigned indexI1Ad = pair.atomI.index(orb1+"Ad");
					unsigned indexJ2Bu = pair.atomJ.index(orb2+"Bu");
					hamElementList.push_back(MatrixElement(indexI1Ad, indexJ2Bu,-val		, pair.bondIJ()));
					
					unsigned indexI1Bu = pair.atomI.index(orb1+"Bu"); //Conjugate part
					unsigned indexJ2Ad = pair.atomJ.index(orb2+"Ad"); //Conjugate part
					hamElementList.push_back(MatrixElement(indexI1Bu, indexJ2Ad, conj(-val)	, pair.bondIJ()));
					
					unsigned indexI1Bd = pair.atomI.index(orb1+"Bd"); //Conjugate part
					unsigned indexJ2Au = pair.atomJ.index(orb2+"Au"); //Conjugate part
					hamElementList.push_back(MatrixElement(indexI1Bd, indexJ2Au, conj(val)	, pair.bondIJ()));
					
					PairingElement::PairIndexMap pIndexMap;
					pIndexMap["iAu"]=indexI1Au;
					pIndexMap["iAd"]=indexI1Ad;
					pIndexMap["iBu"]=indexI1Bu;
					pIndexMap["iBd"]=indexI1Bd;
					pIndexMap["jAu"]=indexJ2Au;
					pIndexMap["jAd"]=indexJ2Ad;
					pIndexMap["jBu"]=indexJ2Bu;
					pIndexMap["jBd"]=indexJ2Bd;
					pairElementList.push_back(PairingElement(pair, pairPotential, orderStr, orderKey, pIndexMap));
					break;
				}
				if( pairingType == "UU" ){
					if( pair.atomI.atomIndex == pair.atomJ.atomIndex and orb1 == orb2 )
						ErrorMessage(preInfo.filename, preInfo.lineNumber,
							 " \""+preInfo.line+"\"\n"+
							 " pairingUU (UU-triplet pairing) cannot be applited for on-site intra-orbital.");
					
					unsigned indexI1Au = pair.atomI.index(orb1+"Au");
					unsigned indexJ2Bu = pair.atomJ.index(orb2+"Bu");
					hamElementList.push_back(MatrixElement(indexI1Au, indexJ2Bu,-val		, pair.bondIJ()));
					
					unsigned indexI1Bu = pair.atomI.index(orb1+"Bu"); //Conjugate part
					unsigned indexJ2Au = pair.atomJ.index(orb2+"Au"); //Conjugate part
					hamElementList.push_back(MatrixElement(indexI1Bu, indexJ2Au, conj(-val)	, pair.bondIJ()));
					
					PairingElement::PairIndexMap pIndexMap;
					pIndexMap["iAu"]=indexI1Au;
					pIndexMap["iBu"]=indexI1Bu;
					pIndexMap["jAu"]=indexJ2Au;
					pIndexMap["jBu"]=indexJ2Bu;
					pairElementList.push_back(PairingElement(pair, pairPotential, orderStr, orderKey, pIndexMap));
					break;
				}
				if( pairingType == "DD" ){
					if( pair.atomI.atomIndex == pair.atomJ.atomIndex and orb1 == orb2 )
						ErrorMessage(preInfo.filename, preInfo.lineNumber,
							 " \""+preInfo.line+"\"\n"+
							 " pairingDD (DD-triplet pairing) cannot be applited for on-site intra-orbital.");
	
					
					unsigned indexI1Ad = pair.atomI.index(orb1+"Ad");
					unsigned indexJ2Bd = pair.atomJ.index(orb2+"Bd");
					hamElementList.push_back(MatrixElement(indexI1Ad, indexJ2Bd, val		, pair.bondIJ()));
					
					unsigned indexI1Bd = pair.atomI.index(orb1+"Bd"); //Conjugate part
					unsigned indexJ2Ad = pair.atomJ.index(orb2+"Ad"); //Conjugate part
					hamElementList.push_back(MatrixElement(indexI1Bd, indexJ2Ad, conj(val)	, pair.bondIJ()));
					
					PairingElement::PairIndexMap pIndexMap;
					pIndexMap["iAd"]=indexI1Ad;
					pIndexMap["iBd"]=indexI1Bd;
					pIndexMap["jAd"]=indexJ2Ad;
					pIndexMap["jBd"]=indexJ2Bd;
					pairElementList.push_back(PairingElement(pair, pairPotential, orderStr, orderKey, pIndexMap));
					break;
				}
			}
		}
	}

	void		addIntraHubbard	(const PreprocessorInfo & preInfo)	{
		/*****************************************
		 * intraHubbard > Fe 2 >  U
		 *****************************************/
		auto & optList = preInfo.optList;
		auto & varList = preInfo.varList;
		
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage(preInfo.filename, preInfo.lineNumber, " \""+preInfo.line+"\"\n"+
						 " Hubbard model cannot be constructed under spin-independent space.");
		}
		if( Lat.HSpace() == NAMBU ){
			ErrorMessage(preInfo.filename, preInfo.lineNumber, " \""+preInfo.line+"\"\n"+
						 " a model for 3-dimensional spin is not applicable with NAMBU space.");
		}
		
		auto atomI = Lat.getAtom();
		if( atomI.atomName != optList[0]){ return; } // Name does not match.
		if( !atomI.hasOrbital(optList[1])){
			ErrorMessage(preInfo.filename, preInfo.lineNumber, " \""+preInfo.line+"\"\n"+
						 " The atom, \'"+atomI.atomName+"\', has no orbital, "+optList[1]+".");
		}

		string orbitalNumber = atomI.getOrbitalNumber(optList[1]);
		string orderKey = "@:"+orbitalNumber+":4den";
		
		auto orderI = order.findOrder(atomI, orderKey);
		if( !orderI.first ){
			ErrorMessage(preInfo.filename, preInfo.lineNumber, " \""+preInfo.line+"\"\n"+
						 " The 4 density order parameter,"+orderKey+" , was not defined.");
		}
		
		r_var hubbardU = parseSiteString(atomI, varList[0])[0].real();
		for( unsigned i=1 ; i<varList.size() ; i++){
			hubbardU = hubbardU * parseSiteString(atomI, varList[i])[0].real();
		}
		
		x_var	cI_uu = hubbardU * 0.5 * (orderI.second[0].real() + orderI.second[3].real());
		x_var	cI_dd = hubbardU * 0.5 * (orderI.second[0].real() - orderI.second[3].real());
		x_var	cI_ud = hubbardU * 0.5 * (orderI.second[1].real() - Im * orderI.second[2].real());
		x_var	cI_du = hubbardU * 0.5 * (orderI.second[1].real() + Im * orderI.second[2].real());
		
		energyMap["3.U Eng"] -= (cI_uu * cI_dd / hubbardU).real();
		energyMap["3.U Eng"] += (cI_ud * cI_du / hubbardU).real();

		////*****************************************************************
		//// Attribute the surrounding Mean-field Coulomb potential to the Hamiltonian.
		////*****************************************************************
		
		auto subIndexI = atomI.orbitalIndexList(optList[1]);
		for( unsigned ii = 0; ii<subIndexI.size() ; ii++){
			for( unsigned jj = 0; jj<subIndexI.size() ; jj++){
				auto & iterI = subIndexI[ii];
				auto & labelI = iterI.first;
				auto & indexI = iterI.second;
				
				auto & iterJ = subIndexI[jj];
				auto & labelJ = iterJ.first;
				auto & indexJ = iterJ.second;
				
				
				// Value "0.5" is to reduce the double counting.
				if( Lat.HSpace() == NORMAL ){
					// Normal part
					if( labelI == "u" and labelJ == "u"){
						hamElementList.push_back(MatrixElement(indexI, indexJ, cI_dd,	vec(0,0,0))); continue;
					}
					if( labelI == "d" and labelJ == "d"){
						hamElementList.push_back(MatrixElement(indexI, indexJ, cI_uu,	vec(0,0,0))); continue;
					}
					if( labelI == "u" and labelJ == "d"){
						hamElementList.push_back(MatrixElement(indexI, indexJ, cI_du,	vec(0,0,0))); continue;
					}
					if( labelI == "d" and labelJ == "u"){
						hamElementList.push_back(MatrixElement(indexI, indexJ, cI_ud,	vec(0,0,0))); continue;
					}
				}
				// Value "0.25" is to reduce the double counting and the particle/hole part.
				if( Lat.HSpace() == EXNAMBU){
					// Particle part
					if( labelI == "Au" and labelJ == "Au"){
						hamElementList.push_back(MatrixElement(indexI, indexJ, cI_dd,	vec(0,0,0))); continue;
					}
					if( labelI == "Ad" and labelJ == "Ad"){
						hamElementList.push_back(MatrixElement(indexI, indexJ, cI_uu,	vec(0,0,0))); continue;
					}
					if( labelI == "Au" and labelJ == "Ad"){
						hamElementList.push_back(MatrixElement(indexI, indexJ, cI_du,	vec(0,0,0))); continue;
					}
					if( labelI == "Ad" and labelJ == "Au"){
						hamElementList.push_back(MatrixElement(indexI, indexJ, cI_ud,	vec(0,0,0))); continue;
					}
					
					// Hole part
					if( labelI == "Bu" and labelJ == "Bu"){
						hamElementList.push_back(MatrixElement(indexI, indexJ,-cI_dd,	vec(0,0,0))); continue;
					}
					if( labelI == "Bd" and labelJ == "Bd"){
						hamElementList.push_back(MatrixElement(indexI, indexJ,-cI_uu,	vec(0,0,0))); continue;
					}
					if( labelI == "Bu" and labelJ == "Bd"){
						hamElementList.push_back(MatrixElement(indexI, indexJ,conj(cI_du),	vec(0,0,0))); continue;
					}
					if( labelI == "Bd" and labelJ == "Bu"){
						hamElementList.push_back(MatrixElement(indexI, indexJ,conj(cI_ud),	vec(0,0,0))); continue;
					}
				}
			}
		}
	}
	void		addIntraDudarevUJ(const PreprocessorInfo & preInfo)	{
		/*****************************************
		 * interU2J > Fe 1:2 >  U2J
		 *****************************************/
		auto & optList = preInfo.optList;
		auto & varList = preInfo.varList;
		
		if( Lat.parameter.STR("spin") == "off"){
			ErrorMessage(preInfo.filename, preInfo.lineNumber, " \""+preInfo.line+"\"\n"+
						 " Inter orbital U-J Dudarev model cannot be constructed under spin-independent space.");
		}
		if( Lat.HSpace() == NAMBU ){
			ErrorMessage(preInfo.filename, preInfo.lineNumber, " \""+preInfo.line+"\"\n"+
						 " a model for 3-dimensional spin is not applicable with NAMBU space.");
		}
		
		auto atomI = Lat.getAtom();
		if( atomI.atomName != optList[0]){ return; } // Name does not match.
		if( !atomI.hasOrbital(optList[1])){
			ErrorMessage(preInfo.filename, preInfo.lineNumber, " \""+preInfo.line+"\"\n"+
						 " The atom, \'"+atomI.atomName+"\', has no orbital, "+optList[1]+".");
		}

		// Find the first matched orbital .... interU2J > Fe '1':2 > U2J
		string orbitalNumber = atomI.getOrbitalNumber(optList[1]);
		string orderKey = "@:"+orbitalNumber+":4den";
		auto orderI = order.findOrder(atomI, orderKey);
		if( !orderI.first ){
			ErrorMessage(preInfo.filename, preInfo.lineNumber, " \""+preInfo.line+"\"\n"+
						 " The 4 density order parameter,"+orderKey+" , was not defined.");
		}
		
		r_var intraUJ = parseSiteString(atomI, varList[0])[0].real();
		for( unsigned i=1 ; i<varList.size() ; i++){
			intraUJ = intraUJ * parseSiteString(atomI, varList[i])[0].real();
		}
		
		intraUJ = 0.5 * intraUJ;
		
		r_var	N  = -intraUJ * orderI.second[0].real();
		r_var	Sx = -intraUJ * orderI.second[1].real();
		r_var	Sy = -intraUJ * orderI.second[2].real();
		r_var	Sz = -intraUJ * orderI.second[3].real();
		
		energyMap["4.DUJ Eng"] -= 0.5 * (N*N+Sx*Sx+Sy*Sy+Sz*Sz) / intraUJ;
		
		////*****************************************************************
		//// Attribute the on-site Mean-field Dudarev UJ potential to the Hamiltonian.
		////*****************************************************************
		switch( Lat.HSpace() ){
				
			case NORMAL: {
				unsigned indexNu = atomI.index(orbitalNumber+"u");
				unsigned indexNd = atomI.index(orbitalNumber+"d");
				
				hamElementList.push_back(MatrixElement(indexNu, indexNd, Sx-Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexNd, indexNu, Sx+Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexNu, indexNu, Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexNd, indexNd,-Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexNu, indexNu, intraUJ+N	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexNd, indexNd, intraUJ+N 	, vec(0,0,0)));
				break;
			}
			case NAMBU: {
				unsigned indexAu = atomI.index(orbitalNumber+"Au");
				unsigned indexBd = atomI.index(orbitalNumber+"Bd");
				hamElementList.push_back(MatrixElement(indexAu, indexAu, Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBd, Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexAu, indexAu, intraUJ+N	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBd,-intraUJ-N 	, vec(0,0,0)));
				if( Sx > 0.0001 or Sy > 0.0001){
					ErrorMessage(preInfo.filename,
								 preInfo.lineNumber,
								 " \""+preInfo.line+"\"\n"+
								 " Cannot be applied for a spin-depende Nambu space Hamiltonian.");
				}
				break;
			}
			case EXNAMBU: {
				unsigned indexAu = atomI.index(orbitalNumber+"Au");
				unsigned indexAd = atomI.index(orbitalNumber+"Ad");
				unsigned indexBu = atomI.index(orbitalNumber+"Bu");
				unsigned indexBd = atomI.index(orbitalNumber+"Bd");
				
				hamElementList.push_back(MatrixElement(indexAu, indexAd, Sx-Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexAd, indexAu, Sx+Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexAu, indexAu, Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexAd, indexAd,-Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexAu, indexAu, intraUJ+N	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexAd, indexAd, intraUJ+N 	, vec(0,0,0)));
				
				hamElementList.push_back(MatrixElement(indexBu, indexBd, Sx+Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBu, Sx-Im*Sy	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBu, indexBu,-Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBd, Sz			, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBu, indexBu,-intraUJ-N	, vec(0,0,0)));
				hamElementList.push_back(MatrixElement(indexBd, indexBd,-intraUJ-N 	, vec(0,0,0)));
				break;
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
		energyMap["3.U Eng"] = 0;
		energyMap["4.DUJ Eng"] = 0;
		energyMap["5.Pair Eng"] = 0;
		
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
		pairElementList.clear();
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
			for( auto & iter : tbm.hamPreprocessor.list_IntraHubbard)	{	addIntraHubbard	(iter);}
			for( auto & iter : tbm.hamPreprocessor.list_IntraDudarevUJ)	{	addIntraDudarevUJ(iter);}
			for( auto & iter : tbm.hamPreprocessor.list_QuantumFieldB)	{	addFieldB		(iter);}
			
			// If space==normal The following part will be ignored.
			if( Lat.HSpace() == NORMAL) continue;
			
			for( auto & iter : tbm.hamPreprocessor.list_PairingS)		{	addPairing		(iter,"S"); }
			for( auto & iter : tbm.hamPreprocessor.list_PairingUU)		{	addPairing		(iter,"UU"); }
			for( auto & iter : tbm.hamPreprocessor.list_PairingDD)		{	addPairing		(iter,"DD"); }
			for( auto & iter : tbm.hamPreprocessor.list_PairingUD)		{	addPairing		(iter,"UD"); }
			
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
	x_var			getTanhDensity			(r_mat bvec, unsigned index_i, unsigned index_j){
		/*
		 This term returns the sum over energy for the tanh leading pair wised wave function:
			\sum_n u^{n}_{index_i} v^{n*}_{index_j} tanh(E/2*Kb*T)
		 */
		
		auto Temperature = Lat.parameter.VAR("Temperature", 0.00001);
		if ( KEigenValVec.size()>0  and index_i < Lat.indexSize() and index_j < Lat.indexSize() ) {
			x_var tanh_den_ij=0;
			for (unsigned k=0; k<KEigenValVec.size(); k++) {
				auto & kp		= KEigenValVec[k].kPoint;
				auto & kEigVal	= KEigenValVec[k].eigenValue;
				auto & kEigVec	= KEigenValVec[k].eigenVector;
				
				x_var	phase = kp[0]*bvec[0]+ kp[1]*bvec[1]+ kp[2]*bvec[2];
				auto	exp_JI=exp(-Im*phase);
    
				for( int i=0 ; i<kEigVal.size() ; i++ ){
					double thEn = tanh(0.5*kEigVal[i]/Temperature.real());
					x_var uI = kEigVec(index_i,i);
					x_var vJ = kEigVec(index_j,i);
					tanh_den_ij += uI*conj(vJ)*thEn*exp_JI;
				}
			}
			return tanh_den_ij/KEigenValVec.size();
		}

		return 0;
	}
	
	void			calculate4DensityOrder	()						{
		if( Lat.parameter.STR("spin")!="on" ){
			ErrorMessage("Error, the spin space suld be turned \'on\' for the 4-density calculation");
		}
		
		double ratio = 1;
		if( Lat.HSpace()==EXNAMBU ) ratio = 0.5;
			
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
					if( parser_I[1] == "u" and parser_J[1] == "d" and Lat.HSpace()!=NAMBU){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[1] += tmpDen.real();
						r4den[2] += tmpDen.imag();
					}
					if( parser_I[1] == "d" and parser_J[1] == "u" and Lat.HSpace()!=NAMBU){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[1] += tmpDen.real();
						r4den[2] -= tmpDen.imag();
					}
					
					if( parser_I[1] == "Au" and parser_J[1] == "Au"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += tmpDen*ratio;
						r4den[3] += tmpDen*ratio;
					}
					if( parser_I[1] == "Ad" and parser_J[1] == "Ad"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += tmpDen*ratio;
						r4den[3] -= tmpDen*ratio;
					}
					if( parser_I[1] == "Au" and parser_J[1] == "Ad" and Lat.HSpace()!=NAMBU){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[1] += tmpDen.real()*ratio;
						r4den[2] += tmpDen.imag()*ratio;
					}                            
					if( parser_I[1] == "Ad" and parser_J[1] == "Au" and Lat.HSpace()!=NAMBU){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[1] += tmpDen.real()*ratio;
						r4den[2] -= tmpDen.imag()*ratio;
					}                            
					if( parser_I[1] == "Bu" and parser_J[1] == "Bu"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += (1-tmpDen)*ratio;
						r4den[3] += (1-tmpDen)*ratio;
					}
					if( parser_I[1] == "Bd" and parser_J[1] == "Bd"){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[0] += (1-tmpDen)*ratio;
						r4den[3] -= (1-tmpDen)*ratio;
					}
					if( parser_I[1] == "Bu" and parser_J[1] == "Bd" and Lat.HSpace()!=NAMBU){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[1] += -tmpDen.real()*ratio;
						r4den[2] += -tmpDen.imag()*ratio;
					}                             
					if( parser_I[1] == "Bd" and parser_J[1] == "Bu" and Lat.HSpace()!=NAMBU){
						x_var tmpDen = getDensityMatrix(index_I, index_J);
						r4den[1] += -tmpDen.real()*ratio;
						r4den[2] -= -tmpDen.imag()*ratio;
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
	void			calculatePairingOrder	()						{
		for( auto & elem : pairElementList ){
			
			double V = elem.V;
			x_var wave = 0;
			
			switch( Lat.HSpace() ){
			case NAMBU:
				if( elem.orderStr == "spair" ){
					wave += getTanhDensity( elem.atomPair.bondIJ(), elem.pIndexMap["jAu"], elem.pIndexMap["iBd"]);
					wave += getTanhDensity(-elem.atomPair.bondIJ(), elem.pIndexMap["iAu"], elem.pIndexMap["jBd"]);
					break;
				}
			case EXNAMBU:
				if( elem.orderStr == "spair" ){
					wave += 0.5*getTanhDensity( elem.atomPair.bondIJ(), elem.pIndexMap["jAd"], elem.pIndexMap["iBu"]);
					wave += 0.5*getTanhDensity(-elem.atomPair.bondIJ(), elem.pIndexMap["iAu"], elem.pIndexMap["jBd"]);
					wave += 0.5*getTanhDensity( elem.atomPair.bondIJ(), elem.pIndexMap["jAu"], elem.pIndexMap["iBd"]);
					wave += 0.5*getTanhDensity(-elem.atomPair.bondIJ(), elem.pIndexMap["iAd"], elem.pIndexMap["jBu"]);
					break;
				}
				if( elem.orderStr == "udpair" ){
					wave += 0.5*getTanhDensity( elem.atomPair.bondIJ(), elem.pIndexMap["jAd"], elem.pIndexMap["iBu"]);
					wave += 0.5*getTanhDensity(-elem.atomPair.bondIJ(), elem.pIndexMap["iAu"], elem.pIndexMap["jBd"]);
					wave -= 0.5*getTanhDensity( elem.atomPair.bondIJ(), elem.pIndexMap["jAu"], elem.pIndexMap["iBd"]);
					wave -= 0.5*getTanhDensity(-elem.atomPair.bondIJ(), elem.pIndexMap["iAd"], elem.pIndexMap["jBu"]);
					break;
				}
				if( elem.orderStr == "uupair" ){
					wave += getTanhDensity( elem.atomPair.bondIJ(), elem.pIndexMap["jAu"], elem.pIndexMap["iBu"]);
					wave -= getTanhDensity(-elem.atomPair.bondIJ(), elem.pIndexMap["iAu"], elem.pIndexMap["jBu"]);
					auto tmp = getTanhDensity(-elem.atomPair.bondIJ(), elem.pIndexMap["iAu"], elem.pIndexMap["jBu"]);
					break;
				}
				if( elem.orderStr == "ddpair" ){
					wave -= getTanhDensity( elem.atomPair.bondIJ(), elem.pIndexMap["jAd"], elem.pIndexMap["iBd"]);
					wave += getTanhDensity(-elem.atomPair.bondIJ(), elem.pIndexMap["iAd"], elem.pIndexMap["jBd"]);
					auto tmp = getTanhDensity(-elem.atomPair.bondIJ(), elem.pIndexMap["iAd"], elem.pIndexMap["jBd"]);
					break;
				}
			}
			
			x_mat pair(1,1);
			pair[0] = V*0.25*wave;
			order.setNew(elem.atomPair.atomI.atomIndex, elem.orderKey, pair);
		}
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



