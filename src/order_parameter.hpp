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
//  order_parameter.hpp
//  TBM^3
//
//  Created by Yuan Yen Tai on 7/19/16.
//

// Define the reverse compare map.
struct reverseStringCompare{
	bool operator()(const std::string& a, const std::string& b) const {
		string A = "", B = "";
		for( auto & ia : a){ A = ia+A;}
		for( auto & ib : b){ B = ib+B;}
		if( A.size() != B.size() ) return A.size() < B.size();
		return A < B;
	}
};
typedef map<string, x_mat, reverseStringCompare> rmap;


class OrderParameter{
private:
	Lattice & Lat;
	
	x_mat	empty_mat_for_return;
	vector<rmap>	optList;
	
	string	parseXVar(x_var val)			{
		if (val.imag() == 0) return DoubleToStr(val.real());
		return val.tostr();
	}

public:
	OrderParameter(Lattice & _Lat): Lat(_Lat){ }
	OrderParameter(const OrderParameter & order): Lat(order.Lat){
		*this = order;
	}
	~OrderParameter(){
		optList.clear();
	}
	
	void	init(){
		empty_mat_for_return = x_mat(1,1);
		
		optList.clear();
		for(unsigned ii=0 ; ii<Lat.latticeSize() ; ii++) optList.push_back(		rmap()		);
	}

	x_mat &				operator()(string opt)					{
		// xxx.order("Fe			cspin") = 1,2,3 (site only operation)
		// xxx.order("Fe:O:+1+0+0	1:2:pairingS") = xxx (bond operation)
		empty_mat_for_return = x_mat(1,1);
		
		/* -------------------------------------------
		 Detecting that if the opt formate is correct.
		 --------------------------------------------- */
		auto first = split(opt, " ");
		if( first.size() != 2){
			ErrorMessage("Error, order parameter operation:\n" +opt +"\n does not match the correct formate.");
		}
		removeSpace(first[1]);
		auto second1 = split(first[0], ":");
		auto second2 = split(first[1], ":");
		if( second1.size() != 1 and second1.size() != 3){
			ErrorMessage("Error, order parameter operation:\n" +opt +"\n does not match the correct formate.");
		}
		if( second2.size() != 1 and second2.size() != 2 and second2.size() != 3){
			ErrorMessage("Error, order parameter operation:\n" +opt +"\n does not match the correct formate.");
		}
		if( second1.size() == 3 and second2.size() != 3){
			ErrorMessage("Error, order parameter operation:\n" +opt +"\n does not match the correct formate.");
		}
		
		/* -------------------------------------------
		 Detecting that if the opt can be performed.
		 --------------------------------------------- */
		AtomPair atomPair;
		if( second1.size() == 1){
			auto atomI = Lat.getAtom();
			atomPair = AtomPair(atomI,atomI,vec(0,0,0));
			if( atomPair.atomI.atomName != second1[0])	{ return empty_mat_for_return; }
		}
		else if( second1.size() == 3){
			atomPair = Lat.getPair(second1[2]);
			if( atomPair.atomI.atomName != second1[0])	{ return empty_mat_for_return; }
			if( atomPair.atomJ.atomName != second1[1])	{ return empty_mat_for_return; }
			if(!atomPair.withinRange() )				{ return empty_mat_for_return; }
		}
		
		if( second2.size() == 2 ){
			if( !atomPair.atomI.hasOrbital(second2[0]))	{ return empty_mat_for_return; }
			if( !atomPair.atomJ.hasOrbital(second2[0]))	{ return empty_mat_for_return; }
		}
		if( second2.size() == 3 ){
			if( !atomPair.atomI.hasOrbital(second2[0]))	{ return empty_mat_for_return; }
			if( !atomPair.atomJ.hasOrbital(second2[1]))	{ return empty_mat_for_return; }
		}
		
		/* -------------------------------------------
		 Formate the optKey and make a storage for the value.
		 --------------------------------------------- */
		string optKey = "";
		if		( second1.size() == 1){ optKey += "@"; }
		else if	( second1.size() == 3){ optKey += vecToStr(atomPair.bondIJ()); }
		
		if( second2.size() == 1){
			optKey += ":"+second2[0];
		}	
		if( second2.size() == 2){
			optKey += ":"+atomPair.atomI.getOrbitalNumber(second2[0]);
			optKey += ":"+second2[1];
		}
		if( second2.size() == 3){
			optKey += ":"+atomPair.atomI.getOrbitalNumber(second2[0]);
			optKey += ":"+atomPair.atomJ.getOrbitalNumber(second2[1]);
			optKey += ":"+second2[2];
		}
		
		// Storage for the value.
		auto  it = optList[atomPair.atomI.atomIndex].find(optKey);
		if( it != optList[atomPair.atomI.atomIndex].end() ) {
			return optList[atomPair.atomI.atomIndex][optKey];
		}
		else{
			optList[atomPair.atomI.atomIndex][optKey] = empty_mat_for_return;
			return optList[atomPair.atomI.atomIndex][optKey];
		}
		
		return empty_mat_for_return;
	}
	
	pair<bool, x_mat>	findOrder(string atomName, string opt)	{
		auto atomI = Lat.getAtom();
		bool has_order =
			(atomI.atomName == atomName) and
			(optList[atomI.atomIndex].find(opt) != optList[atomI.atomIndex].end());
		
		x_mat xvec(1,1);
		if( has_order ){
			return make_pair(has_order, optList[atomI.atomIndex][opt]);
		}
		
		return make_pair(has_order, xvec);
	}
	pair<bool, x_mat>	findOrder(Atom atomI, string opt)		{
		bool has_order = optList[atomI.atomIndex].find(opt) != optList[atomI.atomIndex].end();
		
		x_mat xvec(1,1);
		if( has_order ){
			return make_pair(has_order, optList[atomI.atomIndex][opt]);
		}
		
		return make_pair(has_order, xvec);
	}
	
	void	set(unsigned atomIndex, string opt, x_mat & newVal)	{
		if( optList[atomIndex].find(opt) != optList[atomIndex].end() ){
			optList[atomIndex][opt] = newVal;
		}
	}
	void	setNew(unsigned atomIndex, string opt, x_mat & newVal)	{
		optList[atomIndex][opt] = newVal;
	}
	
	void	load(string sub_filename = "")						{
		init();
		
		string filename = Lat.FileName() + ".ord";
		if( sub_filename != "") filename += "." + sub_filename;
		
		
        ifstream infile(filename);
		
		if ( infile.is_open() ) {
            string line;
			int atomIndex = -1;
			
            while ( getline(infile, line) ) {
				
				istringstream iss(line);
				string header = "";
				
				iss >> header;
				if( header == ">>>"){
					iss >> atomIndex ;
					continue;
				}
				else if( atomIndex >= 0){
					auto lineKeyVal = split(line, "=");
					
					string & optKey = lineKeyVal[0];
					string & optVal = lineKeyVal[1];
					removeSpace(optKey);
					
					auto keyParser = split(optKey, ":");
					AtomPair atomPair;
					if( keyParser[0] == "@")	atomPair = AtomPair(Lat.getAtom(atomIndex), Lat.getAtom(atomIndex), vec(0,0,0));
					else						atomPair = Lat.getPair(atomIndex, keyParser[0]);
					
					optList[atomIndex][optKey] = StrToXVec(optVal);
				}
			}
		}

		infile.close();
	}
	void	load(vector<pair<string, x_mat> > orderOptList)		{
		init();
		
		if ( orderOptList.size() > 0 ) {
			// Cannot find the file and start construct the xxx.lat.ord.
			while( Lat.iterate() ) {
				for(auto & iter: orderOptList){
					// Construct the order parameter
					operator()(iter.first) = iter.second;
				}
			}
		}
	}
	void	save(string sub_filename = "")						{
		
		string filename = Lat.FileName() + ".ord";
		if( sub_filename != "") filename += "." + sub_filename;
		
		ofstream outfile(filename);
		for(unsigned ii=0 ; ii<optList.size() ; ii++){
			auto & atOpt = optList[ii];
			auto atomI = Lat.getAtom(ii);
			outfile	<< fformat(">>> "+IntToStr(ii), 7)<<" "
					<<fformat(atomI.orbitalOriginIndex,4)<<" "
					<<fformat(atomI.atomName,4)<<" "<<atomI.pos<<endl;
			
			for( auto & opt: atOpt){
				outfile<< " "+opt.first <<" = ";
				for(unsigned i=0 ; i<opt.second.size() ; i++){
					outfile << parseXVar(opt.second[i]);
					if(i<opt.second.size() -1) outfile << ",";
				}
				outfile<<endl;
			}
		}
		outfile.close();
	}
	void	clear()												{
		optList.clear();
	}

	vector<rmap> &	getOptList()								{ return optList; }
	void			setOptList(vector<rmap> & optL)				{
		optList = optL;
	}
	
	/*
	 Here is the definition of the i/o file formate:
	 
	 >>> 0	1	Fe		[[ 0	0	0 ]]
	  @:cspin				=	1,2,3
	  @:cspin				=	1,2,3
	  @:1:pspin				=	1,2,3
	  @:2:pspin				=	1,2,3
	  @:1:den				=	1
	  @:2:den				=	1
	  +1+0+0:1:2:pairingS	=	(1,2)
	  +0+1+0:1:2:pairingS	=	(1,2)
	  -1+0+0:1:2:pairingS	=	(1,2)
	  +0-1+0:1:2:pairingS	=	(1,2)
	  +1+0+0:2:1:pairingS	=	(1,2)
	  +0+1+0:2:1:pairingS	=	(1,2)
	  -1+0+0:2:1:pairingS	=	(1,2)
	  +0-1+0:2:1:pairingS	=	(1,2)
	 */
	
	OrderParameter & operator=( const OrderParameter & order)	{
		Lat = order.Lat;
		optList = order.optList;
		return *this;
	}
	
	friend OrderParameter  operator+(OrderParameter  L_val, OrderParameter  R_val)	{
		OrderParameter	ret_order(L_val);
		if( L_val.optList.size() != R_val.optList.size() ){
			ErrorMessage("Error, cannot add two different order parameters.");
		}
		
		for( unsigned i=0 ; i<ret_order.optList.size() ; i++){
			
			auto & retOpt	= ret_order.optList[i];
			auto & LOpt		= L_val.optList[i];
			auto & ROpt		= R_val.optList[i];
			
			for( auto & elem: retOpt ){
				auto & key = elem.first;
				retOpt[key] = LOpt[key] + ROpt[key];
			}
		}
		return ret_order;
	}
	
	friend OrderParameter  operator-(OrderParameter  L_val, OrderParameter  R_val)	{
		OrderParameter	ret_order(L_val);
		if( L_val.optList.size() != R_val.optList.size() ){
			ErrorMessage("Error, cannot subtract two different order parameters.");
		}
		
		for( unsigned i=0 ; i<ret_order.optList.size() ; i++){
			
			auto & retOpt	= ret_order.optList[i];
			auto & LOpt		= L_val.optList[i];
			auto & ROpt		= R_val.optList[i];
			
			for( auto & elem: retOpt ){
				auto & key = elem.first;
				retOpt[key] = LOpt[key] - ROpt[key];
			}
		}
		return ret_order;
	}
	
	friend OrderParameter  operator*(OrderParameter  val, double number)			{
		OrderParameter	ret_order(val);
		
		for( unsigned i=0 ; i<ret_order.optList.size() ; i++){
			
			auto & retOpt	= ret_order.optList[i];
			
			for( auto & elem: retOpt ){
				auto & key = elem.first;
				retOpt[key] = retOpt[key] * number;
			}
		}
		return ret_order;
	}
	
	friend OrderParameter  operator*(double number, OrderParameter  val)			{
		OrderParameter	ret_order(val);
		
		for( unsigned i=0 ; i<ret_order.optList.size() ; i++){
			
			auto & retOpt	= ret_order.optList[i];
			
			for( auto & elem: retOpt ){
				auto & key = elem.first;
				retOpt[key] = retOpt[key] * number;
			}
		}
		return ret_order;
	}
};

























