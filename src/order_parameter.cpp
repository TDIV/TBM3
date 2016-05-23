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
//  order_parameter.cpp
//  TBM^3
//
//  Created by Yuan Yen Tai on 8/04/15.
//


struct SiteOperationUnit{
	Atom		SiteI;
	string		opt;
	x_var		val;
	SiteOperationUnit(	Atom _SiteI,	string _opt,	x_var _val):
					SiteI(_SiteI),	opt(_opt),		val(_val)
	{ }
};

struct PairOperationUnit{
	AtomPair	pit;
	string		opt;
	x_var		val;
	PairOperationUnit(	AtomPair ap,	string _opt,	x_var _val):
						pit(ap),	opt(_opt),		val(_val)
	{ }
};

/* -------------------------------------------------------------------
 The OrderParameter class is designed to create arbitrary "real-space" orders according
 to the input lattice structure.
 It can be used to construct the Hamiltonian when calling the super class in "tbm.cpp".
 -------------------------------------------------------------------*/
class OrderParameter{
protected:
	x_var					empty_value_for_return;
	
	// Data structure for the site
	map<string, unsigned>	siteMap;
	vector<SiteOperationUnit> siteOpterationList;
	
	// Data structure for the pair
	map<string, unsigned>	pairMap;
	vector<PairOperationUnit> pairOpterationList;
	
	
	// The lattice structure will be used to guide the orderParameter
	Lattice Lat;
public:
	unsigned	siteMapSize		()				const	{ return siteMap.size(); }
	unsigned	siteListSize	()				const	{ return siteOpterationList.size(); }
	
	unsigned	pairMapSize		()				const	{ return pairMap.size(); }
	unsigned	pairListSize	()				const	{ return pairOpterationList.size(); }
	
	string		sString			(unsigned ii)	const	{ return siteOpterationList[ii].opt; }
	x_var		sValue			(unsigned ii)	const	{ return siteOpterationList[ii].val; }
	x_var	&	setSiteValue	(unsigned ii)			{ return siteOpterationList[ii].val; }
	SiteOperationUnit	getSiteOperationAt(unsigned ii) {return siteOpterationList[ii];}
	
	string		pString			(unsigned ii)	const	{ return pairOpterationList[ii].opt; }
	x_var		pValue			(unsigned ii)	const	{ return pairOpterationList[ii].val; }
	x_var	&	setPairValue	(unsigned ii)			{ return pairOpterationList[ii].val; }
	PairOperationUnit	getPairOperationAt(unsigned ii) {return pairOpterationList[ii];}
	
public:
	OrderParameter	(Lattice _lat)					{
		Lat = _lat;
		empty_value_for_return = 0;
	}
	OrderParameter	( const OrderParameter & order)	{
		*this = order;
	}
	~OrderParameter	()								{clear();}

	// --------------------------------------------------------
	// - operator()-Access to the order parameters ------------
	// --------------------------------------------------------
	// Example:
	// OrderParameter order;
	// order( At, "Fe 1u"		)=val;
	// order( At, "Fe SR"		)=val;
	x_var &			operator()	(const Atom & At, string opt){
		empty_value_for_return = 0;
		
		auto word = split(opt, " ");
		
		if (word.size()!=2) {
			cout<<"Warning, the option formate is wrong: "+opt;
			return empty_value_for_return;
		}
		
		if (word[0] != At.Name()) {
			return empty_value_for_return;
		}
		
		string skey = IntToStr(At.atomIndex()) +" "+ word[1];
		
		auto it = siteMap.find(skey);
		if (it == siteMap.end()) {
			siteMap[skey] = siteOpterationList.size();
			siteOpterationList.push_back(SiteOperationUnit(At, skey, 0));
			return siteOpterationList[siteMap[skey]].val;
		} else {
			return siteOpterationList[siteMap[skey]].val;
		}
	}
	x_var &			order		(const Atom & At, string opt){ return operator()(At, opt); }
	vector<x_var>	getVars		(const Atom & At, string opt){
		auto word = split(opt, " ");
		vector<x_var> ret_vars;
		if (word.size() >= 1) {
			for (unsigned i=0; i<word.size(); i++) {
				ret_vars.push_back(order(At, At.Name()+" "+word[i]));
			}
		}
		return ret_vars;
	}
	
	// Example:    v-----str-key----v   v---value
	// order( Ap, "Fe:Fe:+1+0+0 pair.1u1d")=val;
	x_var &			operator()	(const AtomPair & Ap, string opt){
		empty_value_for_return = 0;

		auto word = split(opt, " ");

		if (word.size()!=2) {
			cout<<"Warning, the option formate is wrong: "+opt;
			return empty_value_for_return;
		}

		auto ss = split(word[0], ":");
		if (ss.size()!=3) {
			cout<<"Warning, the option formate is wrong: "+opt;
			return empty_value_for_return;
		}
		
		if ( Lat.translateBondString( ss[2] ) != Ap.bond ) {
			return empty_value_for_return;
		}
		
		if ( !(ss[0] == Ap.AtomI.Name() and ss[1] == Ap.AtomJ.Name()) ){
			return empty_value_for_return;
		}
		
		// ---------------------------------------------------------------------
		// All conditions passed.
		// Now, assign/update new pair-order into the pool.
		string skey = IntToStr(Ap.AtomI.atomIndex())+":"+IntToStr(Ap.AtomJ.atomIndex())+":"+Ap.bond+" "+ word[1];
		
		auto it = pairMap.find(skey);
		if (it == pairMap.end()) {
			pairMap[skey] = pairOpterationList.size();
			
			pairOpterationList.push_back(PairOperationUnit(Ap, skey, 0));
			//pairString.push_back(skey);
			//var.push_back(0);
			return pairOpterationList[pairMap[skey]].val;
		} else {
			return pairOpterationList[pairMap[skey]].val;
		}
		// ---------------------------------------------------------------------
		
		return empty_value_for_return;
	}
	x_var &			order		(const AtomPair & Ap, string opt){ return operator()(Ap, opt);}
	
	// --------------------------------------------------------
	
	string		getSiteOrderStr	(Lattice Lat)			{
		/* Firstly -----------------------------------
		 Save the On-site order parameter in the formate:
		 
		 ...
		 #S On-site order parameters
		 
		 |----> Unitcell
		 |    |----> Atom index
		 |    |    |----> Atom name and sub-atom index
		 |    |    |        |------|------|-----> Atom position
		 v    v    v        v      v      v          v-----v-----> order name and value
		 
		 0    0    Bi-0 [[  0.5    0.5    0.5  ]] >>
		 0    1    Fe-1 [[  0      0      0    ]] >> theta:(3.1415927,0)       phi:(0,0)
		 0    2    O-2  [[  0.5    0      0    ]] >>
		 0    3    O-3  [[  0      0.5    0    ]] >>
		 0    4    O-4  [[  0      0      0.5  ]] >>
		 1    5    Bi-0 [[  1.5    0.5    0.5  ]] >>
		 1    6    Fe-1 [[  1      0      0    ]] >> theta:(0,0)               phi:(0,0)
		 1    7    O-2  [[  1.5    0      0    ]] >>
		 1    8    O-3  [[  1      0.5    0    ]] >>
		 1    9    O-4  [[  1      0      0.5  ]] >>
		 2    10   Bi-0 [[  0.5    1.5    0.5  ]] >>
		 */
		auto siter = Lat.site_iteration();
		
		ostringstream oss;
		
		if (siteOpterationList.size()>0) {
			oss<<"#S On-site order parameters."<<endl;
			for (unsigned ii=0 ; ii<siter.size(); ii++) {
				auto si = siter[ii];
				auto label = si.indexLabel;
				r_mat ppos(1,3);
				ppos.setPrintLength(10);
				for (unsigned i=0; i<ppos.size(); i++) { ppos[i]=si.pos[i]; }
				
				// Output the Atom index, name and position.
				oss<< fmt(si.unitcellIndex(),4)<<" "<< fmt(si.atomIndex(),4)<<" "<< fmt(si.subName(),4)<<" "<< ppos<<" >>      " ;
				
				for (unsigned i=0; i<siteOpterationList.size(); i++) {
					auto	opt		= siteOpterationList[i].opt;
					auto	opt_word= split(opt, " ");
					
					string at_name=siter[StrToInt(opt_word[0])].Name();
					string output=order(si, at_name+" "+opt_word[1]).tostr();
					
					if (StrToInt(opt_word[0]) == ii) {
						oss<<opt_word[1]<<":"<<fmt(output, 34)<<" ";
					}
				}
				oss<<endl;
			}
		}

		return oss.str();
	}
	string		getPairOrderStr	(Lattice Lat)			{
		/*
		 Save the bond-pair order parameter in the formate:
		 ...
		 #P atom-pair order parameters
		 |----> Unitcell
		 |    |----> Atom index
		 |    |    |----> Atom name and sub-atom index
		 |    |    |        |------|------|-----> Atom position
		 |    |    |        |      |      |			 |--|----|------|---> bond-direction, order-name, orbital/spin specific(optional) and order-value
		 v    v    v        v      v      v          v  v    v      v
		 
		 0    0    Bi-0 [[  0.5    0.5    0.5  ]] >>
		 0    1    Fe-1 [[  0      0      0    ]] >> +x:pair.1u.1d:(0.1,0)       -x:pair.1u.1d:(0.1,0)
		 0    2    O-2  [[  0.5    0      0    ]] >>
		 0    3    O-3  [[  0      0.5    0    ]] >>
		 0    4    O-4  [[  0      0      0.5  ]] >>
		 1    5    Bi-0 [[  1.5    0.5    0.5  ]] >>
		 1    6    Fe-1 [[  1      0      0    ]] >> +x:pair.1u.1d:(0.1,0)       -x:pair.1u.1d:(0.1,0)
		 1    7    O-2  [[  1.5    0      0    ]] >>
		 1    8    O-3  [[  1      0.5    0    ]] >>
		 1    9    O-4  [[  1      0      0.5  ]] >>
		 2    10   Bi-0 [[  0.5    1.5    0.5  ]] >>
		 */
		//auto siter = Lat.site_iteration();
		
		
		map<int, vector<string> > pairOrderMap;
		
		// Construct the pairOrderMap.
		for (unsigned i=0; i<pairOpterationList.size(); i++) {
			
			/* Begin checking the structure of the pairString */
			auto sp_base0 = split(pairOpterationList[i].opt, " ");
			
			if (sp_base0.size() != 2) { cout<<"Warning, not an effective formate: "<<pairOpterationList[i].opt<<". Ignored!"<<endl; break; }
			
			auto sp_base1 = split(sp_base0[0], ":");
			
			if (sp_base1.size() != 3) { cout<<"Warning, not an effective formate: "<<pairOpterationList[i].opt<<". Ignored!"<<endl; break; }
			/* End checking the structure of the pairString */
			
			auto indexI = StrToInt(sp_base1[0]);
			auto bondStr = Lat.translateBondString( sp_base1[2] );
			auto orderName= sp_base0[1];
			auto orderStr = pairOpterationList[i].val.tostr();
			
			auto pairOrderStr = bondStr+":"+orderName+":"+orderStr;
			
			auto it = pairOrderMap.find(indexI);
			if (it == pairOrderMap.end()) {
				vector<string> tmpStr;
				tmpStr.push_back(pairOrderStr);
				pairOrderMap[indexI] = tmpStr;
			} else {
				pairOrderMap[indexI].push_back(pairOrderStr);
			}
		}
		
		auto siter = Lat.site_iteration();
		ostringstream oss;
		
		if (pairOpterationList.size()>0) {
			oss<<"#P bond-pair order parameters."<<endl;
			for (unsigned ii=0 ; ii<siter.size(); ii++) {
				auto si = siter[ii];
				auto label = si.indexLabel;
				r_mat ppos(1,3);
				ppos.setPrintLength(10);
				for (unsigned i=0; i<ppos.size(); i++) { ppos[i]=si.pos[i]; }
				
				// Output the Atom index, name and position.
				oss<< fmt(si.unitcellIndex(),4)<<" "<< fmt(si.atomIndex(),4)<<" "<< fmt(si.subName(),4)<<" "<< ppos<<" >>      " ;
				
				for (unsigned i=0; i<pairOrderMap[si.atomIndex()].size(); i++) {
					auto	opt		= pairOrderMap[si.atomIndex()][i];
					
					oss<<fmt(opt,34)<<" ";
				}
				oss<<endl;
			}
		}
		
		return oss.str();
	}
	
	x_var		totalSiteOrder	()						{
		x_var total=0;
		for (unsigned i=0; i<siteOpterationList.size(); i++) {
			total += siteOpterationList[i].val;
		}
		return total;
	}
	x_var		totalPairOrder	()						{
		x_var total=0;
		for (unsigned i=0; i<pairOpterationList.size(); i++) {
			total+=pairOpterationList[i].val;
		}
		return total;
	}

	OrderParameter & operator=	(const OrderParameter & order){
		empty_value_for_return = 0;
		siteMap		= order.siteMap;
		siteOpterationList = order.siteOpterationList;
		
		pairMap		= order.pairMap;
		pairOpterationList = order.pairOpterationList;
		
		Lat			= order.Lat;
		return *this;
	}
	
	/* Clear all the stored information and data structures */
	void		clear			()						{
		siteMap.clear();
		siteOpterationList.clear();
		
		pairMap.clear();
		pairOpterationList.clear();
	}
	
	/* Preserve the datastructure and set all values to zero*/
	void		zerolize		()						{
		for (unsigned i=0; i<pairOpterationList.size(); i++) {
			pairOpterationList[i].val = 0;
		}
	}
	
	double		max_abs			()						{
		double mabs = 0;
		for (unsigned i=0; i<siteOpterationList.size(); i++) {
			double val = abs(siteOpterationList[i].val);
			if (mabs<val) mabs=val;
		}
		return mabs;
	}
	
	// ---- Save/load the order parameter in a unified formate.
	bool		save(string order_sub_name)							{

		auto siter = Lat.site_iteration();
		auto piter = Lat.pair_iteration();
		
		if (order_sub_name != "") {
			order_sub_name  = "."+order_sub_name  ;
		}
		
		ofstream out(Lat.filename+".ord"+order_sub_name  );
		
		// Storage of the site-order parameter
		if (siteOpterationList.size()>0) {
			if (out.is_open()) {
				out<<getSiteOrderStr(Lat);
			}
		}

		// Storage of the pair-order parameter
		if (pairOpterationList.size()>0) {
			if (out.is_open()) {
				out<<getPairOrderStr(Lat);
			}
		}
		
		out.close();
		return true;
	}
	bool		load(string order_sub_name, unsigned name_flag = 2)	{
		
		auto siter = Lat.site_iteration();
		auto piter = Lat.pair_iteration();
		
		if (order_sub_name != "") { order_sub_name = "."+order_sub_name; }
		
		string		filename_1 = Lat.filename+order_sub_name;
		string		filename_2 = Lat.filename+".ord"+order_sub_name;
		string		filename = filename_2;
		
		switch (name_flag) {
			case 1 : filename = filename_1; break;
			case 2 : filename = filename_2; break;
		}
	
		
		ifstream infile(filename);
		// ------ Read in the file and construct the parts
        if ( infile.is_open() ) {
			siteMap.clear();
			siteOpterationList.clear();
			pairMap.clear();
			pairOpterationList.clear();
	
            string line;
			string flag="";
			
            while ( getline(infile, line) ) {
				
				auto line_words = split(line, " ");
				
				// --------------------------------------------------------
				// Read in the site order parameter
				if (line_words[0] == "#S") { flag = "#S"; continue;}
				if ( flag == "#S" and line_words.size() > 9 ) {
					
					auto atomI = Lat(StrToInt(line_words[0]));
					
					for	(unsigned i=9; i<line_words.size(); i++) {
						auto str_variable = split(line_words[i], ":");
						auto var_name = str_variable[0];
						auto var = StrToComplex(str_variable[1]);
						
						auto sstr = line_words[1]+" "+var_name;
						siteOpterationList.push_back(SiteOperationUnit(atomI, sstr, var));
						siteMap[sstr] = siteOpterationList.size()-1;
					}
				}
				
				// --------------------------------------------------------
				// Read in the pair order parameter
				if (line_words[0] == "#P") { flag = "#P"; continue;}
				if ( flag == "#P" and line_words.size() > 9 ) {
					// The key-value formate of the bond order
					// v-----str-key----v   v---value
					//"1:2:-x pair.1u1d"    val;
					//string skey = IntToStr(Ap.AtomI.AtomIndex())+":"+IntToStr(Ap.AtomJ.AtomIndex())+":"+Ap.bond+" "+ word[1];
					
					
					
					auto atomI = Lat(StrToInt(line_words[0]));
					
					for	(unsigned i=9; i<line_words.size(); i++) {
						auto str_variable = split(line_words[i], ":");
						
						// Check the structure of line_words[i].
						if (str_variable.size() !=3) {
							cout<<"Warning, not a valid bond-pair order formate: "<<line_words[i]<<" ... ignored!"<<endl;
							continue;
						}
						
						string var_bond = str_variable[0];			// bond variable ex: "-x"
						auto var_name = str_variable[1];			// name of the order ex: "hop"
						auto pairVar = StrToComplex(str_variable[2]);	// value of the order ex: "(1.02,-2.001)"
						
						auto atomJ = Lat(atomI, var_bond);
						
						auto strIndexI = IntToStr(atomI.atomIndex());
						auto strIndexJ = IntToStr(atomJ.atomIndex());
						
						AtomPair atomPair;
						atomPair.AtomI = atomI;
						atomPair.AtomJ = atomJ;
						atomPair.bond = var_bond;
						atomPair.bvec = Lat.BVEC(var_bond);
						
						auto pairKeyStr = strIndexI+":"+strIndexJ+":"+var_bond+" "+var_name;
						
						pairOpterationList.push_back(PairOperationUnit(atomPair, pairKeyStr, pairVar));
						pairMap[pairKeyStr] = pairOpterationList.size()-1;
					}
				}
			}
		}
		return true;
	}
	
};

OrderParameter operator+(const OrderParameter & L_val, const OrderParameter & R_val) {
	OrderParameter ret_val(L_val);
	
	if (
			L_val.siteMapSize()		== R_val.siteMapSize()		and
			L_val.siteListSize()	== R_val.siteListSize()	and
			L_val.pairMapSize()		== R_val.pairMapSize()		and
			L_val.pairListSize()	== R_val.pairListSize()
		) {
		
		for (unsigned i=0; i<ret_val.siteListSize(); i++) {
			ret_val.setSiteValue(i) = L_val.sValue(i)+R_val.sValue(i);
		}
		
		for (unsigned i=0; i<ret_val.pairListSize(); i++) {
			ret_val.setPairValue(i) = L_val.pValue(i)+R_val.pValue(i);
		}
	}
	return ret_val;
}

OrderParameter operator-(const OrderParameter & L_val, const OrderParameter & R_val) {
	OrderParameter ret_val(L_val);
	
	if (
			L_val.siteMapSize()	== R_val.siteMapSize() and
			L_val.siteListSize() == R_val.siteListSize() and
			L_val.pairMapSize()	== R_val.pairMapSize() and
			L_val.pairListSize()	== R_val.pairListSize()
		) {
		
		for (unsigned i=0; i<ret_val.siteListSize(); i++) {
			ret_val.setSiteValue(i) = L_val.sValue(i)-R_val.sValue(i);
		}
		
		for (unsigned i=0; i<ret_val.pairListSize(); i++) {
			ret_val.setPairValue(i) = L_val.pValue(i)-R_val.pValue(i);
		}
	}
	return ret_val;
}

OrderParameter operator*(const OrderParameter & O_val, double d_val) {
	OrderParameter ret_val(O_val);
	
	for (unsigned i=0; i<ret_val.siteListSize(); i++) { ret_val.setSiteValue(i) = ret_val.sValue(i)*d_val; }
	for (unsigned i=0; i<ret_val.pairListSize(); i++) { ret_val.setPairValue(i) = ret_val.pValue(i)*d_val; }
	
	return ret_val;
}

OrderParameter operator*(double d_val, const OrderParameter & O_val) {
	OrderParameter ret_val(O_val);
	
	for (unsigned i=0; i<ret_val.siteListSize(); i++) { ret_val.setSiteValue(i) = ret_val.sValue(i)*d_val; }
	for (unsigned i=0; i<ret_val.pairListSize(); i++) { ret_val.setPairValue(i) = ret_val.pValue(i)*d_val; }
	
	return ret_val;
}














