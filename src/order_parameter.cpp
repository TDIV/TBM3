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

/* -------------------------------------------------------------------
 The OrderParameter class is designed to create arbitrary "real-space" orders according
 to the input lattice structure.
 It can be used to construct the Hamiltonian when calling the super class in "tbm.cpp".
 -------------------------------------------------------------------*/

class OrderParameter{
private:
	x_var					empty_value_for_return;
	
	// Data structure for the site
	map<string, unsigned>	siteMap;
	vector<string>			siteString;
	vector<x_var>			siteValue;
	
	// Data structure for the pair
	map<string, unsigned>	pairMap;
	vector<string>			pairString;
	vector<x_var>			pairValue;
	
	// The lattice structure will be used to guide the orderParameter
	Lattice Lat;
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
		
		string skey = IntToStr(At.AtomIndex()) +" "+ word[1];
		
		auto it = siteMap.find(skey);
		if (it == siteMap.end()) {
			siteMap[skey] = siteString.size();
			siteString.push_back(skey);
			siteValue.push_back(0);
			return siteValue[siteMap[skey]];
		} else {
			return siteValue[siteMap[skey]];
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
	// order( Ap, "Fe:Fe:-x pair.1u1d")=val;
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
		// All the condition passed, assign/update new pair-order into the pool.
		string skey = IntToStr(Ap.AtomI.AtomIndex())+":"+IntToStr(Ap.AtomJ.AtomIndex())+":"+Ap.bond+" "+ word[1];
		
		auto it = pairMap.find(skey);
		if (it == pairMap.end()) {
			pairMap[skey] = pairString.size();
			pairString.push_back(skey);
			pairValue.push_back(0);
			return pairValue[pairMap[skey]];
		} else {
			return pairValue[pairMap[skey]];
		}
		// ---------------------------------------------------------------------
		
		return empty_value_for_return;
	}
	x_var &			order		(const AtomPair & Ap, string opt){ return operator()(Ap, opt);}
	// --------------------------------------------------------
	
	string		getSiteOrderStr	(Lattice Lat)	{
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
		
		if (siteString.size()>0) {
			oss<<"#S On-site order parameters."<<endl;
			for (unsigned ii=0 ; ii<siter.size(); ii++) {
				auto si = siter[ii];
				auto label = si.index_label;
				r_mat ppos(1,3);
				ppos.setPrintLength(10);
				for (unsigned i=0; i<ppos.size(); i++) { ppos[i]=si.pos[i]; }
				
				// Output the Atom index, name and position.
				oss<< fmt(si.unitcell_index(),4)<<" "<< fmt(si.AtomIndex(),4)<<" "<< fmt(si.SubName(),4)<<" "<< ppos<<" >>      " ;
				
				for (unsigned i=0; i<siteString.size(); i++) {
					auto	opt		= siteString[i];
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
	string		getPairOrderStr	(Lattice Lat)	{
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
		for (unsigned i=0; i<pairString.size(); i++) {
			
			/* Begin checking the structure of the pairString */
			auto sp_base0 = split(pairString[i], " ");
			
			if (sp_base0.size() != 2) { cout<<"Warning, not an effective formate: "<<pairString[i]<<". Ignored!"<<endl; break; }
			
			auto sp_base1 = split(sp_base0[0], ":");
			
			if (sp_base1.size() != 3) { cout<<"Warning, not an effective formate: "<<pairString[i]<<". Ignored!"<<endl; break; }
			/* End checking the structure of the pairString */
			
			auto indexI = StrToInt(sp_base1[0]);
			auto bondStr = Lat.translateBondString( sp_base1[2] );
			auto orderName= sp_base0[1];
			auto orderStr = pairValue[i].tostr();
			
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
		
		if (pairString.size()>0) {
			oss<<"#P bond-pair order parameters."<<endl;
			for (unsigned ii=0 ; ii<siter.size(); ii++) {
				auto si = siter[ii];
				auto label = si.index_label;
				r_mat ppos(1,3);
				ppos.setPrintLength(10);
				for (unsigned i=0; i<ppos.size(); i++) { ppos[i]=si.pos[i]; }
				
				// Output the Atom index, name and position.
				oss<< fmt(si.unitcell_index(),4)<<" "<< fmt(si.AtomIndex(),4)<<" "<< fmt(si.SubName(),4)<<" "<< ppos<<" >>      " ;
				
				for (unsigned i=0; i<pairOrderMap[si.AtomIndex()].size(); i++) {
					auto	opt		= pairOrderMap[si.AtomIndex()][i];
					
					oss<<fmt(opt,34)<<" ";
				}
				oss<<endl;
			}
		}
		
		return oss.str();
	}
	
	x_var		totalSiteOrder	(){
		x_var total=0;
		for (unsigned i=0; i<siteValue.size(); i++) {
			total+=siteValue[i];
		}
		return total;
	}
	x_var		totalPairOrder	(){
		x_var total=0;
		for (unsigned i=0; i<pairString.size(); i++) {
			total+=pairValue[i];
		}
		return total;
	}

	OrderParameter & operator=	(const OrderParameter & order){
		empty_value_for_return = 0;
		siteMap		= order.siteMap;
		siteString	= order.siteString;
		siteValue	= order.siteValue;
		
		pairMap		= order.pairMap;
		pairString	= order.pairString;
		pairValue	= order.pairValue;
		
		Lat			= order.Lat;
		return *this;
	}
	
	unsigned	siteMapSize		()				const	{ return siteMap.size(); }
	unsigned	siteStringSize	()				const	{ return siteString.size(); }
	unsigned	siteValueSize	()				const	{ return siteValue.size(); }
	
	unsigned	pairMapSize		()				const	{ return pairMap.size(); }
	unsigned	pairStringSize	()				const	{ return pairString.size(); }
	unsigned	pairValueSize	()				const	{ return pairValue.size(); }
	
	string		sString			(unsigned ii)	const	{ return siteString[ii]; }
	x_var		sValue			(unsigned ii)	const	{ return siteValue[ii]; }
	x_var	&	setSiteValue	(unsigned ii)			{ return siteValue[ii]; }
	
	string		pString			(unsigned ii)	const	{ return pairString[ii]; }
	x_var		pValue			(unsigned ii)	const	{ return pairValue[ii]; }
	x_var	&	setPairValue	(unsigned ii)			{ return pairValue[ii]; }

	/* Clear all the stored information and data structures */
	void		clear			()						{
		siteMap.clear();
		siteString.clear();
		siteValue.clear();
		
		pairMap.clear();
		pairString.clear();
		pairValue.clear();
	}
	
	/* Preserve the datastructure and set all values to zero*/
	void		zerolize		()						{
		for (unsigned i=0; i<siteValue.size(); i++) { siteValue[i]=0; }
		for (unsigned i=0; i<pairValue.size(); i++) { pairValue[i]=0; }
	}
	
	double		max_abs			()						{
		double mabs = 0;
		for (unsigned i=0; i<siteValue.size(); i++) {
			double val = abs(siteValue[i]);
			if (mabs<val) mabs=val;
		}
		return mabs;
	}
	
	// ---- Save/load the order parameter in a unified formate.
	bool		save(string order_sub_name){

		auto siter = Lat.site_iteration();
		auto piter = Lat.pair_iteration();
		
		if (order_sub_name != "") {
			order_sub_name  = "."+order_sub_name  ;
		}
		
		ofstream out(Lat.filename+".ord"+order_sub_name  );
		
		// Storage of the site-order parameter
		if (siteString.size()>0) {
			if (out.is_open()) {
				out<<getSiteOrderStr(Lat);
			}
		}

		// Storage of the pair-order parameter
		if (pairString.size()>0) {
			if (out.is_open()) {
				out<<getPairOrderStr(Lat);
			}
		}
		
		out.close();
		return true;
	}
	bool		load(string order_sub_name, unsigned name_flag = 2){
		
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
			siteString.clear();
			siteValue.clear();
			pairMap.clear();
			pairString.clear();
			pairValue.clear();
	
            string line;
			string flag="";
			
            while ( getline(infile, line) ) {
				
				auto line_words = split(line, " ");
				
				// Read in the site order parameter
				if (line_words[0] == "#S") { flag = "#S"; continue;}
				if ( flag == "#S" and line_words.size() > 9 ) {
					
					for	(unsigned i=9; i<line_words.size(); i++) {
						auto str_variable = split(line_words[i], ":");
						auto var_name = str_variable[0];
						auto var = StrToComplex(str_variable[1]);
						
						auto sstr = line_words[1]+" "+var_name;
						siteString.push_back(sstr);
						siteValue.push_back(var);
						siteMap[sstr] = siteString.size()-1;
					}
				}
				
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
						
						auto var_bond = str_variable[0];			// bond variable ex: "-x"
						auto var_name = str_variable[1];			// name of the order ex: "hop"
						auto var = StrToComplex(str_variable[2]);	// value of the order ex: "(1.02,-2.001)"
						
						auto atomJ = Lat(atomI, var_bond);
						
						auto strIndexI = IntToStr(atomI.AtomIndex());
						auto strIndexJ = IntToStr(atomJ.AtomIndex());
						
						auto pairKeyStr = strIndexI+":"+strIndexJ+":"+var_bond+" "+var_name;
						//cout<<pairKeyStr<<endl;
						pairString.push_back(pairKeyStr);
						pairValue.push_back(var);
						pairMap[pairKeyStr] = pairString.size()-1;
					}
				}
			}
		}
		return true;
	}
	bool		importHoppingFromWannier90(string filename)	{
		
		ifstream infile(filename.c_str());
		
		if (infile.is_open()) {
			
			string line;
			while ( getline(infile, line)) {
				
				auto wordInLine = split(line," ");
				if (wordInLine.size() == 7) {
					
					string bondKey = "";
					for (int i=0 ; i<3 ; i++){
						if (wordInLine[i][0] != '-'){
							wordInLine[i] = "+" + wordInLine[i];
						}
						bondKey += wordInLine[i];
					}
					
					// Check if the bondKey exist in the Lat (Lattice) formate.
					//if (Lat.hasBondKey(bondKey)) {
					//	cout<<bondKey<<endl;
					//}
				}
				else{
					string warningStr= "Warning, formate not matched for w90:"+line+" of file:"+filename+"!";
					cout<<warningStr<<endl;
				}
			}
		}
		else {
			// Faild to open the file
			string errorStr = "Error, faild to open the file:"+filename+" !";
			cout<<errorStr<<endl;
			throw errorStr;
			return false;
		}
		
		infile.close();
		return true;
	}
};

OrderParameter operator+(const OrderParameter & L_val, const OrderParameter & R_val) {
	OrderParameter ret_val(L_val);
	
	if (
			L_val.siteMapSize()		== R_val.siteMapSize()		and
			L_val.siteStringSize()	== R_val.siteStringSize()	and
			L_val.siteValueSize()	== R_val.siteValueSize()	and
			L_val.pairMapSize()		== R_val.pairMapSize()		and
			L_val.pairStringSize()	== R_val.pairStringSize()	and
			L_val.pairValueSize()	== R_val.pairValueSize()
		) {
		
		for (unsigned i=0; i<ret_val.siteValueSize(); i++) {
			ret_val.setSiteValue(i) = L_val.sValue(i)+R_val.sValue(i);
		}
		
		for (unsigned i=0; i<ret_val.pairValueSize(); i++) {
			ret_val.setPairValue(i) = L_val.pValue(i)+R_val.pValue(i);
		}
	}
	return ret_val;
}

OrderParameter operator-(const OrderParameter & L_val, const OrderParameter & R_val) {
	OrderParameter ret_val(L_val);
	
	if (
			L_val.siteMapSize()	== R_val.siteMapSize() and
			L_val.siteStringSize() == R_val.siteStringSize() and
			L_val.siteValueSize()	== R_val.siteValueSize() and
			L_val.pairMapSize()	== R_val.pairMapSize() and
			L_val.pairStringSize() == R_val.pairStringSize() and
			L_val.pairValueSize()	== R_val.pairValueSize()
		) {
		
		for (unsigned i=0; i<ret_val.siteValueSize(); i++) {
			ret_val.setSiteValue(i) = L_val.sValue(i)-R_val.sValue(i);
		}
		
		for (unsigned i=0; i<ret_val.pairValueSize(); i++) {
			ret_val.setPairValue(i) = L_val.pValue(i)-R_val.pValue(i);
		}
	}
	return ret_val;
}

OrderParameter operator*(const OrderParameter & O_val, double d_val) {
	OrderParameter ret_val(O_val);
	
	for (unsigned i=0; i<ret_val.siteValueSize(); i++) { ret_val.setSiteValue(i) = ret_val.sValue(i)*d_val; }
	for (unsigned i=0; i<ret_val.pairValueSize(); i++) { ret_val.setPairValue(i) = ret_val.pValue(i)*d_val; }
	
	return ret_val;
}

OrderParameter operator*(double d_val, const OrderParameter & O_val) {
	OrderParameter ret_val(O_val);
	
	for (unsigned i=0; i<ret_val.siteValueSize(); i++) { ret_val.setSiteValue(i) = ret_val.sValue(i)*d_val; }
	for (unsigned i=0; i<ret_val.pairValueSize(); i++) { ret_val.setPairValue(i) = ret_val.pValue(i)*d_val; }
	
	return ret_val;
}














