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
//  tbm.cpp
//  TBM^3
//
//  Created by Yuan-Yen Tai on 7/23/15.
//

/* -------------------------------------------------------------------
 This is the entry class to construct any kind of TB model that depends
 on a specific lattice input file.

 1. The virtual method "init_order" will be called to initialize the 
    input order_parameters of a given lattice structure:
        [file name].lat
	with a order parameter file name:
        [file name].lat.ord.[order name]

 2. The virtual method "Hamiltonian" will be called to construct a given
   hopping, exchange interaction and chemical potentials in k-space.

 3. The virtual method "render" will be called for the designate calculations
   for example:
   a. Total electron calculation with the adjusted Chemical potential.
   b. Total energy.
   c. Band structure.
   d. Local density of states.
   e. others ...
 
 Please refer to the online documentation to see more specific examples.
 -------------------------------------------------------------------*/

class TBModel{
public:
	TBModel	(Lattice _L):		Lat(_L),
				initOrder			(Lat),	   
				newOrder			(Lat),	    
				initDenOrder		(Lat),	    
				newDenOrder			(Lat),	    
				hoppingOrder		(Lat){

		symmetry = _L.getSymmetry();
		Temperature = 0.0001;
		Mu = 0;
		alpha_Mu = 1;
		newDenMix = 0.5;
		
		variational_dt = 0.1;	// The time step for the variational method.
		isChemicalPotentialAdded = false;
		isAtomSiteSetted = false;
		isAtomPairSetted = false;
		site_iteration_index = -1; // The initial value '-1' is to avoid the indexing bug
		pair_iteration_index = -1; // The initial value '-1' is to avoid the indexing bug
		site_list = Lat.site_iteration();
		pair_list = Lat.pair_iteration();
		
		Ham	= x_mat(Lat.index_size()).get();
		
		//initOrder = OrderParameter(Lat);
		
		// Handle the parameters from the inputfile (filename).
		string filename = Lat.filename+".lif";
        ifstream infile(filename.c_str());
		if ( infile.is_open() ) {
			int flag=-1;
            string line;
            while ( getline(infile, line) ) {
                istringstream iss(line);
                string head;
				
                iss>>head;
				if (head[0]=='#')	{ flag=-1;}
				if (head=="#0")		{
					flag= 0;
					continue;
				}
				
				if (flag == 0) {
					replaceAll(line, " ", "");
					auto variable = split(line, "=");
					if (variable.size()==2) { parameter[variable[0]] = StrToComplex(variable[1]); }
				}
			}
		}
		infile.close();
	}
	
private:
	H_SYMMETRY			symmetry;			// Tell the system and Lattice which kind of symmetry to be used.
											//(1 NORMAL, 2 NAMBU, 3 EXNAMBU).
	// -----------------------------------
	// --- some dummy operational variable.
	// -----------------------------------
	Atom				temporaryAtomSite;
	AtomPair			temporaryAtomPair;
	bool				isAtomSiteSetted;
	bool				isAtomPairSetted;
	int					site_iteration_index;
	int					pair_iteration_index;
	vector<Atom>		site_list;
	vector<AtomPair>	pair_list;

	vector<Atom>		LDOSsites;

	// --- Data sets for band structure calculation.
	vector<r_mat>		kSpaceHighSymmetryPoints;
	vector<string>		kSpaceHighSymmetryLabels;
	
	// --- Data sets for periodic-supercell calculation.
	vector<r_mat>		kPoints;
	
	// --- Data sets for calculate gived total electron and chemical potential.
	// --- These variables will be called by : setElectronCarrier(...) and getCarrier(...)
	map<string,double>	ElectronCarrier;
	string				carrierList;
	
	bool	isChemicalPotentialAdded;
	
	// --------------------------------
	// To set the iteration informations within the site_iterate/pair_iterate function
	// --------------------------------
	void			setSite		(const Atom		& at)	{temporaryAtomSite = at; isAtomSiteSetted = true;}
	void			setPair		(const AtomPair	& ap)	{temporaryAtomPair = ap; isAtomPairSetted = true;}
	
	// --------------------------------
	// Use this VAR(...) function to investigate parameters from the user defined input file.
	// --------------------------------
private:	map<string, x_var>	parameter;
protected:
	x_var			VAR		(string str_key)			{
		double _sign=+1;
		
		if	(str_key[0]=='+') {
			_sign=+1;
			replaceAll(str_key, "+", "");
		}
		else if	(str_key[0]=='-'){
			_sign=-1;
			replaceAll(str_key, "-", "");
		}
		
		auto it=parameter.find(str_key);
		if (it!=parameter.end()){ return _sign*parameter[str_key]; }
		else					{
			string error_msg = "Error, cannot find parameter: "+str_key+", from input file: "+Lat.filename+".";
			cout<<error_msg<<endl;
			cout<<"Program terminated!"<<endl<<endl;
			throw error_msg; // Throw an error, if not catched, the program will be stoped.
		}
		return 0.0;
	} /*
	   Return the user defined variable. If user not defined this variable, return the default value.
	*/
	x_var			VAR		(string str_key, x_var val)	{
		double _sign=+1;
		
		if	(str_key[0]=='+') {
			_sign=+1;
			replaceAll(str_key, "+", "");
		}
		else if	(str_key[0]=='-'){
			_sign=-1;
			replaceAll(str_key, "-", "");
		}
		
		auto it=parameter.find(str_key);
		if (it!=parameter.end()){ return _sign*parameter[str_key]; }
		else					{
			//string warning_msg = "Warning, cannot find parameter: "+str_key+", from input file: "+Lat.filename+".";
			//cout<<warning_msg<<endl;
		}
		return val;
	} /*
	   Return the user defined variable. If user not defined this variable, return the default value.
	*/


protected:

	// ############################################
	//	Core part for the construction Hamiltonian and related operations
	// ############################################
	Lattice			Lat;			// Storage of the lattice structure
	x_mat			Ham;			// The matrix for the Hamiltonian
	vector<pair<r_mat, x_mat> >	KEigenVec;
	r_mat			k_space;		// The k-space vector will be passed by other methods.
	double			Temperature;	// Temperature, default value = 0.00001.
	double			Mu;				// Chemical potential, default value = 0.
	double			alpha_Mu;		// Fraction of increment of Chemical potential, default value = 1.
	double			variational_dt;	// The time step for the variational method.

	/* Use these iterate functions in the while loop to iterate through Sites/Pairs
	// Example:
	//			// site iterate
	//			while(site_iterate()){
	//				add_site("O 1d 1d +delta");
	//			}
	//			// pair iterate
	//			while(pair_iterate()){
	//				add_bond("O:1d    O:1d +x   hopping");
	//				add_bond("O-2:1d  O-2:1d +x hopping");
	//				add_bond("Fe:1d   O:1d +x   hopping");
	//			}
	*/
	bool			site_iterate()						{
		site_iteration_index ++ ;
		
		if (site_iteration_index < site_list.size()) {
			setSite(site_list[site_iteration_index]);
			return true;
		} else {
			site_iteration_index=-1;
		}
		return false;
	}
	bool			pair_iterate()						{
		pair_iteration_index ++ ;
		
		if (pair_iteration_index < pair_list.size()) {
			setPair(pair_list[pair_iteration_index]);
			return true;
		} else {
			pair_iteration_index=-1;
		}
		return false;
	}

	// --------------------------------
	// To get the iteration informations within the while(site_iterate())/while(pair_iterate()) loops
	// --------------------------------
	Atom			getSite		()						{return temporaryAtomSite;}
	AtomPair		getPair		()						{return temporaryAtomPair;}
	
	// --------------------------------
	// These function emulate operations to build the Hamiltonian
	// add_site		: Fe 1d 2d	(on-site operation
	// add_bond		: Mn:2d O:1d +x (non-hermition conjugate operation)
	// add_bond_hc	: Mn:2d O:1d +x (hermitian conjgate operation)
	// --------------------------------
	Element	get_site_element(const Atom	 & at,	  string opt, x_var val){
		
		Element e;
		
		e.opt=opt;
		auto word=split(opt, " ");
		
		if (word.size() !=3)	 return e(-1,-1,0, string("Wrong arguments in operation: ")+opt);
		
		if (word[0] != at.Name()) { return e(-1,-1,0,"non"); }
		if (at.find(word[1]) and at.find(word[2])) {
			int		I=at[word[1]];
			int		J=at[word[2]];
			e(I,J,val);
		} else {
			return e(-1,-1,0,"non");
		}
		
		return e;
	}
	Element	get_bond_element(const AtomPair & ap, string opt, x_var val){
		
		Element e;
		
		e.opt=opt;
		auto word=split(opt, " ");
		
		if (word.size() !=3)	 return e(-1,-1,0, string("Wrong arguments in operation: ")+opt);
		if (ap.bond!= Lat.translateBondString(word[2]) )	 return e(-1,-1,0, "non");
		
		
		auto AtI=split(word[0], ":");
		auto AtJ=split(word[1], ":");
		
		if ( AtI.size() !=2 or AtJ.size() !=2 )							return e(-1,-1,0, "non");
		
		// ------------
		bool flagAtom=false;
		
		if		( AtI[0] == ap.AtomI.Name()		and AtJ[0] == ap.AtomJ.Name())		flagAtom=true;
		else if ( AtI[0] == ap.AtomI.SubName()	and AtJ[0] == ap.AtomJ.SubName())	flagAtom=true;
		else if ( AtI[0] == ap.AtomI.Name()		and AtJ[0] == ap.AtomJ.SubName())	flagAtom=true;
		else if ( AtI[0] == ap.AtomI.SubName()	and AtJ[0] == ap.AtomJ.Name())		flagAtom=true;
		
		if ( flagAtom and ap.AtomI.find(AtI[1]) and ap.AtomJ.find(AtJ[1]) ){
			int		I=ap.AtomI[ AtI[1] ];
			int 	J=ap.AtomJ[ AtJ[1] ];
			auto	bvec=ap.bvec;
			
			// Apply value and the phase from k-space
			x_var phase = k_space[0]*bvec[0]+ k_space[1]*bvec[1]+ k_space[2]*bvec[2];
			auto	exp_IJ=exp(-Im*phase);
			
			e(I,J,val*exp_IJ);
			
		} else {
			return e(-1,-1,0,"non");
		}
		
		return e;
	}
	
	Element	add_site		(const Atom		& at, string opt, x_var val){
		
		Element e = get_site_element(at, opt, val);
		if (e.info == "success") {
			Ham(e.I,e.J)+=e.val;
		}

		return e;
	}
	Element	add_bond		(const AtomPair & ap, string opt, x_var val){
		
		Element e=get_bond_element(ap, opt, val);
		if (e.info=="success") {
			Ham(e.I,e.J)+=e.val;
		}
	
		return e;
	}
	Element	add_bond_hc		(const AtomPair & ap, string opt, x_var val){
	
		Element e=get_bond_element(ap, opt, val);
		if (e.info=="success") {
			Ham(e.I,e.J)+=e.val;
			Ham(e.J,e.I)+=e.val_conj;
		}
		
		return e;
	}
	Element	add_site		(const Atom		& at, string opt)			{
		Element e;
		
		e.opt=opt;
		
		auto word=split(opt, " ");
		
		if (word.size() !=4){ return e(-1,-1,0, string("Wrong arguments in operation: ")+opt); }
		else{
			string new_opt=word[0]+" "+word[1]+" "+word[2];
			e=add_site(at, new_opt, VAR(word[3]));
		}
		
		return e;
	}
	Element	add_bond		(const AtomPair & ap, string opt)			{
		Element e;
		
		e.opt=opt;
		
		auto word=split(opt, " ");
		
		if (word.size() !=4){ return e(-1,-1,0, string("Wrong arguments in operation: ")+opt); }
		else{
			string new_opt=word[0]+" "+word[1]+" "+word[2];
			e=add_bond(ap, new_opt, VAR(word[3]));
		}
		
		return e;
	}
	Element	add_bond_hc		(const AtomPair & ap, string opt)			{
		Element e;
		
		e.opt=opt;
		
		auto word=split(opt, " ");
		
		if (word.size() !=4){ return e(-1,-1,0, string("Wrong arguments in operation: ")+opt); }
		else{
			string new_opt=word[0]+" "+word[1]+" "+word[2];
			e=add_bond_hc(ap, new_opt, VAR(word[3]));
		}
		
		return e;
	}

	Element	add_site		(string opt, x_var val)						{return add_site   (getSite(), opt, val);}
	Element	add_bond		(string opt, x_var val)						{return add_bond   (getPair(), opt, val);}
	Element	add_bond_hc		(string opt, x_var val)						{return add_bond_hc(getPair(), opt, val);}
	Element	add_site		(string opt)								{return add_site   (getSite(), opt);}
	Element	add_bond		(string opt)								{return add_bond   (getPair(), opt);}
	Element	add_bond_hc		(string opt)								{return add_bond_hc(getPair(), opt);}
	
	x_var	q_var_site		(const Atom		& at, string opt, x_var val){ setSite(at); return q_var_site(opt, val); }
	x_var	q_var_bond		(const AtomPair	& ap, string opt, x_var val){ setPair(ap);  return q_var_bond(opt, val); }
	x_var	q_var_site_hc	(const Atom		& at, string opt, x_var val){ setSite(at); return q_var_site_hc(opt, val); }
	x_var	q_var_bond_hc	(const AtomPair	& ap, string opt, x_var val){ setPair(ap);  return q_var_bond_hc(opt, val); }
	r_var	c_var_site		(const Atom		& at, string opt, r_var val){ setSite(at); return c_var_site(opt, val); }
	r_var	c_var_bond		(const AtomPair	& ap, string opt, r_var val){ setPair(ap);  return c_var_bond(opt, val); }
	
	x_var	q_var_site		(string opt, x_var val)						{
		Element e=get_site_element(getSite(), opt, val);
		
		x_var	total=0;
		
		if (e.info == "success" and KEigenVec.size() > 0) {
			for (unsigned k=0; k<KEigenVec.size(); k++){	// Scan through k-space
				auto kEigen	= KEigenVec[k].first;
				auto kVec	= KEigenVec[k].second;
				for (unsigned i=0; i<kEigen.size(); i++) {
					double EigF = 1.0/(1+exp( kEigen[i]/Temperature));
					x_var v_I = kVec(e.I, i);
					x_var v_J = kVec(e.J, i);
					total+= EigF*conj(v_I)*v_J*val;
				}
			}
		}
		return total/KEigenVec.size();
	}
	x_var	q_var_bond		(string opt, x_var val)						{
		Element e=get_bond_element(getPair(), opt, val);
		
		x_var	total=0;
		
		if (e.info == "success" and KEigenVec.size() > 0) {
			for (unsigned k=0; k<KEigenVec.size(); k++){	// Scan through k-space
				auto kEigen	= KEigenVec[k].first;
				auto kVec	= KEigenVec[k].second;
				for (unsigned i=0; i<kEigen.size(); i++) {
					double EigF = 1.0/(1+exp( kEigen[i]/Temperature));
					x_var v_I = kVec(e.I, i);
					x_var v_J = kVec(e.J, i);
					total+= EigF*conj(v_I)*v_J*val;
				}
			}
		}
		return total;
	}
	r_var	q_var_site_hc	(string opt, x_var val)						{
		x_var tmp_val = q_var_site(opt, val);
		return 2*tmp_val.real();
	}
	r_var	q_var_bond_hc	(string opt, x_var val)						{
		x_var tmp_val = q_var_bond(opt, val);
		return 2*tmp_val.real();
	}
	r_var	c_var_site		(string opt, r_var val)						{
		auto si=getSite();
		if (opt == si.Name()) { return val; }
		return 0;
	}
	r_var	c_var_bond		(string opt, r_var val)						{
		auto pi=getPair();
		auto sI=pi.AtomI;
		auto sJ=pi.AtomJ;
		auto opt_word = split(opt, " ");
		if (opt_word[0] == sI.Name() and opt_word[1] == sJ.Name() and opt_word[2] == pi.bond) {
			return val;
		}
		return 0;
	}

	void	add_Chemical_Potential	()									{
		//---site iteration---
		auto siter = Lat.site_iteration();
		for (unsigned ii=0; ii<siter.size(); ii++) {
			auto si = siter[ii];
			
			// Add chemical potential
			auto label = si.index_label;
			for (unsigned i=0; i<label.size(); i++) {
				string opt = si.Name()+" "+label[i]+" "+label[i];
				add_site(si, opt, -Mu);
			}
		}
		isChemicalPotentialAdded = true;
	}

	// The operation list handling the correspond operations.
	vector<SiteOperationUnit> HundCouplingOperationList;
	vector<PairOperationUnit> SuperExchangeOperationList;
	vector<PairOperationUnit> DMInteractionOperationList;
	
	void	add_hund_spin			(string opt, double Jh, vector<x_var> SpinI)				{
		auto si		= getSite();
		auto word	= split(opt, " ");
		
		if (SpinI.size() == 3 and word[0] == si.Name()) {
			double Sx = Jh*SpinI[0].real();
			double Sy = Jh*SpinI[1].real();
			double Sz = Jh*SpinI[2].real();
			
			string sub_opt_ud = si.Name()+" "+word[1]+"u "+" "+word[1]+"d ";
			string sub_opt_du = si.Name()+" "+word[1]+"d "+" "+word[1]+"u ";
			string sub_opt_uu = si.Name()+" "+word[1]+"u "+" "+word[1]+"u ";
			string sub_opt_dd = si.Name()+" "+word[1]+"d "+" "+word[1]+"d ";
			
			add_site(sub_opt_ud, Sx-Im*Sy);
			add_site(sub_opt_du, Sx+Im*Sy);
			add_site(sub_opt_uu, Sz);
			add_site(sub_opt_dd,-Sz);
			HundCouplingOperationList.push_back(
										SiteOperationUnit(si, opt, Jh)
												);
		}
	}
	void	add_classical_spin		(string opt, r_var Js, vector<x_var> Si, vector<x_var> Sj)	{
		auto pit = getPair();
		auto siteI = pit.AtomI;
		auto siteJ = pit.AtomJ;
		
		auto word = split(opt, " ");
		
		if (siteI.Name() == word[0] and siteJ.Name() == word[1] and pit.bond == word[2]) {
			SuperExchangeOperationList.push_back(
										PairOperationUnit(pit, opt, Js)
												 );
		}
	}
	void	add_dm_spin		(string opt, r_var Js, vector<x_var> Si, vector<x_var> Sj)	{
		auto pit = getPair();
		auto siteI = pit.AtomI;
		auto siteJ = pit.AtomJ;
		
		auto word = split(opt, " ");
		
		if (siteI.Name() == word[0] and siteJ.Name() == word[1] and pit.bond == word[2]) {
			DMInteractionOperationList.push_back(
										PairOperationUnit(pit, opt, Js)
												 );
		}
	}
	
	// ############################################
	//	Extended calculation level
	// ############################################
	
	// --------------------------------
	// The variables for variational method.
	// --------------------------------
	OrderParameter	initOrder;		// The initialization of Order parameter.
	OrderParameter	newOrder;		// The output of new order parameter.
	OrderParameter	initDenOrder;	// The storage of electron density information.
	OrderParameter	newDenOrder;	// The storage of new-electron density information.
	double			newDenMix;		// The mix-weight for the new density.
	
	OrderParameter	hoppingOrder;	/* A special order that allow users to storage hopping terms
									 as an static order parameter */
	
	virtual	void	init_order		()=0;
	virtual void	Hamiltonian		()=0;							// The interface for user defined Hamiltonian.
	pair<r_mat,x_mat>	HamEvd		(r_mat kp)				{
		
		HundCouplingOperationList.clear();
		SuperExchangeOperationList.clear();
		
		k_space	= kp;
		Ham	= x_mat(Lat.index_size()).get();
		
		temporaryAtomSite = Atom();
		temporaryAtomPair = AtomPair();
		
		Hamiltonian();
		
		if (!isChemicalPotentialAdded) {
			cout<<"Warning ! You did not add the Chemical to the Hamiltonian,\
				please put \"add_Chemical_Potential()\" in Hamiltonian()."<<endl;
		}
		
		// Diagonalize the Hamiltonian and store in the value in Eigen and Vec.
		r_mat eig;
		x_mat vec;
		Ham.evd(eig, vec);
		return make_pair(eig,vec);
	
	}	// The executable function for diagonalization.

	// --------------------------------
	// Setup for periodic-supercell calculation.
	// --------------------------------
	void			add_k_point		(r_mat kpoint)			{ kPoints.push_back(kpoint); }
	void			clear_k_point	()						{
		kPoints.clear();
		KEigenVec.clear();
	}
	vector<double>	KHamEvd			()						{
		KEigenVec.clear();
		k_space = r_mat(1,3);
		Hamiltonian();
		
		double	max_eigen = -1000;
		double	min_eigen =  1000;
		
		time_t r_time;
		time(&r_time);
		printf ("The current local time is: %s", ctime (&r_time));
		
		for (unsigned i=0; i<kPoints.size(); i++) {
			auto eigvec = HamEvd(kPoints[i]);
			KEigenVec.push_back(eigvec);
			if (eigvec.first[0]						< min_eigen) { min_eigen = eigvec.first[0];				}
			if (eigvec.first[eigvec.first.size()-1]	> max_eigen) { max_eigen = eigvec.first[eigvec.first.size()-1];	}
		}
		
		time(&r_time);
		printf ("The current local time is: %s", ctime (&r_time));
		
		vector<double> return_variable_list;
		return_variable_list.push_back(min_eigen);
		return_variable_list.push_back(max_eigen);
		return_variable_list.push_back(calcTotalEnergy());
		
		return return_variable_list;
	}
	virtual	void	render			()=0;	// Use this virtual member to render the desired operation.
	
	// --------------------------------
	//	Calculate and get the density matrix D_ij
	// --------------------------------
	pair<int, int>	indexSite				(const Atom &At, string opt)			{
		auto word = split(opt, " ");
		
		int	index1=-2;
		int	index2=-2;
		
		if		(word.size()==2 and word[0] == At.Name()) {
			index1 = At[word[1]];
			index2 = index1;
		}else if(word.size()==3 and word[0] == At.Name()) {
			index1 = At[word[1]];
			index2 = At[word[2]];
		} else if (word.size()<2 or word.size() > 3) {
			cout<<"Warning, the option formate is wrong: "+opt;
		}
		return make_pair(index1, index2);
	}
	x_var			getDensityMatrix		(unsigned I, unsigned J)				{
		if ( KEigenVec.size()>0 )
		if (symmetry == NORMAL and I<Ham.cols() and J<Ham.rows()){
			x_var den_ij=0;
			for (unsigned k=0; k<KEigenVec.size(); k++) {
				auto kEig = KEigenVec[k].first;
				auto kVec = KEigenVec[k].second;
    
				for( int i=0 ; i<kEig.size() ; i++ ){
					double EigF = 1.0/(1+exp( kEig[i]/Temperature));
					x_var vI = kVec(I,i);
					x_var vJ = kVec(J,i);
					den_ij += conj(vI)*vJ*EigF;
				}
			}
			return den_ij/KEigenVec.size();
		}
		
		return 0;
	}
	x_var			getSiteDensityMatrix	(const Atom &At, string opt)			{
		auto site_index = indexSite(At, opt);
		auto val = getDensityMatrix(site_index.first, site_index.second);
		return val;
	}
	void			selectLDOSsite			(unsigned index)						{
		while (site_iterate()) {
			auto si = getSite();
			if (si.AtomIndex() == index and
				(si.Name()!="VA" or si.Name()!="VC" or si.Name()!="BD")
				)
			{
				cout<<"Selected ldos atom-index: "<<si.Name()<<"-"<<index<<", ";
				LDOSsites.push_back(si);
				break;
			}
		}
		cout<<endl;
	}
	void			select_All_LDOSsite		()										{
		while (site_iterate()) {
			auto si = getSite();
			if (si.index_label.size()>0) {
				LDOSsites.push_back(si);
			}
		}
	}
	void			calcLDOS				(double dE, double Gamma)				{
		
		auto var_list= KHamEvd();

		
		if (KEigenVec.size() > 0) {
			cout<<endl;
			cout<<"Calculating the LDOS, Mu = "<<Mu<<endl;
			
			double	minE = var_list[0]-4*dE;
			double	maxE = var_list[1]+4*dE;
			cout<<"# Energy range"<<endl;
			cout<<"From:"<<minE<<"    To:"<<maxE<<endl;
			
			
			vector<pair<double,double> > LDOS;
			for (double w=minE; w<maxE; w+=dE){
				auto _ldos = make_pair(w, 0.0);
				LDOS.push_back(_ldos);
				
				for (unsigned i=0; i<LDOSsites.size(); i++) {
					LDOSsites[i].LDOS.push_back(_ldos);
				}
			}
			
			for (unsigned iw=0; iw<LDOS.size(); iw++) {
				double w		=LDOS[iw].first;
				cout<<fmt(w,8)<<" ";
				if ( (iw+1) % 10 == 0) { cout<<endl; }
				
				for (unsigned k=0; k<KEigenVec.size(); k++) {
					auto kEig = KEigenVec[k].first;
					auto kVec = KEigenVec[k].second;
					
					for ( int n=0; n<kEig.size(); n++) {
						double e_w		= kEig[n]-w;
						double delta	= Gamma/(pi*(e_w*e_w + Gamma*Gamma));
						
						for (unsigned site_index=0; site_index<LDOSsites.size(); site_index++) {
							
							auto si = LDOSsites[site_index];
							auto index_label = LDOSsites[site_index].index_label;
							for (unsigned ll_index=0; ll_index<index_label.size(); ll_index++) {
								auto site_label_index = indexSite(si, si.Name()+" "+index_label[ll_index]);
								
								x_var u_n		= kVec(site_label_index.first, n);
								LDOSsites[site_index].LDOS[iw].second	+= (conj(u_n)*u_n).real() * delta/KEigenVec.size();
								
							}
							
						}
					}
				}
			}
		}
		cout<<endl;
		
		ofstream out(Lat.filename+".ldos.py");
		if (LDOSsites.size()>0) {
			auto si=LDOSsites[0];
			out<<"import numpy as np"<<endl;
			out<<"import matplotlib.pyplot as plt"<<endl;
			out<<"x=[";
			for (unsigned i=0 ; i<si.LDOS.size(); i++) { out<<fmt(si.LDOS[i].first)<<", "; }
			out<<"]"<<endl<<endl;

			// Output of each LDOS component
			for (unsigned ii=0; ii<LDOSsites.size(); ii++) {
				auto si=LDOSsites[ii];

				out<<"y"+IntToStr(si.AtomIndex())+"=[";
				for (unsigned i=0 ; i<si.LDOS.size(); i++) { out<<fmt(si.LDOS[i].second)<<", "; }
				out<<"]"<<endl;
			}

			out<<endl;
			// Plot, using matplotlib.
			for (unsigned ii=0; ii<LDOSsites.size(); ii++) {
				auto si=LDOSsites[ii];
				out<<"plt.plot(x, y"+IntToStr(si.AtomIndex())+", '-', linewidth=2)"<<endl;
			}
			out<<"plt.show()"<<endl;
		}
		out.close();
	}

	// --------------------------------
	// -Calculate the chemical potential and electron filling.
	// --------------------------------
	void			setElectronCarrier		(string inputs)							{
		ElectronCarrier.clear();
		ElectronCarrier["VA"]=0; // VA is for virtual atom.
		ElectronCarrier["VC"]=0; // VC is for vacancy.
		ElectronCarrier["BD"]=0; // BD is for boundary.
		
		carrierList = inputs;
		
		auto word = split(inputs, " ");
		for (unsigned i=0; i<word.size(); i++) {
			auto ss = split(word[i], ":");
			ElectronCarrier[ss[0]] = StrToDouble(ss[1]);
		}
	}
	void			initElectronDensity		(string inputs)							{
		map<string,double>	InitDensityList;
		InitDensityList.clear();
		InitDensityList["VA"]=0; // VA is for virtual atom.
		InitDensityList["VC"]=0; // VC is for vacancy.
		InitDensityList["BD"]=0; // BD is for boundary.
		
		auto denlist = inputs;
		
		auto word = split(denlist, " ");
		for (unsigned i=0; i<word.size(); i++) {
			auto ss = split(word[i], ":");
			InitDensityList[ss[0]] = StrToDouble(ss[1]);
		}
		
		while (site_iterate()) {
			auto si = getSite();
			initDenOrder(si, si.Name()+" den") = InitDensityList[si.Name()];
		}
		//cout<<initDenOrder.getSiteOrderStr(Lat)<<endl;
		newDenOrder = initDenOrder;
	}
	double			getCarrier				(string inputs)							{
		
		auto it = ElectronCarrier.find(inputs);
		if (it != ElectronCarrier.end()){ return ElectronCarrier[inputs]; }
		else							{ cout<<"Warning, no definition for: \""+inputs+"\" in ["+ carrierList+"]"<<endl; }
		return 0;
	}
	double			calculateTotalElectron	()										{
		
		if (ElectronCarrier.size() == 0) {
			cout<<"Warning, carrier density not defined!\nPlease add setElectronCarrier(\"Atom1:N1 Atom2:N2 ...\") in your code for the use of calculateTotalElectron(). \n";
			return 0;
		}

		double totalElectron=0;
		auto siter = Lat.site_iteration();
		
		for (unsigned i=0; i<siter.size(); i++) {
			totalElectron += getCarrier( siter[i].Name() );
		}
		
		return totalElectron;
	}

	// --------------------------------
	//	Calculate the electron density for each site
	// --------------------------------
	x_var			calcDensity(const Atom & At, string line, double _mu=0)			{
		auto index = indexSite(At, line);
		if ( KEigenVec.size() > 0){
			if (symmetry == NORMAL and index.first >= 0)
			if (At.Name() != "VA" and At.Name() != "VC" and At.Name() != "BD")
			{
				x_var den_ii=0;
				
				for (unsigned k=0; k<KEigenVec.size(); k++) {
					
					auto	kEig = KEigenVec[k].first;
					auto	kVec = KEigenVec[k].second;
					for( int i=0 ; i<kEig.size() ; i++ ){
						double EigF = 1.0/(1+exp( (kEig[i]-_mu)/Temperature));
						x_var v = kVec(index.first,i);
						den_ii += conj(v)*v*EigF;
					}
					
				}
				return den_ii/KEigenVec.size();
			}
		}
		
		return 0;
	}
	OrderParameter	calcElectronFilling		(double _mu=0)							{
		
		OrderParameter order(Lat);
		while (site_iterate()) {
			auto si = getSite();
			auto label = si.index_label;
			order(si,si.Name()+" den")=0;
			
			for (unsigned i=0; i<label.size(); i++) {
				string opt = si.Name()+" "+label[i];
				order(si,si.Name()+" den")+=calcDensity(si, opt, _mu);
			}
		}
		
		return order;
	}
	OrderParameter	calcPolarizedSpin		()										{
		//// Construct the frame of the electron density
		auto siter = Lat.site_iteration();
		
		OrderParameter order(Lat);
		for (unsigned i=0; i<siter.size(); i++) {
			auto si = siter[i];
			
			auto atom_orb_info = si.OrbitalInfo();
			int	orb_number = atom_orb_info.first;
			string spin_degree = atom_orb_info.second;
			if (spin_degree == "s" and
				(si.Name()!="VA" or si.Name()!="VC" or si.Name()!="BD")
				)
			for (unsigned orb_i=0; orb_i<orb_number; orb_i++) {
				string orb = IntToStr(orb_i+1);
				x_var C1uu = getSiteDensityMatrix(si, si.Name()+" "+orb+"u "+orb+"u");
				x_var C1dd = getSiteDensityMatrix(si, si.Name()+" "+orb+"d "+orb+"d");
				x_var C1ud = getSiteDensityMatrix(si, si.Name()+" "+orb+"u "+orb+"d");
				x_var C1du = getSiteDensityMatrix(si, si.Name()+" "+orb+"d "+orb+"u");
				            
				
				x_var psx = C1ud.real()+C1du.real();
				x_var psy = C1ud.imag()-C1du.imag();
				x_var psz = C1uu.real()-C1dd.real();
				order(si, si.Name()+" "+orb+"psx")=psx;
				order(si, si.Name()+" "+orb+"psy")=psy;
				order(si, si.Name()+" "+orb+"psz")=psz;
				
				double n = abs(C1uu+C1dd);
				double m = abs(psx*psx + psy*psy + psz*psz);
				order(si, si.Name()+" "+orb+"pol.elec1") = (n-m)/2;
				order(si, si.Name()+" "+orb+"pol.elec2") = (n+m)/2;
			}
		}
		
		order.save("p-spin");
		
		return order;
	}

	// --------------------------------
	//	Calculate the LLG for the spin configuration
	// --------------------------------
	void			setProgramControler		()										{
		ofstream	out(Lat.filename+".control");
		out<<"continue";
		out.close();
	}
	string			getProgramControler		()										{
		ifstream	infile(Lat.filename+".control");
		if (infile.is_open()) {
			string line;
			getline(infile, line);
			return line;
		}
		infile.close();
		return "continue";
	}
	
	void			setVariational			()										{
		
		newOrder.clear();			// Zerolize the newOrder and use it for Variational() calculation
		variational_dt = VAR("vdt").real();
		
		map<int, r_mat> Field1, Field2, Field3;
		while (site_iterate()) {
			Field1[getSite().AtomIndex()] = r_mat(1,3);
			Field2[getSite().AtomIndex()] = r_mat(1,3);
			Field3[getSite().AtomIndex()] = r_mat(1,3);
		} // Construct the force term
		
		// Field term for Hund coupling
		// Construct the Field-Term for variational method of Quantum Hund coupling sping.
		for (unsigned i=0; i<HundCouplingOperationList.size(); i++)		{
			auto HundOpt = HundCouplingOperationList[i];
			auto si = HundOpt.SiteI;
			auto Jh = HundOpt.val.real();
			auto opt= HundOpt.opt;
			
			auto word = split(opt, " ");
			
			string sub_opt_ud = word[0]+" "+word[1]+"u "+" "+word[1]+"d ";
			string sub_opt_du = word[0]+" "+word[1]+"d "+" "+word[1]+"u ";
			string sub_opt_uu = word[0]+" "+word[1]+"u "+" "+word[1]+"u ";
			string sub_opt_dd = word[0]+" "+word[1]+"d "+" "+word[1]+"d ";
			
			r_mat F(1,3);
			F(0,0)  += Jh*q_var_site(si, sub_opt_ud, 1  ).real(); // \sigma_x
			F(0,0)  += Jh*q_var_site(si, sub_opt_du, 1  ).real(); // \sigma_x
			F(0,1)  += Jh*q_var_site(si, sub_opt_ud,-Im ).real(); // \sigma_y
			F(0,1)  += Jh*q_var_site(si, sub_opt_du, Im ).real(); // \sigma_y
			F(0,2)  += Jh*q_var_site(si, sub_opt_uu, 1  ).real(); // \sigma_z
			F(0,2)  += Jh*q_var_site(si, sub_opt_dd,-1  ).real(); // \sigma_z
			Field1[si.AtomIndex()] = Field1[si.AtomIndex()] + F;
		}
		
		// Field term for Super-Exchange
		// Construct the Field-Term for variational method of Super Exchange coupling sping.
		for (unsigned i=0; i<SuperExchangeOperationList.size() ; i++)	{
			auto SuperExOpt = SuperExchangeOperationList[i];
			auto si = SuperExOpt.pit.AtomI;
			auto sj = SuperExOpt.pit.AtomJ;
			auto Js = SuperExOpt.val.real();
			auto opt= SuperExOpt.opt;
			
			r_mat Sj(1,3);
			auto FSi = initOrder.getVars(si, "Sx Sy Sz");
			auto FSj = initOrder.getVars(sj, "Sx Sy Sz");
			if (FSi.size() == 3 and FSj.size() == 3){
				Sj(0,0) = Js*FSj[0].real();
				Sj(0,1) = Js*FSj[1].real();
				Sj(0,2) = Js*FSj[2].real();
				Field2[si.AtomIndex()] = Field2[si.AtomIndex()] + Sj;
			}
		}
		
		// Field term for DM-interaction
		/*
		for ( ... ) {
			...
		}
		*/
		
		// To add up the Field term and translate to Force term and update the spin configuration.
		while (site_iterate()) {
			auto siteI = getSite();
			r_mat F1 =-1.0*Field1[siteI.AtomIndex()];
			r_mat F2 =-1.0*Field2[siteI.AtomIndex()];
			//r_mat F3 =-1.0*Field3[siteI.AtomIndex()];
			
			r_mat Si(1,3), Snew(1,3);
			r_mat ForceAll;
			
			Si(0,0) = initOrder(siteI, siteI.Name()+" Sx").real();
			Si(0,1) = initOrder(siteI, siteI.Name()+" Sy").real();
			Si(0,2) = initOrder(siteI, siteI.Name()+" Sz").real();
			
			auto dummy_check= sqrt(cdot(Si,Si));
			
			// Dummy detection for the spin operation.
			// If dummy_check is not zero then this is a real spin site.
			if ( dummy_check > 0.00000001) {
				
				// --- Backup for the Gaussian-white-noise operation.
				//r_mat damp(1,3);
				//damp[0] = sqrt(2*Temperature)*gaussian_white();
				//damp[1] = sqrt(2*Temperature)*gaussian_white();
				//damp[2] = sqrt(2*Temperature)*gaussian_white();
				
				ForceAll = curl(Si, F1+F2);				// Force
				//ForceAll = curl(Si, F1+F2+F3);				// Force
				
				r_mat FdampAll=curl( 1.0*ForceAll, Si); // damping
				
				Snew =Si+ForceAll*variational_dt		// The force term
						+2*FdampAll*variational_dt;	// The damping term
				
				auto denumerator = sqrt(cdot(Snew,Snew));
			
				Snew =Snew * (1/denumerator);	// Renormalize the new spin
			}
			
			newOrder(siteI,siteI.Name()+" Sx") = Snew(0,0);
			newOrder(siteI,siteI.Name()+" Sy") = Snew(0,1);
			newOrder(siteI,siteI.Name()+" Sz") = Snew(0,2);
		}
	}
	void			calcVariation			(unsigned iterations=3000)				{
		cout<<" Calculating LLG Spin Dynamics: "<<endl;

		double		den_diff			=100;
		double		diff				=100;
		unsigned int iter= 0;
		
		double		alpha = VAR("alpha", 0).real();
		
		setProgramControler(); // Initialize the flow of the program.
		
		double		__den_diff	= VAR("den_diff",0.01).real();
		double		__spin_diff	= VAR("spin_diff",0.002).real();

		cout<<"Den_tolerance: "<<__den_diff<<endl
			<<"Spin_tolerance: "<<__spin_diff<<endl;

		
		bool is_init_mu=true;
		while( abs(den_diff) > __den_diff or diff > __spin_diff ){
			
			
			// Read in the controle flow of the program.
			if (getProgramControler() != "continue") { break; }
			if (iter>iterations) break;
			if (iter!=0) is_init_mu = false;
			
			initDenOrder.save("den.ini");
			
			auto var_list = calcChemicalPotential(is_init_mu);

			// Calculation for electron density.
			auto diffDenOrder = initDenOrder - newDenOrder;
			den_diff = diffDenOrder.max_abs();
			
			if (abs(alpha) < 0.000001) den_diff = 0;
			
			auto total_energy = var_list[2];
			
			// ***********************************************
			// Calculation for LLG spin.
			// ***********************************************
			if (abs(den_diff) < __den_diff) {	// See if the chemical potential (Mu) has been adjusted to the right filling.
				var_list = KHamEvd();			// Diagonalize the Hamiltonian in k-space
				total_energy = var_list[2];
				setVariational();				// Calculate and store new value in the newOrder
				
				calcPolarizedSpin();

				auto diffOrder = initOrder - newOrder; // Calculate the difference of order
				diff = diffOrder.max_abs();
				initOrder = newOrder;		// Assign newOrder to the initial Order
				
				// Print information
				cout<<" Loop:"		 <<fmt(iter, 4)
					<<" Mu:"		 <<fmt(Mu,16)
					<<" den_diff:"	 <<fmt(den_diff,16)
					<<" diff:"		 <<fmt(diff, 20)
					<<" TotalE:"	 <<fmt(total_energy)<<endl;
			}
			else{
				// Print information
				cout<<" Loop:"		 <<fmt(iter, 4)
					<<" Mu:"		 <<fmt(Mu,16)
					<<" den_diff:"	 <<fmt(den_diff,16)
					<<" diff:"		 <<fmt(diff, 20)<<endl;
			}
			
			initDenOrder = newDenMix*newDenOrder+(1-newDenMix)*initDenOrder;
			
			saveElectronFilling(initDenOrder);
			initOrder.save("");
			iter++;
		}
	}
	vector<double>	calcChemicalPotential	(bool is_uncertain_Mu=true)				{
		//cout<<endl;
		double	totalElectron		= calculateTotalElectron();
		auto	efOrder				= calcElectronFilling();
		double	n_diff				= 1;
		double	__mu_diff			= VAR("mu_diff",0.0001).real();
		
		
		double	destMu = Mu;
		Mu = 0;
		auto	var_lis = KHamEvd();
		double	E_min = var_lis[0];
		double	E_max = var_lis[1];
		if (is_uncertain_Mu) { destMu = (E_max+E_min)/2; }
		
		
		unsigned len = 10;
		cout<<"  -Calculating chemical potential \"Mu\" in range:"<<E_min<<" "<<E_max<<endl;
		cout<<"    "<<fmt("Mu",len)<<" "<<fmt("Dest e-den",len)<<"   "<<fmt("True e-den",len)<<"   "<<fmt("tot_den_diff",len)<<endl;
		

		unsigned steps = 0;
		
		do{
			efOrder = calcElectronFilling(destMu);
			n_diff  = totalElectron - efOrder.totalSiteOrder().real();
			
			cout<<"    "<<fmt(destMu,len)<<" "<<fmt(totalElectron,len) <<" - "<<fmt( totalElectron-n_diff,len )<<" = "<<fmt(n_diff,len) <<endl;
			
			if (abs(n_diff) > abs(__mu_diff) ) {
				if (n_diff > 0) { E_min = destMu; }
				if (n_diff < 0) { E_max = destMu; }
				destMu = (E_max+E_min)/2;
				steps++;
			}
		
		}while ( abs(n_diff) > abs(__mu_diff) );
		
		
		Mu = destMu;
		cout<<endl;
		newDenOrder = efOrder;
		newDenOrder.save("den.calc");
		return var_lis;
	}
	void			saveElectronFilling		(OrderParameter eOrder)					{
		ofstream out(Lat.filename+".den");
		out<<"#0 Chemical Potential"<<endl;
		out<<fmt(Mu)<<endl;
		out<<eOrder.getSiteOrderStr(Lat);
		out.close();
	}
	void			loadElectronFilling		()										{
		ifstream infile(Lat.filename+".den");
		if (infile.is_open()) {
			string line, control_flag="";
			while ( getline(infile, line) ) {
				auto words = split(line, " ");
				if (words[0] == "#0")	{ control_flag = "Chemical_Potential"; continue;}
				else if (words[0][0] == '#') { break; }
				if (control_flag == "Chemical_Potential") { Mu = StrToDouble(line); continue; }
			}
			initDenOrder.load("den", 1);
		}
		infile.close();
	}

	// --------------------------------
	// -Calculate Total energy
	// --------------------------------
	double			calcTotalEnergy			()										{

		double totalE=0;
		for (unsigned k=0; k<KEigenVec.size(); k++) {
			auto kEig = KEigenVec[k].first;
			for (unsigned i=0; i<kEig.size(); i++) {
				double EigF = 1.0/(1+exp( kEig[i]/Temperature));
				totalE+=EigF*kEig[i];
			}
		}
		
		totalE = totalE / KEigenVec.size();
		
		// Construct the Field-Term for variational method of Super Exchange coupling sping.
		for (unsigned i=0; i<SuperExchangeOperationList.size() ; i++)	{
			auto SuperExOpt = SuperExchangeOperationList[i];
			auto si = SuperExOpt.pit.AtomI;
			auto sj = SuperExOpt.pit.AtomJ;
			auto Js = SuperExOpt.val.real();
			auto opt= SuperExOpt.opt;
			
			r_mat Sj(1,3);
			auto FSi = initOrder.getVars(si, "Sx Sy Sz");
			auto FSj = initOrder.getVars(sj, "Sx Sy Sz");
			if (FSi.size() == 3 and FSj.size() == 3){
				
				for (unsigned k=0; k<FSi.size(); k++) {
					totalE += Js*FSi[k].real()*FSj[k].real();
				}
			}
		}
		
		// Construct the DM Exchange coupling sping.
		//for (unsigned i=0; i<DMInteractionOperationList.size() ; i++)	{
		//	auto DMInteraction= DMInteractionOperationList[i];
		//	auto si = DMInteraction.pit.AtomI;
		//	auto sj = DMInteraction.pit.AtomJ;
		//	auto Js = DMInteraction.val.real();
		//	auto opt= DMInteraction.opt;
		//	
		//	r_mat Sj(1,3);
		//	auto FSi = initOrder.getVars(si, "Sx Sy Sz");
		//	auto FSj = initOrder.getVars(sj, "Sx Sy Sz");
		//	if (FSi.size() == 3 and FSj.size() == 3){
		//		
		//		for (unsigned k=0; k<FSi.size(); k++) {
		//			//totalE += Js*FSi[k].real()*FSj[k].real();
		//		}
		//	}
		//}
		
		// Construct the DM interaction coupling sping.
		//for ..
		
		string	filename = Lat.filename+".eng";
		ofstream out(filename);
		out<<totalE;
		out.close();
		
		return totalE;
	}
	
	// --------------------------------
	// -Calculate and save the band structure.
	// --------------------------------
	void			clear_ksymm_point		()										{ kSpaceHighSymmetryPoints.clear(); }
	void			add_ksymm_point			(string kPointLabel, r_mat kpoint)		{
		if (kpoint.cols()==1 and kpoint.rows()==3) {
			kSpaceHighSymmetryPoints.push_back(kpoint);
			kSpaceHighSymmetryLabels.push_back(kPointLabel);
		}
	}
	void			calculateBandStructure	(string filename="", unsigned Nsteps=10){
		cout<<"Calculating the Band structure ..."<<endl;
		
		// Manage filename
		if (filename.size()==0) {
			filename = Lat.filename+".ban";
		} else{
			filename = filename+".ban";
		}
		
		vector<r_mat>	highSymmetryLine;
		
		// Construct the k-space high-symmetry-line from kSpaceHighSymmetryPoints.
		if (kSpaceHighSymmetryPoints.size()>=2) {
			
			for (unsigned i=0; i<kSpaceHighSymmetryPoints.size()-1; i++) {
				
				unsigned j=i+1;
				auto	partOfHighSymmetryLine =
					make_line(kSpaceHighSymmetryPoints[i], kSpaceHighSymmetryPoints[j], Nsteps);
				
				for (unsigned ii=0; ii<partOfHighSymmetryLine.size(); ii++) {
					highSymmetryLine.push_back(
						partOfHighSymmetryLine[ii]
					);
				}
			}
		}
		
		if (highSymmetryLine.size()>0 and filename.size()>0) {
			ofstream out(filename);
			
			for (unsigned i=0; i<highSymmetryLine.size(); i++) {
				auto kp = highSymmetryLine[i];
				auto eigvec = HamEvd(kp);
				eigvec.first.setPrintLength(16);
				out<< kp <<" "<<eigvec.first<<endl;
			}
		
			out.close();
		}
	}
};




