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
//  TBModelBase.hpp
//  TBM^3
//

// -----------------------------------
// Make a linear points between two different k-points.
// -----------------------------------
vector<r_mat> make_line(r_mat point_a, r_mat point_b, unsigned Nsteps=10){
	vector<r_mat> ll;
	for (double t=0; t<1; t+=1.0/Nsteps) {
		r_mat position= (point_b-point_a)*t+point_a;
		ll.push_back(position);
	}
	
	return ll;
}

class TBModelBase{
	
public:
	TBModelBase(string filename):
		Lat				(filename),
		spinNormalLat	(filename, "on", "normal"),
		TBD				(Lat),
		spinNormalTBD	(spinNormalLat),
		tbd				(TBD),
		stbd			(spinNormalTBD	)
	{	}
	
protected:
	Lattice	Lat;
	Lattice	spinNormalLat;
	
private:
	TBDataSource TBD;
	TBDataSource spinNormalTBD;
	
protected:
	TBDataSource & tbd;
	TBDataSource & stbd;
	
	
	vector<pair<string,r_mat> >	kSpaceHighSymmetryPoints;
	
	bool iterate(){ return Lat.iterate(); }
	
	virtual void initOrder()	= 0;
	virtual void Hamiltonian()	= 0;
	virtual void run()			= 0;
	
	void		constructHam	(TBDataSource & rtbd, bool withMu = true)					{
		// Construct the Hamiltonian from the system-build structure: TBDataSource (tbd).
		
		if( Lat.parameter.VAR("disable_quantum", 0).real() == 0 ){
			rtbd.constructTBMHam();
			if( withMu ) rtbd.addChemicalPotential(rtbd.Lat.parameter.VAR("Mu").real());
		}
		
		// Construct the customized Hamiltonian.
		Hamiltonian();
	}

	/*-----------------------------------------------------
	 A set of operations for the Band structure calculation.
	 ------------------------------------------------------*/
	void	clear_ksymm_point		()														{ kSpaceHighSymmetryPoints.clear(); }
	void	add_ksymm_point			(string kPointLabel, r_mat kpoint)						{
		if (kpoint.cols()==1 and kpoint.rows()==3) {
			kSpaceHighSymmetryPoints.push_back(make_pair(kPointLabel, kpoint));
		}
	}
	void	calculateBandStructure	(TBDataSource & rtbd, unsigned Nsteps=50)				{
		if( Lat.parameter.VAR("disable_quantum", 0).real() != 0 ){
			cout<< "Warning, due to flag 'disable_quantum' turned on."
				<< "The band structure calculation will be ignored."<<endl;
			return;
		}
		
		auto B = Lat.basisVector.getBVec();
		auto & b1 = B[0]*0.5;
		auto & b2 = B[1]*0.5;
		auto & b3 = B[2]*0.5;
			
		for(auto & kp: Lat.kSymmPointParser.kSymmPointList){
			add_ksymm_point(kp.first, +kp.second[0]*b1 +kp.second[1]*b2 +kp.second[2]*b3);
		}
		
		vector<pair<string, r_mat> >	highSymmetryLine;
		
		// Construct the k-space high-symmetry-line from kSpaceHighSymmetryPoints.
		if (kSpaceHighSymmetryPoints.size()>=2) {
			
			for (unsigned i=0; i<kSpaceHighSymmetryPoints.size()-1; i++) {
				
				unsigned j=i+1;
				auto	partOfHighSymmetryLine =
					make_line(kSpaceHighSymmetryPoints[i].second, kSpaceHighSymmetryPoints[j].second, Nsteps);
				
				for (unsigned ii=0; ii<partOfHighSymmetryLine.size(); ii++) {
					string label = ii==0 ? kSpaceHighSymmetryPoints[i].first : "-";
					highSymmetryLine.push_back(
						make_pair(label, partOfHighSymmetryLine[ii])
					);
				}
			}
		}
		
		/* Store the band structure into 'xxx.lat.ban' */
		string	filename = Lat.FileName()+".ban";
		ofstream out(filename);
		
		constructHam(rtbd);
		for ( auto kp: highSymmetryLine){
			auto & evv = rtbd.HamEvd(kp.second);
			evv.eigenVector.setPrintLength(16);
			out<< kp.first<<" "<< kp.second <<" "<<evv.eigenValue<<endl;
		}
		
		out.close();
	}
	
	/*-----------------------------------------------------
	 Diagonalize the full Hamiltonian and store the results in KEigenValVec.
	 ------------------------------------------------------*/
	void	KHamEvd					(TBDataSource & rtbd, bool withMu = true)				{
		
		constructHam(rtbd, withMu);
		
		if( Lat.parameter.VAR("disable_quantum", 0).real() != 0 ) return;
		
		auto Nb = Lat.parameter.VEC("Nb");
		
		if( Nb.size() == 1){
			auto B = Lat.basisVector.getBVec();
			auto & b1 = B[0]*0.5;
			auto N1 = Nb[0];
			
			for (double i1=-N1; i1<N1; i1++){
				auto kPoint = (i1/N1)*b1;
				rtbd.KEigenValVec.push_back(rtbd.HamEvd(kPoint));
			}
		}
		if( Nb.size() == 2){
			auto B = Lat.basisVector.getBVec();
			auto & b1 = B[0]*0.5;
			auto & b2 = B[1]*0.5;
			auto N1 = Nb[0], N2 = Nb[1];
			for (double i1=-N1; i1<N1; i1++)
			for (double i2=-N2; i2<N2; i2++){
				auto kPoint = (i1/N1)*b1 + (i2/N2)*b2;
				rtbd.KEigenValVec.push_back(rtbd.HamEvd(kPoint));
			}
		}
		if( Nb.size() == 3){
			auto B = Lat.basisVector.getBVec();
			auto & b1 = B[0]*0.5;
			auto & b2 = B[1]*0.5;
			auto & b3 = B[2]*0.5;
			auto N1 = Nb[0], N2 = Nb[1], N3 = Nb[2];
			for (double i1=-N1; i1<N1; i1++)
			for (double i2=-N2; i2<N2; i2++)
			for (double i3=-N3; i3<N3; i3++) {
				auto kPoint = (i1/N1)*b1 + (i2/N2)*b2 + (i3/N3)*b3;
				rtbd.KEigenValVec.push_back(rtbd.HamEvd(kPoint));
			}
		}
	}
	
	/*-----------------------------------------------------
	 Calculate total electron / chemical potential
	 ------------------------------------------------------*/
	double	countTotalElectron		()														{
		double totalElectron = 0;
		while(iterate()){
			auto atomI = Lat.getAtom();
			totalElectron += Lat.coreCharge.getCharge(atomI.atomName);
		}
		return totalElectron;
	}
	r_var	calculateElectronFilling(TBDataSource & rtbd, double Mu=0)						{
		
		r_var totalDen = 0;
		while( rtbd.Lat.iterate() ){
			
			auto atomI = rtbd.Lat.getAtom();
			auto allSubIndex = atomI.allIndexList();
			
			r_var atomDen = 0;
			for(auto & subI: allSubIndex){
				
				auto label = split(subI.first, ".");
				char space_id = label[1][0];
				x_var den = 0;
				if( space_id == 'u' or space_id == 'd' or space_id == 'n' or space_id == 'A' )
					den = rtbd.getDensityMatrix(subI.second, subI.second, Mu);
				
				if( space_id == 'B' )
					den = 1-rtbd.getDensityMatrix(subI.second, subI.second, Mu);
				
				atomDen += den.real();
			}
			switch( rtbd.Lat.HSpace() ){
				case NORMAL: if( rtbd.Lat.parameter.STR("spin") == "off")	atomDen = atomDen * 2;
					break;
				case EXNAMBU: if( rtbd.Lat.parameter.STR("spin") == "on")	atomDen = atomDen / 2;
					break;
			}
			totalDen += atomDen;
			
			if( allSubIndex.size() > 0){
				rtbd.order(atomI.atomName+" den") = xvec(atomDen);
			}
		}
		
		
		return totalDen;
	}
	void	calculateChemicalPotential(bool printResult = false)							{

		
		KHamEvd(stbd, false);
		double	destDen = countTotalElectron();
		
		double minE = stbd.minE;
		double maxE = stbd.maxE;
		
		double		dest_den_diff = 0.01;
		double		destMu = 0.5*(minE + maxE);
		unsigned	printLen = 10;
		double		n_diff = 1;
		
		
		if( printResult )
		cout<<" "<<fmt("Mu",printLen)
			<<" "<<fmt("Dest e-den",printLen)
			<<"   "<<fmt("True e-den",printLen)
			<<"   "<<fmt("total_n_diff",printLen)<<endl;
		
		do{
			double totalDen = calculateElectronFilling(stbd, destMu);
			n_diff  = destDen - totalDen;
			
			if( printResult )
			cout<<" "<<fmt(destMu,printLen)
				<<" "<<fmt(destDen,printLen)
				<<" - "<<fmt( destDen-n_diff,printLen )
				<<" = "<<fmt(n_diff,printLen) <<endl;
			
			if (abs(n_diff) > abs(dest_den_diff) ) {
				if (n_diff > 0) { minE = destMu; }
				if (n_diff < 0) { maxE = destMu; }
				destMu = (maxE+minE)/2;
			}
		
		}while ( abs(n_diff) > abs(dest_den_diff) );
		
		stbd.order.save("");
		
		Lat.parameter.VAR("Mu") = destMu;
		spinNormalLat.parameter.VAR("Mu") = destMu;
	}
	
public:
	/* Expand the lattice into a larger size. */
	void saveExpandedLattice(unsigned N1, unsigned N2, unsigned N3)							{
		auto	parseFilename = split(Lat.FileName(),".");
		parseFilename.pop_back();
		
		string	filename_expand = "";
		for( auto & parts : parseFilename){ filename_expand += parts+"."; }
		
		filename_expand += IntToStr(N1)+"x"+IntToStr(N2)+"x"+IntToStr(N3)+".lat";
		cout<<endl<<">>>'"<<filename_expand<<"'"<<endl;
		
		vector<Atom>	expandedAtomList;
		auto atomList	= Lat.getAtomList();
		
		vector<rmap>	expandedOptList;
		OrderParameter	tmpOrder(Lat);
		tmpOrder.load();
		auto optList = tmpOrder.getOptList();
		
		// Expand the lattice in here.
		auto AVec = Lat.basisVector.getAVec();
		if		( AVec.size() == 3){
			for( unsigned i = 0 ; i<N1 ; i+=1)
			for( unsigned j = 0 ; j<N2 ; j+=1)
			for( unsigned k = 0 ; k<N3 ; k+=1){
				for(unsigned index=0 ; index<atomList.size() ; index++){
					Atom tmpAtom;
					tmpAtom = atomList[index];
					tmpAtom.pos = tmpAtom.pos+ i*AVec[0] +j*AVec[1] +k*AVec[2];
					expandedAtomList.push_back(tmpAtom);
					expandedOptList.push_back(optList[index]);
				}
			}
		}
		else if	( AVec.size() == 2){
			for( unsigned i = 0 ; i<N1 ; i+=1)
			for( unsigned j = 0 ; j<N2 ; j+=1){
				for(unsigned index=0 ; index<atomList.size() ; index++){
					Atom tmpAtom;
					tmpAtom = atomList[index];
					tmpAtom.pos = tmpAtom.pos+ i*AVec[0] +j*AVec[1];
					expandedAtomList.push_back(tmpAtom);
					expandedOptList.push_back(optList[index]);
				}
			}
		}
		else if	( AVec.size() == 1){
			for( unsigned i = 0 ; i<N1 ; i+=1){
				for(unsigned index=0 ; index<atomList.size() ; index++){
					Atom tmpAtom;
					tmpAtom = atomList[index];
					tmpAtom.pos = tmpAtom.pos+ i*AVec[0];
					expandedAtomList.push_back(tmpAtom);
					expandedOptList.push_back(optList[index]);
				}
			}
		}
		
		// Save the expanded file.
		ofstream outfile(filename_expand);
		outfile<<Lat.parameter.getFileString();
		outfile<<Lat.kSymmPointParser.getFileString();
		outfile<<Lat.bondVector.getFileString();
		outfile<<Lat.basisVector.getFileString(N1,N2,N3);
		outfile<<Lat.orbitalProfile.getFileString();
		outfile<<Lat.atomParser()<<endl;
		for(unsigned ii=0 ; ii<expandedAtomList.size() ; ii++){
			outfile<<expandedAtomList[ii].posToStr()<<endl;
		}
		outfile.close();
		
		// Expand the order parameter.
		Lattice expLat(filename_expand);
		OrderParameter expOrder(expLat);
		expOrder.setOptList(expandedOptList);
		expOrder.save();
		
		
        ifstream infile(Lat.FileName()+".tbm");
		ofstream tbmOutfile(filename_expand+".tbm");
		if ( infile.is_open() ) {
			string line = "";
			while(getline(infile, line)){
				tbmOutfile<<line<<endl;
			}
		}
		tbmOutfile.close();
		infile.close();
	}
	void render(){
		run();
	}
};























