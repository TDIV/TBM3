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
		
		rtbd.order.save();
		
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
			<<" "<<fmt("True e-den",printLen)
			<<" "<<fmt("total_n_diff",printLen)<<endl;
		
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
		
		Lat.parameter.VAR("Mu") = destMu;
		spinNormalLat.parameter.VAR("Mu") = destMu;
	}
	double	iterateDenOrder			(OrderParameter & newOrder, double mix = 0.1)			{
		
		if( Lat.parameter.VAR("disable_quantum", 1).real() == 1 ){ return 0; }
		
		tbd.order_old = tbd.order;
		calculateChemicalPotential(false);
		
		double ratio_a = 0;
		double ratio_b = 1;
		
		if( mix > 0 and mix < 1){
			ratio_a = 1-mix;
			ratio_b = mix;
		}
		
		double max_den_diff = 0;
		while( iterate() ){
			auto parameter_old = tbd.order_old.findOrder(Lat.getAtom(),	"@:den");
			auto parameter_new = stbd.order.findOrder(Lat.getAtom(),"@:den");
			if( parameter_old.first and parameter_new.first ){
				auto den_diff = abs(parameter_old.second[0].real() - parameter_new.second[0].real());
				if( max_den_diff < den_diff ) max_den_diff = den_diff;
			}
			
			auto mixOrder = ratio_a * parameter_old.second + ratio_b * parameter_new.second;
			
			newOrder.set(Lat.getAtom().atomIndex, "@:den", mixOrder);
		}
		tbd.calculateEnergy();
		newOrder.save();
		
		return max_den_diff;
	}

public:
	/* Execute the calculations. */
	void render				()																{
		run();
	}
	
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

	/* Convert the lattice file formate to VESTA formate. */
	void convertTo_VESTA	(vector<string> argsForOrder)									{
		string filename = Lat.FileName()+".vesta";
		ofstream outfile(filename);
		outfile<<"#VESTA_FORMAT_VERSION 3.3.0"<<endl<<endl<<endl;
		outfile<<"CRYSTAL"<<endl<<endl;
		outfile<<"TITLE"<<endl;
		outfile<<Lat.FileName()<<endl<<endl;
		
		// Convert the unitcell structure.
		outfile<<"CELLP"<<endl;
		auto AVec = Lat.basisVector.getAVec();
		double len_a = sqrt(cdot(AVec[0], AVec[0]));
		double len_b = sqrt(cdot(AVec[1], AVec[1]));
		double len_c = sqrt(cdot(AVec[2], AVec[2]));
		double ang_bc= acos( cdot(AVec[1], AVec[2])/(len_b*len_c) )*360/(2*pi);
		double ang_ca= acos( cdot(AVec[2], AVec[0])/(len_c*len_a) )*360/(2*pi);
		double ang_ab= acos( cdot(AVec[0], AVec[1])/(len_a*len_b) )*360/(2*pi);
		outfile<<" "<<fformat(len_a, 10);
		outfile<<" " <<fformat(len_b,10);
		outfile<<" " <<fformat(len_c, 9);
		outfile<<" " <<fformat(abs(ang_bc), 10);
		outfile<<" " <<fformat(abs(ang_ca), 10);
		outfile<<" " <<fformat(abs(ang_ab), 10);
		outfile<<endl<<"  0.000000   0.000000   0.000000   0.000000   0.000000   0.000000"<<endl;
		
		// Convert the atom position.
		outfile<<"STRUC"<<endl;
		auto atomList = Lat.getAtomList();
		for(auto & atom: atomList){
			outfile	<<"  "<<fformat(atom.atomIndex+1,5)
					<<" "<<fformat(atom.atomName,6)
					<<" "<<fformat(atom.orbitalOriginIndex,2)
					<<"  1.0000  ";
			double pos_a = cdot(atom.pos, AVec[0])/(len_a*len_a);
			double pos_b = cdot(atom.pos, AVec[1])/(len_b*len_b);
			double pos_c = cdot(atom.pos, AVec[2])/(len_c*len_c);
			outfile<<fformat(pos_a, 11);
			outfile<<fformat(pos_b, 11);
			outfile<<fformat(pos_c, 11);
			outfile<<"  1a       1";
			outfile<<endl;
			outfile<<"                            0.000000   0.000000   0.000000  0.00"<<endl;
		}
		outfile<<"  0 0 0 0 0 0 0"<<endl;
		
		// Handle the arguments for plotting the order parameters.
		// vec=@:cspin>1,2,3
		// vec=@:1:4den>2,3,4
		vector<boost::tuple<string, Atom, x_mat> > orderAtomVecRColor;
		
		for( unsigned i = 2 ; i<argsForOrder.size() ; i++ ){
			string orderArg = argsForOrder[i];
			cout<<" >> Handling order argument: ";
			cout<<orderArg<<"."<<endl;
			auto argParser = split(orderArg, "=");
			if( argParser.size() <= 1 )		{
				cout<<"Warning, \'"<<orderArg<<"\' is not an effective order argument. Ignored!"<<endl;
				continue;
			}
			
			vestaOrderHandler(argParser[0], argParser[1], orderAtomVecRColor);

		}
		
		outfile<<"VECTR"<<endl;
		for( unsigned ii=0 ; ii < orderAtomVecRColor.size() ; ii++ ){
			
			auto & ordElement = orderAtomVecRColor[ii];
			auto arg = ordElement.get<0>();
			auto vec = ordElement.get<2>();
			auto atom = ordElement.get<1>();
			
			if( arg == "vec"){
				outfile<<"   "<<fformat(ii+1,4);
				for( unsigned i=0; i<vec.size() ; i++){
					outfile<<fformat(vec[i].real());
				}
				outfile<<endl;
				outfile<<"    "<<fformat(atom.atomIndex+1,3)<<" 0    0    0    0"<<endl;
				outfile<<" 0 0 0 0 0"<<endl;
			}
		}
		outfile<<" 0 0 0 0 0"<<endl;
		
		outfile<<"VECTT"<<endl;
		for( unsigned ii=0 ; ii < orderAtomVecRColor.size() ; ii++ ){
			
			auto & ordElement = orderAtomVecRColor[ii];
			auto arg = ordElement.get<0>();
			auto vec = ordElement.get<2>();
			
			if( arg == "vec"){
				outfile<<"   "<<fformat(ii+1,4)<<"0.200 255   0   0 0"<<endl;
			}

		}
		outfile<<" 0 0 0 0 0"<<endl;
		
		cout<<endl;
		                              
		outfile.close();
	}
private:
	void vestaOrderHandler( string ordArg, string ordNameArg, vector<boost::tuple<string, Atom, x_mat> > & orderAtomVec){
		
		auto ordNameParser = split(ordNameArg, "~");
		string ordName = ordNameParser[0];
		vector<int> ordVecIndex;
		
		if( ordNameParser.size() == 1 ){
			ordVecIndex.push_back(1);
			ordVecIndex.push_back(2);
			ordVecIndex.push_back(3);
		}
		else if( ordNameParser.size() == 2 ){
			for( auto indexStr: split(ordNameParser[1], ",")){
				ordVecIndex.push_back(StrToInt(indexStr));
			}
		}
		
		tbd.order.load();
		while(Lat.iterate()){
			auto ordMat = tbd.order.findOrder(Lat.getAtom(), ordName);
			if( ordMat.first ){
				x_mat tmpOrd(1, ordVecIndex.size());
				for(unsigned i=0 ; i<ordVecIndex.size() ; i++){
					unsigned index = ordVecIndex[i];
					if( index <= ordMat.second.size() )
						tmpOrd[i] = ordMat.second[index-1];
				}
				//cout<<ordMat.second<<" "<<tmpOrd<<endl;
				orderAtomVec.push_back(boost::make_tuple(ordArg, Lat.getAtom(), tmpOrd));
			}
		}
	
	}
};























