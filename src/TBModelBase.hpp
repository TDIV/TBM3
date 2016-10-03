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
//  TBModelBase.hpp
//  TBM^3
//
//  Created by Yuan-Yen Tai on 7/19/16.
//

// -----------------------------------
// Make a linear points between two different k-points.
// -----------------------------------
vector<r_mat> make_line(r_mat point_a, r_mat point_b, unsigned Nsteps=10){
	vector<r_mat> ll;
	for (double t=0; t<1; t+=1.0/Nsteps) {
		r_mat position= point_b-point_a;
		position = position*t;
		position = position+point_a;
		ll.push_back(position);
	}
	
	return ll;
}

class TBModelBase{
	
public:
	TBModelBase(string filename):
		Lat				(filename),
		spinNormalLat	(filename),
		TBD				(Lat, tbm),
		spinNormalTBD	(spinNormalLat, tbm),
		tbd				(TBD),
		stbd			(spinNormalTBD	)
	{
		kMeshSelection = "Nb";
	}
	
protected:
	Lattice		Lat;
	Lattice		spinNormalLat;
	TBMParser	tbm;
	void		loadTBM(bool isWithHamiltonianBlock = true){
		tbm.open(Lat.FileName()+".tbm", Lat.FileName(), isWithHamiltonianBlock);
		
		Lat.parameter	= tbm.parameter;
		Lat.bondVector	= tbm.bondVector;
		
		spinNormalLat.parameter	= tbm.parameter;
		spinNormalLat.bondVector= tbm.bondVector;
	}
	
private:
	TBDataSource TBD;
	TBDataSource spinNormalTBD;
	
protected:
	TBDataSource & tbd;		// A reference to object -> TBD
	TBDataSource & stbd;	// A reference to object -> spinNormalTBD
	
	bool iterate(){ return Lat.iterate(); }
	
	virtual void initOrder()	= 0;
	virtual void Hamiltonian()	= 0;
	virtual void run()			= 0;
	
	void		constructHam	(TBDataSource & rtbd, bool withMu = true)					{
		// Construct the Hamiltonian from the system-build structure: TBDataSource (tbd).
		
		if( Lat.parameter.VAR("disable_quantum", 0).real() == 0 )
			rtbd.constructTBMHam(withMu);
		
		// Construct the customized Hamiltonian.
		Hamiltonian();
	}

	/*-----------------------------------------------------
	 The Band structure calculation.
	 ------------------------------------------------------*/
	void	calculateBandStructure	(TBDataSource & rtbd, unsigned Nsteps=50)				{
		if( Lat.parameter.VAR("disable_quantum", 0).real() != 0 ){
			cout<< "Warning, due to flag 'disable_quantum' turned on."
				<< "The band structure calculation will be ignored."<<endl;
			return;
		}
		
		auto B = Lat.basisVector.getBVec();
		auto & b1 = B[0];
		auto & b2 = B[1];
		auto & b3 = B[2];
		
		vector<pair<string,r_mat> >	kSpaceHighSymmetryPoints;
		for(auto & kp: tbm.kSymmPointParser.kSymmPointList){
			r_mat B1 = kp.second[0]*b1;
			r_mat B2 = kp.second[1]*b2;
			r_mat B3 = kp.second[2]*b3;
			r_mat kPoint = B1+B2+B3;
			
			if ( kPoint.cols()==1 and kPoint.rows()==3) {
				kSpaceHighSymmetryPoints.push_back(make_pair(kp.first, kPoint));
			}
			
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
				
				if ( j== kSpaceHighSymmetryPoints.size() -1 ){
					string label = kSpaceHighSymmetryPoints[j].first;
					highSymmetryLine[highSymmetryLine.size()-1] =
						make_pair(label, kSpaceHighSymmetryPoints[j].second);
				}
			}
		}
		
		/* Store the band structure into 'xxx.lat.ban' */
		string	filename = Lat.FileName()+".ban";
		ofstream out(filename);
		
		constructHam(rtbd);
		for ( auto kp: highSymmetryLine){
			auto evv = rtbd.HamEvd(kp.second);
			evv.eigenVector.setPrintLength(16);
			out<< kp.first<<" "<< kp.second <<" "<<evv.eigenValue<<endl;
		}
		
		out.close();
	}
	void	calculateKWannierCenter	(TBDataSource & rtbd, unsigned Nsteps=50)				{
		if( Lat.parameter.VAR("disable_quantum", 0).real() != 0 ){
			cout<< "Warning, due to flag 'disable_quantum' turned on."
				<< "The Wannier center calculation will be ignored."<<endl;
			return;
		}
		
		//cout<<Nsteps<<endl;
		
		auto B = Lat.basisVector.getBVec();
		auto & b1 = B[0];
		auto & b2 = B[1];
		auto & b3 = B[2];
		
		// Translate the K-Points in the Bravis Vector.
		vector<boost::tuple<string, r_mat, unsigned, r_mat, r_mat, unsigned> > kWannierPointList;
		for(auto & unpacker: tbm.kWannierParser.kWannierPointList){
			auto label				= unpacker.get<0>();
			auto nKP				= unpacker.get<1>();
			
			r_mat B1 = nKP[0]*b1;
			r_mat B2 = nKP[1]*b2;
			r_mat B3 = nKP[2]*b3;
			r_mat kPoint = B1+B2+B3;
			
			auto nBand				= unpacker.get<2>();
			
			auto integralFrom		= unpacker.get<3>();
			r_mat F1 = integralFrom[0]*b1;
			r_mat F2 = integralFrom[1]*b2;
			r_mat F3 = integralFrom[2]*b3;
			r_mat FromPoint = F1+F2+F3;
			
			auto integralTo			= unpacker.get<4>();
			r_mat T1 = integralTo[0]*b1;
			r_mat T2 = integralTo[1]*b2;
			r_mat T3 = integralTo[2]*b3;
			r_mat ToPoint = T1+T2+T3;
			
			auto integralSequence	= unpacker.get<5>();
			
			if ( kPoint.size() == 3 and FromPoint.size() == 3 and ToPoint.size() == 3) {
				kWannierPointList.push_back(boost::make_tuple(label, kPoint, nBand, FromPoint, ToPoint, integralSequence));
			}
			
		}
		
		// Making high symmetry lines throughout high-symmetry points.
		vector<boost::tuple<string, r_mat, unsigned, r_mat, r_mat, unsigned> > kWannierHighSymmetryLine;
		if (kWannierPointList.size()>=2) {
			
			for (unsigned i=0; i<kWannierPointList.size()-1; i++) {
				unsigned j=i+1;
				
				auto I_label			= kWannierPointList[i].get<0>();
				auto I_KP				= kWannierPointList[i].get<1>();
				auto I_nBand			= kWannierPointList[i].get<2>();
				auto I_From				= kWannierPointList[i].get<3>();
				auto I_To				= kWannierPointList[i].get<4>();
				auto I_Sequence			= kWannierPointList[i].get<5>();
				
				auto J_label			= kWannierPointList[j].get<0>();
				auto J_KP				= kWannierPointList[j].get<1>();
				auto J_Sequence			= kWannierPointList[j].get<5>();
				
				if( I_Sequence == J_Sequence ){
				
					auto	partOfHighSymmetryLine = make_line(I_KP, J_KP, Nsteps);
					
					for (unsigned ii=0; ii<partOfHighSymmetryLine.size(); ii++) {
						string label = ii==0 ? I_label : "-";
						
						kWannierHighSymmetryLine.push_back(
							boost::make_tuple(label, partOfHighSymmetryLine[ii], I_nBand, I_From, I_To, I_Sequence)
						);
					}
				
					if ( j== kWannierPointList.size() -1 ){
						kWannierHighSymmetryLine[kWannierHighSymmetryLine.size()-1] =
							boost::make_tuple(J_label, J_KP, I_nBand, I_From, I_To, I_Sequence);
					}
				}
				else{
					kWannierHighSymmetryLine[kWannierHighSymmetryLine.size()-1] =
						boost::make_tuple(I_label, I_KP, I_nBand, I_From, I_To, I_Sequence);
				}
			}
		}
		
		string	filename = Lat.FileName()+".kwcenter";
		ofstream out(filename);
		vector<boost::tuple<string, r_mat, r_mat> > kWannierCenter;
		for( auto & unpacker: kWannierHighSymmetryLine ){
			auto &	label	= unpacker.get<0>();
			auto &	kpoint	= unpacker.get<1>();
			auto &	nBand	= unpacker.get<2>();
			auto &	IntFrom	= unpacker.get<3>();
			auto &	IntTo	= unpacker.get<4>();
			auto &	Sequence= unpacker.get<5>();
			
			auto integralLine = make_line(IntFrom, IntTo, Nsteps);
			
			x_mat Dk(nBand,nBand), F(nBand,nBand);
			for( unsigned i=0 ; i<Dk.cols() ; i++){
				Dk(i,i) = 1;
			}
			
			x_mat vec0, vecI, vecJ;
			auto evvI = rtbd.HamEvd(kpoint);
			vec0 = evvI.eigenVector;
			vecI = vec0;
			
			for( unsigned j=1; j<integralLine.size() ; j++ ){
				
				auto kJ = kpoint + integralLine[j];
				auto evvJ = rtbd.HamEvd(kJ);
				
				vecJ = evvJ.eigenVector;
				if( j==integralLine.size()-1 ){ vecJ = vec0; }
				
				for( unsigned m=0 ; m<nBand ; m++)
				for( unsigned n=0 ; n<nBand ; n++){
					F(m,n) = 0;
					for( unsigned ii=0 ; ii<vecI.cols() ; ii++){
						F(m,n) += conj( vecI(ii,m) ) *  vecJ(ii,n)  ;
					}
				}
				
				Dk = Dk * F;
				vecI = vecJ;
			}
			x_mat Eig,Vec;
			Dk.gevd(Eig, Vec);
			
			out<<label<<" "<<kpoint<<" ";
			
			vector<double> thetaList;
			for( unsigned i = 0 ; i<Eig.size() ; i++){
				thetaList.push_back(complexToArg(Eig[i].real(), Eig[i].imag()));
			}
			sort(thetaList.begin(), thetaList.end());
			
			out<<"[[ ";
			for( auto & tt: thetaList){
				out<<fformat(tt)<<" ";
			}
			out<<" ]]";
			out<<endl;
		}
		
		out.close();
		
	}
	
	/*-----------------------------------------------------
	 Diagonalize the full Hamiltonian and store the results in KEigenValVec.
	 ------------------------------------------------------*/
	string	kMeshSelection;
	void	KHamEvd					(TBDataSource & rtbd, bool withMu = true)				{
		
		constructHam(rtbd, withMu);
		
		if( Lat.parameter.VAR("disable_quantum", 0).real() != 0 ) return;
		
		auto Nb = Lat.parameter.VEC(kMeshSelection, Lat.parameter.VEC("Nb"));
		auto B = Lat.basisVector.getBVec();
		
		rtbd.KEigenValVec.clear();
		
		//if( Nb.size() == 1){
		//	auto & b1 = B[0];
		//	auto N1 = Nb[0].real();
		//	
		//	for (double i1=0; i1<N1; i1++){
		//		auto kPoint = (i1/N1)*b1;
		//		auto  tmpEVV = rtbd.HamEvd(kPoint);
		//		if( tmpEVV.message == "Success")
		//			rtbd.KEigenValVec.push_back(tmpEVV);
		//	}
		//}
		//if( Nb.size() == 2){
		//	auto & b1 = B[0];
		//	auto & b2 = B[1];
		//	auto N1 = Nb[0].real(), N2 = Nb[1].real();
		//	for (double i1=0; i1<N1; i1++)
		//	for (double i2=0; i2<N2; i2++){
		//		auto kPoint = (i1/N1)*b1 + (i2/N2)*b2;
		//		auto  tmpEVV = rtbd.HamEvd(kPoint);
		//		if( tmpEVV.message == "Success")
		//			rtbd.KEigenValVec.push_back(tmpEVV);
		//	}
		//}
		if( Nb.size() == 3){
			auto & b1 = B[0];
			auto & b2 = B[1];
			auto & b3 = B[2];
			auto N1 = Nb[0].real(), N2 = Nb[1].real(), N3 = Nb[2].real();
			for (double i1=0; i1<N1; i1++)
			for (double i2=0; i2<N2; i2++)
			for (double i3=0; i3<N3; i3++) {
				r_mat B1 = (i1/N1)*b1;
				r_mat B2 = (i1/N1)*b2;
				r_mat B3 = (i1/N1)*b3;
				auto kPoint = B1+B2+B3;
				auto  tmpEVV = rtbd.HamEvd(kPoint);
				if( tmpEVV.message == "Success")
					rtbd.KEigenValVec.push_back(tmpEVV);
			}
		}
		
		
	}
	
	/*-----------------------------------------------------
	 Calculate total electron
	 ------------------------------------------------------*/
	double	countTotalElectron			()													{
		double totalElectron = 0;
		while(iterate()){
			auto atomI = Lat.getAtom();
			totalElectron += tbm.coreCharge.getCharge(atomI.atomName);
		}
		return totalElectron;
	}
	r_var	calculateElectronFilling	(TBDataSource & rtbd, double Mu=0)					{
		
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
				case NAMBU:
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
	void	calculateChemicalPotential	(bool printResult = false)							{
		
		stbd.order  = tbd.order;

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
	void	calculateLDOS				(TBDataSource & rtbd)								{
		if( Lat.parameter.VAR("disable_quantum", 0).real() != 0 ){
			cout<< "Warning, due to flag 'disable_quantum' turned on."
				<< "The LDOS calculation will be ignored."<<endl;
			return;
		}
		
		double stepE = abs(Lat.parameter.VAR("ldos_dE", 0.01).real());
		double Gamma = abs(Lat.parameter.VAR("ldos_Gamma", 0.05).real());

		// Define the LDOS calculation content (atomLDOSList).
		typedef vector< boost::tuple<string, string, unsigned> > OrbitalIndexLabel;
		typedef vector< pair<double, double> > LDOS;
		vector<boost::tuple<Atom, OrbitalIndexLabel, LDOS > > atomLDOSList;
		
		
		kMeshSelection = "ldos_Nb";
		KHamEvd(rtbd);
		kMeshSelection = "Nb";
		
		x_mat ldos_E(1,2);
		ldos_E[0] = rtbd.minE;
		ldos_E[1] = rtbd.maxE;
		ldos_E = Lat.parameter.VEC("ldos_E", ldos_E);
		
		if( ldos_E[0].real() > ldos_E[1].real()){
			auto tmpE = ldos_E[1];
			ldos_E[1] = ldos_E[0];
			ldos_E[0] = tmpE;
		}
		
		cout<<ldos_E<<endl;
		
		// Handling the LDOS label for the up coming calculation.
		for( auto & elem: tbm.ldosList.LDOSSelector ){
			
			auto isAtom = Lat.getAtom(elem.first);
			auto & atom = isAtom.second;
			if( isAtom.first ){
				OrbitalIndexLabel	orbitalIndexLabel;
				LDOS				ldos;
				
				if( elem.second.size() == 0 ){ // If size = 0, this will sum over all index for that site.
					auto & indexList = atom.allIndexList();
					for( auto & ilelem: indexList){
						auto e = split(ilelem.first,".");
						auto orbitalNum = atom.getOrbitalNumber(e[0]);
						orbitalIndexLabel.push_back(boost::make_tuple(orbitalNum, e[1], ilelem.second));
					}
				}
				else {
				
					for(auto & refStr: elem.second){
						
						if( IsIntStr(refStr) ){ // Input format of : 1 2 3
							auto indexList = isAtom.second.orbitalIndexList(refStr);
							for( auto & ilelem: indexList){
								orbitalIndexLabel.push_back(boost::make_tuple(refStr, ilelem.first, ilelem.second));
							}
						}
						else{		// Input format of : 1u
							auto indexList = isAtom.second.spinIndexList(refStr);
							char spin = refStr[refStr.size()-1];
							refStr.pop_back();
							for( auto & ilelem: indexList){
								string ABLabel = ilelem.first;
								if( ABLabel[ABLabel.size()-1] != spin) ABLabel.push_back(spin);
								orbitalIndexLabel.push_back(boost::make_tuple(refStr, ABLabel, ilelem.second));
							}
						}
					}
				}
				
				for( double eng = ldos_E[0].real(); eng <= ldos_E[1].real()+4*stepE ; eng += stepE ){
					ldos.push_back(make_pair(eng, 0));
				}
				
				atomLDOSList.push_back(boost::make_tuple(isAtom.second, orbitalIndexLabel, ldos));
			}
		}
		
		/* Store the band structure into 'xxx.lat.ldos' */
		string	filename = Lat.FileName()+".ldos";
		ofstream out(filename);
		
		out<<"import numpy as np"<<endl;
		out<<"import matplotlib.pyplot as plt"<<endl<<endl;
		
		unsigned outputIndex = 0;
		
		for( auto & elem: atomLDOSList ){
			
			out<<"# "<<elem.get<0>().atomName<<" "<<elem.get<0>().pos<<" >>> ";
			for( auto & in: elem.get<1>()){
				
				out<<" "<<in.get<0>()<<"."<<in.get<1>();
				
				int sign = in.get<1>()[0] == 'B' ? -1 :  1;
				
				for( auto & ldosIter: elem.get<2>() ){
					double w = ldosIter.first;
					for (unsigned k=0; k<rtbd.KEigenValVec.size(); k++) {
						//auto & kpoint = rtbd.KEigenValVec[k].kPoint;
						auto & kEig = rtbd.KEigenValVec[k].eigenValue;
						auto & kVec = rtbd.KEigenValVec[k].eigenVector;
						
						for ( int n=0; n<kEig.size(); n++) {
							double e_w		= kEig[n]-sign * w;
							double delta	= Gamma/(pi*(e_w*e_w + Gamma*Gamma));
							
							x_var u_n		= kVec(in.get<2>(), n);
							ldosIter.second += (conj(u_n)*u_n).real() * delta/tbd.KEigenValVec.size();
						}
					}
				}
			}
			
			out<<endl;
			string atomIndexStr = IntToStr( elem.get<0>().atomIndex );
			
			out<<"x"<<outputIndex<<"=np.array([";
			for( auto & ldosIter: elem.get<2>() ){
				out<<fformat(ldosIter.first, 12)<<",";
			}
			out<<"])";
			out<<endl;
			
			out<<"y"<<outputIndex<<"=np.array([";
			for( auto & ldosIter: elem.get<2>() ){
				out<<fformat(ldosIter.second, 12)<<",";
			}
			out<<"])"<<endl;
			out<<endl;
			outputIndex++;
		}
		
		outputIndex = 0;
		for( auto & elem: atomLDOSList ){
			out<<"plt.plot(x"<<outputIndex<<", y"<<outputIndex<<", '-', linewidth=2)"<<endl;
			outputIndex++;
		}
		out<<"plt.show()"<<endl;
		

		out.close();
	}
	double	iterateDenOrder				(OrderParameter & newOrder, double mix = 0.1)		{
		
		if( Lat.parameter.VAR("disable_quantum", 1).real() == 1 ){ return 0; }
		
		tbd.order_old = tbd.order;
		// *********************************
		calculateChemicalPotential(false);
		// *********************************
		
		if( tbd.Lat.parameter.STR("spin") == "on" ) {
			tbd.calculate4DensityOrder();
		}
		
		double ratio_a = 0;
		double ratio_b = 1;
		
		if( mix > 0 and mix < 1){
			ratio_a = 1-mix;
			ratio_b = mix;
		}
		
		double max_den_diff = 0;
		
		// Mixing and difference checking for @:den
		while( iterate() ){
			
			auto parameter_old = tbd.order_old.findOrder(Lat.getAtom(),	"@:den");
			auto parameter_new = stbd.order.findOrder(Lat.getAtom(),"@:den");
			
			if( parameter_old.first and parameter_new.first ){
				auto den_diff = abs(parameter_old.second[0].real() - parameter_new.second[0].real());
				if( max_den_diff < den_diff ) max_den_diff = den_diff;
				
				// Set up updated @:den order parameter
				auto para_old = ratio_a * parameter_old.second;
				auto para_new = ratio_b * parameter_new.second;
				auto mixOrder = para_old + para_new;
				newOrder.set(Lat.getAtom().atomIndex, "@:den", mixOrder);
			}
			
			// Set up updated @:coulomb order parameter
			auto parameter_coulomb = stbd.order.findOrder(Lat.getAtom(),	"@:coulomb");
			if( parameter_coulomb.first ){
				newOrder.setNew(Lat.getAtom().atomIndex, "@:coulomb", parameter_coulomb.second);
			}
		}
		
		// Mixing and difference checking for @:n:4den
		if( tbd.Lat.parameter.STR("spin") == "on" )
		while( iterate() ){
			auto atom = Lat.getAtom();
			auto orbitals = atom.allOrbitalLabel();
			
			for( auto orb: orbitals ){
				string orbN = atom.getOrbitalNumber(orb);
				string orderKey = "@:"+orbN+":4den";
				
				auto parameter_old = tbd.order_old.findOrder(Lat.getAtom(),	orderKey);
				auto parameter_new = tbd.order.findOrder(Lat.getAtom(),orderKey);
			
				if( parameter_old.first and parameter_new.first ){
					auto vdiff = parameter_old.second - parameter_new.second;
					auto den_diff0 = abs(vdiff[0]);
					auto den_diff1 = abs(vdiff[1]);
					auto den_diff2 = abs(vdiff[2]);
					auto den_diff3 = abs(vdiff[3]);
					if( max_den_diff < den_diff0 ) max_den_diff = den_diff0;
					if( max_den_diff < den_diff1 ) max_den_diff = den_diff1;
					if( max_den_diff < den_diff2 ) max_den_diff = den_diff2;
					if( max_den_diff < den_diff3 ) max_den_diff = den_diff3;
					
					// Set up updated @:n:4den order parameter
					auto para_old = ratio_a * parameter_old.second;
					auto para_new = ratio_b * parameter_new.second;
					auto mixOrder = para_old + para_new;
					newOrder.setNew(Lat.getAtom().atomIndex, orderKey, mixOrder);
				}
			}
		}
		
		tbd.calculateEnergy();
		newOrder.save();
		
		return max_den_diff;
	}
	
	void	calculateSpinSusceptibility	(TBDataSource & rtbd);

public:
	/* Execute the calculations. */
	void render				()																{

		loadTBM();
		Lat.createAtomList( tbm.parameter.STR("spin", "on"), tbm.parameter.STR("space","normal"));
		spinNormalLat.createAtomList();
		
		string solver = tbm.parameter.STR("SOLVER", "GPU");
		if		( solver == "CPU")		{ gmt::SOLVER = CPU; }
		else if	( solver == "GPU")		{ gmt::SOLVER = GPU; }
		else if	( solver == "GPUx2")	{ gmt::SOLVER = GPUx2; }
		else if	( solver == "GPUx3")	{ gmt::SOLVER = GPUx3; }
		else if	( solver == "GPUx4")	{ gmt::SOLVER = GPUx4; }
		else if	( solver == "GPUx5")	{ gmt::SOLVER = GPUx5; }
		else if	( solver == "GPUx6")	{ gmt::SOLVER = GPUx6; }
		else if	( solver == "GPUx7")	{ gmt::SOLVER = GPUx7; }
		else if	( solver == "GPUx8")	{ gmt::SOLVER = GPUx8; }
		else if	( solver == "GPUx9")	{ gmt::SOLVER = GPUx9; }
		else if	( solver == "GPUx10")	{ gmt::SOLVER = GPUx10;}
		
		run();
	}
	
	/* Expand the lattice into a larger size. */
	void saveExpandedLattice(unsigned N1, unsigned N2, unsigned N3)							{
		
		loadTBM(false);
		Lat.createAtomList();
		
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
					r_mat A0 = i*AVec[0];
					r_mat A1 = j*AVec[1];
					r_mat A2 = k*AVec[2];
					tmpAtom.pos = tmpAtom.pos+ A0 + A1 + A2;
					expandedAtomList.push_back(tmpAtom);
					expandedOptList.push_back(optList[index]);
				}
			}
		}
		else{
			ErrorMessage("Error, should be given a 3D basis vector.");
		}

		
		// Save the expanded file.
		ofstream outfile(filename_expand);
		outfile<<Lat.basisVector.getFileString(N1,N2,N3);
		outfile<<Lat.orbitalProfile.getFileString();
		outfile<<Lat.atomParser()<<endl;
		for(unsigned ii=0 ; ii<expandedAtomList.size() ; ii++){
			outfile<<expandedAtomList[ii].posToStr()<<endl;
		}
		outfile.close();
		
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
		
		// Expand the order parameter.
		Lattice expLat(filename_expand);
		expLat.createAtomList();
		
		OrderParameter expOrder(expLat);
		expOrder.setOptList(expandedOptList);
		expOrder.save();
	}

	/* Convert the lattice file formate to VESTA formate. */
	void convertTo_VESTA	(vector<string> argsForOrder)									{
		
		loadTBM(false);
		Lat.createAtomList("on", "normal", false);
		
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
		double V123 = cdot(curl(AVec[0], AVec[1]), AVec[2]);
		double V213 = cdot(curl(AVec[1], AVec[0]), AVec[2]);
		double V312 = cdot(curl(AVec[2], AVec[0]), AVec[1]);
		
		outfile<<"STRUC"<<endl;
		auto atomList = Lat.getAtomList();
		for(auto & atom: atomList){
			auto atomName = atom.atomName;
			
			replaceAll(atomName, "-", ".");
			auto atomNameParser = split(atomName, ".");
			
			outfile	<<"  "<<fformat(atom.atomIndex+1,5)
					<<" "<<fformat(atomNameParser[0],6)
					<<" "<<fformat(atom.atomName+"_"+IntToStr(atom.orbitalOriginIndex),8)
					<<"  1.0000  ";
			
			double pos_a = cdot(curl(atom.pos, AVec[1]), AVec[2])*(1.0/V123);
			double pos_b = cdot(curl(atom.pos, AVec[0]), AVec[2])*(1.0/V213);
			double pos_c = cdot(curl(atom.pos, AVec[0]), AVec[1])*(1.0/V312);
			outfile<<fformat(pos_a, 11)<<" ";
			outfile<<fformat(pos_b, 11)<<" ";
			outfile<<fformat(pos_c, 11)<<" ";
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
			
			if( arg == "vec" || arg == "rvec" || arg == "gvec" || arg == "bvec" || arg == "rgbvec"){
				outfile<<"   "<<fformat(ii+1,4);
				for( unsigned i=0; i<vec.size() ; i++){
					outfile<<fformat(vec[i].real())<<" ";
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
			
			if( arg == "vec" )	{ outfile<<"   "<<fformat(ii+1,4)<<" 0.200   0   0   0 0"<<endl; }
			if( arg == "rvec")	{ outfile<<"   "<<fformat(ii+1,4)<<" 0.200 255   0   0 0"<<endl; }
			if( arg == "gvec")	{ outfile<<"   "<<fformat(ii+1,4)<<" 0.200   0 255   0 0"<<endl; }
			if( arg == "bvec")	{ outfile<<"   "<<fformat(ii+1,4)<<" 0.200   0   0 255 0"<<endl; }
			if( arg == "rgbvec"){
				auto vecLen = sqrt(cdot(vec,vec).real());
				unsigned colR = (unsigned) abs(255* vec[0].real() / vecLen);
				unsigned colG = (unsigned) abs(255* vec[1].real() / vecLen);
				unsigned colB = (unsigned) abs(255* vec[2].real() / vecLen);
				outfile<<"   "<<fformat(ii+1,4)<<" 0.200 "<<colR<<" "<<colG<<" "<<colB<<" "<<" 0"<<endl;
			}

		}
		outfile<<" 0 0 0 0 0"<<endl;
		
		cout<<endl;
		                              
		outfile.close();
	}
	
	/* Shift the cordinate of each atom. */
	void shiftXYZ			(double X, double Y, double Z)									{
		
		loadTBM(false);
		Lat.createAtomList("on", "normal", false);
		
		tbd.order.load();
		Lat.shiftXYZ(X, Y, Z);
		
		// Save the expanded file.
		ofstream outfile(Lat.FileName());
		outfile<<Lat.basisVector.getFileString();
		outfile<<Lat.orbitalProfile.getFileString();
		outfile<<Lat.atomParser()<<endl;
		while(Lat.iterate()){
			outfile<<Lat.getAtom().posToStr()<<endl;
		}
		outfile.close();
		
		tbd.order.save();
	}
	
	/* Change the atom name. */
	void changeAtom(vector<string> optList)													{
		
		loadTBM(false);
		Lat.createAtomList("on", "normal", false);
		
		tbd.order.load();
		Lat.changeAtomName(optList);
		
		// Save the expanded file.
		ofstream outfile(Lat.FileName());
		outfile<<Lat.basisVector.getFileString();
		outfile<<Lat.orbitalProfile.getFileString();
		outfile<<Lat.atomParser()<<endl;
		while(Lat.iterate()){
			outfile<<Lat.getAtom().posToStr()<<endl;
		}
		outfile.close();
		
		tbd.order.save();
		
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


#include "TBModelBaseTmp.hpp"




















