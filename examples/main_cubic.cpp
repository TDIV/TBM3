
#include <stdlib.h>
#include "tbm3.h"

using namespace lift;

// Define the model
class MyModel: public TBModel, public HoppingModelConstructor{
public:
	MyModel(Lattice lat):	TBModel(lat), HoppingModelConstructor(lat){
		init_order();
		render();
	}

	void init_order()		{
		cout<<"Initialize the calculation ..."<<endl<<endl;
		
		Lat.bondStringMap["+x"] = "+1..";
		Lat.bondStringMap["-x"] = "-1..";
		Lat.bondStringMap["+y"] = ".+1.";
		Lat.bondStringMap["-y"] = ".-1.";
		Lat.bondStringMap["+z"] = "..+1";
		Lat.bondStringMap["-z"] = "..-1";
		
		while( pair_iterate() ) {
			auto pit = getPair();
			hoppingOrder(pit, "Fe:Fe:+1.. hop.1.1") = 1;
			hoppingOrder(pit, "Fe:Fe:-1.. hop.1.1") = 1;
			hoppingOrder(pit, "Fe:Fe:.+1. hop.1.1") = 1;
			hoppingOrder(pit, "Fe:Fe:.-1. hop.1.1") = 1;
			hoppingOrder(pit, "Fe:Fe:..+1 hop.1.1") = 1;
			hoppingOrder(pit, "Fe:Fe:..-1 hop.1.1") = 1;
		}
		
		hoppingOrder.save("hop");
		
		initHoppingTerms();

		// Initialize the AFM order
		while (site_iterate()) {
			auto si = getSite();
			int rx=-1, ry=-1, rz=-1;
			if (int(si.pos[0])%2 ) rx=1;
			if (int(si.pos[1])%2 ) ry=1;
			if (int(si.pos[2])%2 ) rz=1;
			
			int mag_sign = rx*ry*rz;
			
			double theta=0, phi=0;
			if (mag_sign > 0) {
				theta	= 0;
				phi		= 0;
			} else if (mag_sign < 0) {
				theta	= pi;
				phi		= 0;
			}
			
			double	Sx = sin(theta)*cos(phi);
			double	Sy = sin(theta)*sin(phi);
			double	Sz = cos(theta);
			
			initOrder(si, "Fe Sx")	=Sx;
			initOrder(si, "Fe Sy")	=Sy;
			initOrder(si, "Fe Sz")	=Sz;
		}
		initOrder.save("");
		
		/*
		OrderParameter pairOrder(Lat);
		while (pair_iterate()) {
			auto pit = getPair();
			pairOrder(pit, "Fe:Fe:+1.. hop.1u.1u") = 2;
			pairOrder(pit, "Fe:Fe:.+1. hop.1u.1u") = 2;
		}
		pairOrder.save("pair");
		*/
	}
	
	void Hamiltonian()		{
		add_Chemical_Potential();
		
		OrderParameter order = initOrder;
		double Jh = -VAR("Jh").real();

		//---site iteration---
		while (site_iterate()) {
			auto si = getSite();
			add_hund_spin("Fe 1", Jh, order.getVars(si, "Sx Sy Sz"));
		}

		constructHoppingHamiltonian(Ham, k_space);
		
		//---pair iteration---
		//while (pair_iterate()) {
		//	auto pit = getPair();
		//	cout<<pit.AtomI.atomIndex()<<" "<<pit.AtomJ.atomIndex()<<endl;
		//	/*
		//	auto si = pit.AtomI;
		//	auto sj = pit.AtomJ;
		//	*/
		//	
		//	add_bond( "Fe:1u  Fe:1u  +x", VAR("t") );
		//	add_bond( "Fe:1u  Fe:1u  +y", VAR("t") );
		//	add_bond( "Fe:1u  Fe:1u  +z", VAR("t") );
		//	add_bond( "Fe:1d  Fe:1d  +x", VAR("t") );
		//	add_bond( "Fe:1d  Fe:1d  +y", VAR("t") );
		//	add_bond( "Fe:1d  Fe:1d  +z", VAR("t") );
		//	                   
		//	add_bond( "Fe:1u  Fe:1u  -x", VAR("t") );
		//	add_bond( "Fe:1u  Fe:1u  -y", VAR("t") );
		//	add_bond( "Fe:1u  Fe:1u  -z", VAR("t") );
		//	add_bond( "Fe:1d  Fe:1d  -x", VAR("t") );
		//	add_bond( "Fe:1d  Fe:1d  -y", VAR("t") );
		//	add_bond( "Fe:1d  Fe:1d  -z", VAR("t") );
		//}
		
		
	}

	void render()			{
		// Set the electron carrier for different atoms
		setElectronCarrier	(" Fe:1 ");
		initElectronDensity	(" Fe:1 ");
		
		Mu			= VAR("Mu",0).real();			// Setup chemical potential from input file.
		Temperature = VAR("Temperature",0).real();// Setup Temperature from input file.

		// Get the reciprocal lattice vector
		auto	B = Lat.get_reciprocal();
		auto	b1=B.row(0)*0.5;
		auto	b2=B.row(1)*0.5;
		auto	b3=B.row(2)*0.5;
		
		// Set k-points for general purpose calculation
		clear_k_point();
		unsigned N1=VAR("Nb1").real(), N2=VAR("Nb2").real(), N3=VAR("Nb3").real();
		for (double i1=0; i1<N1; i1++)
		for (double i2=0; i2<N2; i2++)
		for (double i3=0; i3<N3; i3++) {
			add_k_point( (i1/N1)*b1 + (i2/N2)*b2 + (i3/N3)*b3 );
		}

		if (VAR("isCalculateMu").real() == 1) { calcChemicalPotential(); }
		if (VAR("isCalculateVar").real() == 1){ calcVariation(); }
		if (VAR("isCalculateTotE", 0).real() == 1){ calcTotalEnergy(); }

		if (VAR("isCalculateLDOS").real() == 1) {
			selectLDOSsite(0);
			selectLDOSsite(1);
			selectLDOSsite(2);
			calcLDOS(0.01, 0.1);
		}

		// Calculate the band structure through high symmetry points
		if (VAR("isCalculateBand").real() == 1) {
			add_ksymm_point( "", +0.0*b1 +0.0*b2 +0.0*b3);
			add_ksymm_point( "", +0.0*b1 +0.0*b2 +1.0*b3);
			add_ksymm_point( "", +0.0*b1 +1.0*b2 +1.0*b3);
			add_ksymm_point( "", +1.0*b1 +1.0*b2 +1.0*b3);
			add_ksymm_point( "", +0.0*b1 +0.0*b2 +0.0*b3);
			
			calculateBandStructure();
		}
	}
};

int main(int argc, char** argv) {
	// Handle program arguments
	vector<string> args;
	for (int i=0; i<argc; i++) args.push_back(string(argv[i]));
	
	// Construct lattice structure(class) from input.
	string	filename="";
	if		(args.size()==1) { filename = "cubic.lat"; }
	else if (args.size()==2) { filename = args[1]; }
	
	if (filename.size()>0) {
		string cmd = "~/GitHub/tbm-cube/bin/lift.py "+filename;
		cout<<system(cmd.c_str())<<endl;
		
		Lattice Lat(filename, NORMAL);	//Load the lattice file.
		MyModel tbm(Lat);				//render the model
		
	}

	return 0;
}

