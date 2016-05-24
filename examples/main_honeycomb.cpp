
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
		
		Lat.bondStringMap["+a"] = "+1..";
		Lat.bondStringMap["+b"] = ".+1.";
		Lat.bondStringMap["+c"] = "-1-1.";
		Lat.bondStringMap["-a"] = "-1..";
		Lat.bondStringMap["-b"] = ".-1.";
		Lat.bondStringMap["-c"] = "+1+1.";
	
	}
	
	void Hamiltonian()		{
		add_Chemical_Potential();
		
		OrderParameter order = initOrder;

		//constructHoppingHamiltonian(Ham, k_space);
		
		//---pair iteration---
		while (pair_iterate()) {
			auto pit = getPair();
			
			add_bond_hc( "C:1u  C:1u  +a", VAR("hop") );
			add_bond_hc( "C:1u  C:1u  +b", VAR("hop") );
			add_bond_hc( "C:1u  C:1u  +c", VAR("hop") );
			
			add_bond_hc( "C:1d  C:1d  +a", VAR("hop") );
			add_bond_hc( "C:1d  C:1d  +b", VAR("hop") );
			add_bond_hc( "C:1d  C:1d  +c", VAR("hop") );
		}
		
	}

	void render()			{
		// Set the electron carrier for different atoms
		setElectronCarrier	(" C:1 ");
		initElectronDensity	(" C:1 ");
		
		Mu			= VAR("Mu",0).real();			// Setup chemical potential from input file.
		Temperature = VAR("Temperature",0).real();// Setup Temperature from input file.

		// Get the reciprocal lattice vector
		auto	B = Lat.get_reciprocal();
		auto	b1=B.row(0);
		auto	b2=B.row(1);
		auto	b3=B.row(2);
		auto	k = b1;
		k[1] = k[0]/sqrt(3);
		cout<<k<<endl;
		
		// Set k-points for general purpose calculation
		clear_k_point();
		unsigned N1=VAR("Nb1").real(), N2=VAR("Nb2").real(), N3=VAR("Nb3").real();
		for (double i1=0; i1<N1; i1++)
		for (double i2=0; i2<N2; i2++)
		for (double i3=0; i3<N3; i3++) {
			add_k_point( (i1/N1)*b1 + (i2/N2)*b2 + (i3/N3)*b3 );
		}

		if (VAR("isCalculateMu", 0).real() == 1) { calcChemicalPotential(); }
		if (VAR("isCalculateVar", 0).real() == 1){ calcVariation(); }
		if (VAR("isCalculateTotE", 0).real() == 1){ calcTotalEnergy(); }

		if (VAR("isCalculateLDOS", 0).real() == 1) {
			selectLDOSsite(0);
			selectLDOSsite(1);
			selectLDOSsite(2);
			calcLDOS(0.01, 0.1);
		}

		// Calculate the band structure through high symmetry points
		if (VAR("isCalculateBand").real() == 1) {
			add_ksymm_point( "",  0.0*b1 + 0.0*b2);
			add_ksymm_point( "",  0.5*b1 + 0.0*b2);
			add_ksymm_point( "",  1.0*k + 0.0*b2);
			add_ksymm_point( "",  0.0*b1 + 0.0*b2);
			
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
	if		(args.size()==1) { filename = "honeycomb.lat"; }
	else if (args.size()==2) { filename = args[1]; }
	
	if (filename.size()>0) {
		string cmd = "~/GitHub/tbm-cube/bin/lift.py "+filename;
		cout<<system(cmd.c_str())<<endl;
		
		Lattice Lat(filename, NORMAL);	//Load the lattice file.
		MyModel tbm(Lat);				//render the model
		
	}

	return 0;
}

