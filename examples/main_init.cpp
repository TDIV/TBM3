
#include <stdlib.h>
#include "lift.h"

using namespace lift;

// Define the model
class MyModel: public TBModel{
public:
	MyModel(Lattice lat):	TBModel(lat){ render(); }

	void init_order()		{
		// Initialize the AFM order
		
		Mu = VAR("Mu").real();
		cout<<Mu<<endl;
		initElectronDensity	(" Bi:0  LS:0    Fe:2  Mn:0.7 O:2");
		saveElectronFilling(initDenOrder);
		
		srand(time(NULL));
		while (site_iterate()) {
			auto si = getSite();
			
			initOrder(si, si.Name()+" Sx")	= 0;
			initOrder(si, si.Name()+" Sy")	= 0;
			initOrder(si, si.Name()+" Sz")	= 0;
			
			double Fe_canting_y_to_x = 0;
			double cant_spinA = VAR("cant_spinA").real();
			double cant_spinB = VAR("cant_spinB").real();
			if (si.Name() == "Fe" and Lat(si, "-x-x").Name() == "Mn") { Fe_canting_y_to_x = cant_spinA; }
			if (si.Name() == "Fe" and Lat(si, "+x+x").Name() == "Mn") { Fe_canting_y_to_x = cant_spinB; }
			
			int Mn_spin = (int)VAR("Mn_spin").real();
			int Fe_spin = (int)VAR("Fe_spin").real();
			
			if (VAR("has_mag").real() != 0)
				if (si.Name() == "Fe" or si.Name() == "Mn") {
					int rx=-1, ry=-1, rz=-1;
					if (int(si.pos[0])%2 ) rx=1;
					if (int(si.pos[1])%2 ) ry=1;
					if (int(si.pos[2])%2 ) rz=1;
					
					int mag_sign = rx*ry*rz;
					
					double theta=0, phi=0;
					double mn_theta_up=0, mn_phi_up=0;
					double mn_theta_dn=0, mn_phi_dn=0;
					
					double fe_theta_up=0, fe_phi_up=0;
					double fe_theta_dn=0, fe_phi_dn=0;
					
					if (Fe_spin == 1) {
						fe_theta_up =pi/2;
						fe_phi_up	=0;
						fe_theta_dn =pi/2;
						fe_phi_dn	=pi;
					}
					
					if (Fe_spin == 2) {
						fe_theta_up = pi/2;
						fe_phi_up	= pi/2+Fe_canting_y_to_x;
						fe_theta_dn = pi/2;
						fe_phi_dn	=-pi/2-Fe_canting_y_to_x;
					}
					
					if (Fe_spin == 3) {
						fe_theta_up =0;
						fe_phi_up	=0;
						fe_theta_dn =pi;
						fe_phi_dn	=0;
					}
					
					// ----------- Mn ------------
					if (Mn_spin == 1) {
						mn_theta_up =pi/2;
						mn_phi_up	=0;
						mn_theta_dn =pi/2;
						mn_phi_dn	=0;
					}
					
					if (Mn_spin == 2) {
						mn_theta_up =pi/2;
						mn_phi_up	=pi/2;
						mn_theta_dn =pi/2;
						mn_phi_dn	=pi/2;
					}
					
					if (Mn_spin == 3) {
						mn_theta_up =0;
						mn_phi_up	=0;
						mn_theta_dn =0;
						mn_phi_dn	=0;
					}
					
					if (si.Name() == "Fe" != 0){
						if (mag_sign > 0) {
							theta	= fe_theta_up;
							phi		= fe_phi_up;
						} else if (mag_sign < 0) {
							theta	= fe_theta_dn;
							phi		= fe_phi_dn;
						}
					}
					if (si.Name() == "Mn" != 0){
						if (mag_sign > 0) {
							theta	= mn_theta_up;
							phi		= mn_phi_up;
						} else if (mag_sign < 0) {
							theta	= mn_theta_dn;
							phi		= mn_phi_dn;
						}
					}

					double	Sx = sin(theta)*cos(phi);
					double	Sy = sin(theta)*sin(phi);
					double	Sz = cos(theta);
					
					initOrder(si, si.Name()+" Sx")	=Sx;
					initOrder(si, si.Name()+" Sy")	=Sy;
					initOrder(si, si.Name()+" Sz")	=Sz;
				}
			
		}
		
		cout<<Lat.filename<<endl;
		initOrder.save(Lat, "init");
		initOrder.save(Lat, "");
	}
	
	void Hamiltonian()		{ }

	void render()			{
		init_order();
	}
};

int main(int argc, char** argv) {
	
	// Handle program arguments
	vector<string> args;
	for (int i=0; i<argc; i++) args.push_back(string(argv[i]));
	
	// Construct lattice structure(class) from input.
	string	filename="";
	if		(args.size()==1) { filename = "../input"; }
	else if (args.size()==2) { filename = args[1]; }
	
	if (filename.size()>0) {
		
		string cmd = "lift.py "+filename;
		cout<<system(cmd.c_str())<<endl;
		
		Lattice Lat(filename, NORMAL);
		
		//render the model
		MyModel tbm(Lat);
		
	}
	
	
	return 0;
}





















