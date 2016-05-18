
#include <stdlib.h>
#include "lift.h"

using namespace lift;

// Define the model
class MyModel: public TBModel{
public:
	MyModel(Lattice lat):	TBModel(lat){ render(); }
	
	vector<AtomPair>	coulomb_iteration;
	
	void init_order()		{
		// Initialize the AFM order
		initOrder.load(Lat, "");
		loadElectronFilling();
		
		// Construct the full pack of lattice bond around an Atom.
		//LatticeBond		lat_bnd(Lat, 2, 2.01);
		
		// Construct the long-range iteration
		while (pair_iterate()) {
			auto si = getPair().AtomI;
			auto sj = getPair().AtomJ;
			auto bond = getPair().bvec;
			
			r_mat bb(1,3);
			bb[0] = bond[0];
			bb[1] = bond[1];
			bb[2] = bond[2];
				AtomPair ap;
				ap.AtomI = si;
				ap.AtomJ = sj;
				ap.bvec	 = bb;
				if (	si.Name() != "Bi" and si.Name() != "LS" and si.Name() != "VA"
					and sj.Name() != "VA") {
					coulomb_iteration.push_back(ap);
					//cout<<bb<<" "<<si.Name()<<"-"<<si.AtomIndex()<<" "<<sj.Name()<<"-"<<sj.AtomIndex()<<endl;
				}
		}
	}
	
	void Hamiltonian()		{
		// Orbital d_{3z^2-r^2}	: 1
		// Orbital d_{x^2-y^2}	: 2
		
		add_Chemical_Potential();
		
		OrderParameter order = initOrder;
		
		map<string, int> OO_bond_sign;
		OO_bond_sign["+x+y"] = -1;
		OO_bond_sign["+x-y"] = +1;
		OO_bond_sign["-x+y"] = +1;
		OO_bond_sign["-x-y"] = -1;
		
		OO_bond_sign["+x+z"] = -1;
		OO_bond_sign["+x-z"] = +1;
		OO_bond_sign["-x+z"] = +1;
		OO_bond_sign["-x-z"] = -1;
		
		OO_bond_sign["+y+z"] = -1;
		OO_bond_sign["+y-z"] = +1;
		OO_bond_sign["-y+z"] = +1;
		OO_bond_sign["-y-z"] = -1;
		
		map<string, double> Coulomb_Z;
		Coulomb_Z["Fe"]	= 5;
		Coulomb_Z["Mn"]	= 4;
		Coulomb_Z["O"]	= 0;
		Coulomb_Z["Bi"]	= 3;
		Coulomb_Z["LS"]	= 2.7;
		
		//---site iteration---
		while (site_iterate()) {
			auto si = getSite();
			
			double Jh=0;
			if (si.Name() == "Fe" ) Jh=-VAR("Jh_fe").real();
			if (si.Name() == "Mn" ) Jh=-VAR("Jh_mn").real();
			
			add_hund_spin("Fe 1", Jh, order.getVars(si, "Sx Sy Sz"));
			add_hund_spin("Fe 2", Jh, order.getVars(si, "Sx Sy Sz"));
			add_hund_spin("Mn 1", Jh, order.getVars(si, "Sx Sy Sz"));
			add_hund_spin("Mn 2", Jh, order.getVars(si, "Sx Sy Sz"));
			
			add_site("O  1u 1u delta");
			add_site("O  1d 1d delta");
			add_site("Mn 1u 1u onMn");
			add_site("Mn 1d 1d onMn");
			add_site("Mn 2u 2u onMn");
			add_site("Mn 2d 2d onMn");
		}
		
		OrderParameter	OnSiteCoulombOrder;
		//---pair iteration---
		while (pair_iterate()) {
			auto pit = getPair();
			auto si = pit.AtomI;
			auto sj = pit.AtomJ;
			
			//---bond iteration---
			for ( map<string, int>::iterator it = OO_bond_sign.begin(); it!=OO_bond_sign.end(); it++) {
				add_bond("O:1u  O:1u "+it->first, it->second*VAR("tpp") );
				add_bond("O:1d  O:1d "+it->first, it->second*VAR("tpp") );
			}
			
			if (si.Name() == "Fe" or si.Name() == "Mn") {
				x_var adjust_t = 0;
				if (si.Name() == "Fe") adjust_t = VAR("t_fe");
				if (si.Name() == "Mn") adjust_t = VAR("t_mn");
				
				add_bond_hc( si.Name()+":1u  O:1u  +x", +VAR("t1x")*adjust_t );
				add_bond_hc( si.Name()+":1u  O:1u  +y", +VAR("t1y")*adjust_t );
				add_bond_hc( si.Name()+":1u  O:1u  +z", +VAR("t1z")*adjust_t );
				add_bond_hc( si.Name()+":1d  O:1d  +x", +VAR("t1x")*adjust_t );
				add_bond_hc( si.Name()+":1d  O:1d  +y", +VAR("t1y")*adjust_t );
				add_bond_hc( si.Name()+":1d  O:1d  +z", +VAR("t1z")*adjust_t );
				
				add_bond_hc( si.Name()+":2u  O:1u  +x", +VAR("t2x")*adjust_t );
				add_bond_hc( si.Name()+":2u  O:1u  +y", +VAR("t2y")*adjust_t );
				add_bond_hc( si.Name()+":2u  O:1u  +z", +VAR("t2z")*adjust_t );
				add_bond_hc( si.Name()+":2d  O:1d  +x", +VAR("t2x")*adjust_t );
				add_bond_hc( si.Name()+":2d  O:1d  +y", +VAR("t2y")*adjust_t );
				add_bond_hc( si.Name()+":2d  O:1d  +z", +VAR("t2z")*adjust_t );
				
				add_bond_hc( si.Name()+":1u  O:1u  -x", -VAR("t1x")*adjust_t );
				add_bond_hc( si.Name()+":1u  O:1u  -y", -VAR("t1y")*adjust_t );
				add_bond_hc( si.Name()+":1u  O:1u  -z", -VAR("t1z")*adjust_t );
				add_bond_hc( si.Name()+":1d  O:1d  -x", -VAR("t1x")*adjust_t );
				add_bond_hc( si.Name()+":1d  O:1d  -y", -VAR("t1y")*adjust_t );
				add_bond_hc( si.Name()+":1d  O:1d  -z", -VAR("t1z")*adjust_t );
				
				add_bond_hc( si.Name()+":2u  O:1u  -x", -VAR("t2x")*adjust_t );
				add_bond_hc( si.Name()+":2u  O:1u  -y", -VAR("t2y")*adjust_t );
				add_bond_hc( si.Name()+":2u  O:1u  -z", -VAR("t2z")*adjust_t );
				add_bond_hc( si.Name()+":2d  O:1d  -x", -VAR("t2x")*adjust_t );
				add_bond_hc( si.Name()+":2d  O:1d  -y", -VAR("t2y")*adjust_t );
				add_bond_hc( si.Name()+":2d  O:1d  -z", -VAR("t2z")*adjust_t );
			}
			
			// Classical spin part
			r_var	Js=0;
			if (si.Name() == "Fe" and sj.Name()== "Fe") Js = VAR("Js_fe").real();
			if (si.Name() == "Mn" and sj.Name()== "Mn") Js = VAR("Js_mn").real();
			if (si.Name() == "Fe" and sj.Name()== "Mn") Js = VAR("Js_fe_mn").real();
			if (si.Name() == "Mn" and sj.Name()== "Fe") Js = VAR("Js_fe_mn").real();
			
			add_classical_spin( "Fe Fe "+pit.bond, Js, order.getVars(si, "Sx Sy Sz"), order.getVars(sj, "Sx Sy Sz"));
			add_classical_spin( "Mn Mn "+pit.bond, Js, order.getVars(si, "Sx Sy Sz"), order.getVars(sj, "Sx Sy Sz"));
			add_classical_spin( "Fe Mn "+pit.bond, Js, order.getVars(si, "Sx Sy Sz"), order.getVars(sj, "Sx Sy Sz"));
			add_classical_spin( "Mn Fe "+pit.bond, Js, order.getVars(si, "Sx Sy Sz"), order.getVars(sj, "Sx Sy Sz"));
		}
		
		// Construct local Coulomb potential profile.
		for (unsigned i=0; i<coulomb_iteration.size(); i++) {
			auto si = coulomb_iteration[i].AtomI;
			auto sj = coulomb_iteration[i].AtomJ;
			auto distance = coulomb_iteration[i].bvec;
			
			double	nDen  = initDenOrder(sj, sj.Name()+" den").real();
		  	double	alpha = VAR("alpha").real();
		  	double	ion_Z = Coulomb_Z[sj.Name()];
		  	
		  	double coulomb = Hatree_Coulomb_Potential(alpha, nDen - ion_Z, distance);
		  	OnSiteCoulombOrder(si, si.Name()+" V")+=coulomb;
		}
		
		OnSiteCoulombOrder.save(Lat, "V");
		
		// Using the On-site Coulomb order to adjust local chemical potential
		while (site_iterate()) {
			auto si = getSite();
			add_site(si.Name()+" 1u 1u", OnSiteCoulombOrder(si, si.Name()+" V") );
			add_site(si.Name()+" 2u 2u", OnSiteCoulombOrder(si, si.Name()+" V") );
			add_site(si.Name()+" 1d 1d", OnSiteCoulombOrder(si, si.Name()+" V") );
			add_site(si.Name()+" 2d 2d", OnSiteCoulombOrder(si, si.Name()+" V") );
		}
	}
	
	void render()			{
		// Set the electron carrier for different atoms
		setElectronCarrier	(" Bi:3  LS:2.7  Fe:5  Mn:4   O:0");
		initElectronDensity	(" Bi:0  LS:0    Fe:2  Mn:0.7 O:2");
		
		// Get the reciprocal lattice vector
		auto	B = Lat.get_reciprocal();
		auto	b1=B.row(0)*0.5;
		auto	b2=B.row(1)*0.5;
		auto	b3=B.row(2)*0.5;
		Mu			= VAR("Mu").real();			// Setup chemical potential from input file.
		Temperature = VAR("Temperature").real();// Setup Temperature from input file.
		
		init_order();
		
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
		
		auto var_list = KHamEvd();
		auto min_eig = var_list[0];
		auto max_eig = var_list[1];
		auto totalE  = var_list[2];
		cout<<min_eig<<" "<<max_eig<<" "<<totalE<<endl;
		
		if (VAR("isCalculateLDOS").real() == 1) {
			selectLDOSsite(1);
			selectLDOSsite(2);
			calcLDOS(0.01, 0.04);
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




