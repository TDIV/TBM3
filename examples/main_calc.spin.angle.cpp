
#include "lift.h"

using namespace lift;

double gaussian_white(){
	
	const double epsilon = std::numeric_limits<double>::min();
	const double tau = 2.0*3.14159265358979323846;
	
	static double z0;
	
	double u1, u2;
	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	}
	while ( u1 <= epsilon );

	z0 = sqrt(-2.0 * log(u1)) * cos(tau * u2);
	return z0;
}

// Define the model
class MyModel: public TBModel{
public:
	MyModel(Lattice lat):	TBModel(lat){ render(); }
	
	void init_order()	{ }

	void Hamiltonian() { }

	void setVariational(){ }
	
	void render(){
		OrderParameter orderA, orderB, orderC, orderAngle;
		
		orderA.load(Lat, "init");
		orderB.load(Lat, "");
		
		orderC = orderA-orderB;
		
		while (site_iterate()) {
			auto si = getSite();
			if (si.Name() == "Fe" or si.Name() == "Mn") {
				double dSx=0, dSy=0, dSz=0;
				dSx = orderC(si, si.Name()+" Sx").real();
				dSy = orderC(si, si.Name()+" Sy").real();
				dSz = orderC(si, si.Name()+" Sz").real();
				
				double abs_spin = sqrt(dSx*dSx + dSz*dSz + dSz*dSz);
				double angle = 2*asin(abs_spin/2);
				angle = angle*360/(2*pi);
				orderAngle(si, si.Name()+" ang") = angle;
			}
		}
		orderAngle.save(Lat, "angle");
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
		
		Lattice Lat(filename, NORMAL);
		
		//render the model
		MyModel tbm(Lat);
		
	}
	
	
	return 0;
}





















