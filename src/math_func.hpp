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
//  math_func.h
//  TBM^3
//
//  Created by Yuan-Yen Tai on 9/10/15.
//

#ifndef _math_func_h
#define _math_func_h

// -----------------------------------
// Gaussian white noise for LLG spin dynamics
double gaussian_white(){
	// Reference : Box-Muller algorithm on Wiki.
	
	//const double epsilon = std::numeric_limits<double>::min();
	const double tau = 2.0*3.14159265358979323846;
	
	static double z0;
	
	double u1, u2;
	//do
	//{
	u1 = rand() * (1.0 / RAND_MAX);
	u2 = rand() * (1.0 / RAND_MAX);
	//}
	//while ( u1 <= epsilon );
	
	z0 = sqrt(-2.0 * log(u1)) * cos(tau * u2);
	return z0;
}

double Hatree_Coulomb_Potential(double alpha, double n_Z, r_mat R){
	// Reference : PRB 84, 024422
	
	if (R.size() == 3) {
		double distance = sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
		if (distance > 0.000001) { return alpha*n_Z/distance; }
		else	return 0;
	}
	return 0;
}

double Coulomb_Screening(double r0, r_mat R){
	// Reference : PRB 84, 024422
	
	if (R.size() == 3) {
		double distance = sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
		if (distance > 0.000001) { return exp(-distance/r0); }
		else	return 0;
	}
	return 0;
}

r_mat vec(r_var x, r_var y, r_var z){
	r_mat bond(1,3);
	bond[0] = x;
	bond[1] = y;
	bond[2] = z;
	return bond;
}

x_mat xvec(x_var x)						{
	x_mat bond(1,1);
	bond[0] = x;
	return bond;
}
x_mat xvec(x_var x, x_var y)			{
	x_mat bond(1,2);
	bond[0] = x;
	bond[1] = y;
	return bond;
}
x_mat xvec(x_var x, x_var y, x_var z)	{
	x_mat bond(1,3);
	bond[0] = x;
	bond[1] = y;
	bond[2] = z;
	return bond;
}

void normalizeXVec(x_mat & xvec){
	double length = sqrt(cdot(xvec,xvec).real());
	if( length != 0 )
		xvec = xvec * (1/length);
}

string vecToStr(r_mat rvec){
	string vstr = "";
	for(unsigned i=0 ; i<rvec.size() ; i++){
		if(rvec[i] >= 0){ vstr += "+"; }
		else			{ vstr += "-"; }
		vstr += DoubleToStr(abs(rvec[i]));
	}
	return vstr;
}

vector<r_mat> getUnitVectors(vector<r_mat> AVec){
	vector<r_mat> retVec;
	for( auto & vec: AVec){
		double vec_len = sqrt(cdot(vec,vec));
		vec = vec * (1/vec_len);
		retVec.push_back(vec);
	}
	return retVec;
}

#endif

