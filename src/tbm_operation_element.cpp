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
//  tbm_operation_element.cpp
//  TBM^3
//
//  Created by Yuan Yen Tai on 9/29/15.
//

/* -------------------------------------------------------------------
 The tbm_operation_element is a buffer and will be called inside the super class "tbm.cpp"
 to help the construction of the Hamiltonain.
 -------------------------------------------------------------------*/



// -----------------------------------
// Make a linear points between two different k-points.
vector<r_mat> make_line(r_mat point_a, r_mat point_b, unsigned Nsteps=10){
	vector<r_mat> ll;
	for (double t=0; t<1; t+=1.0/Nsteps) {
		r_mat position= (point_b-point_a)*t+point_a;
		ll.push_back(position);
	}
	
	return ll;
}
// -----------------------------------

class	Element{
public:
	int		I;
	int		J;
	x_var	val;
	x_var	val_conj;
	string	info;
	string	opt;
	
	Element& operator()(int i, int j, x_var _val, string _info="success"){
		I=i;
		J=j;
		val=_val;
		val_conj = conj(_val);
		info=_info;
		return *this;
	}
};

class	SiteOperationUnit{
public:
	Atom		SiteI;
	string		opt;
	x_var		val;
	SiteOperationUnit(	Atom _SiteI,	string _opt,	x_var _val):
					SiteI(_SiteI),	opt(_opt),		val(_val)
	{ }
};

class	PairOperationUnit{
public:
	AtomPair	pit;
	string		opt;
	x_var		val;
	PairOperationUnit(	AtomPair ap,	string _opt,	x_var _val):
						pit(ap),	opt(_opt),		val(_val)
	{ }
};



