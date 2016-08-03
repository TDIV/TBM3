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
//  header.h
//  TBM^3
//
//  Created by Yuan Yen Tai on 9/10/15.
//


#ifndef _header_h
#define _header_h


#include <vector>
#include <deque>
#include <map>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

//#pragma once

#if defined __GNUC__
#   pragma GCC system_header
#elif defined __SUNPRO_CC
#   pragma disable_warn
#elif defined _MSC_VER
#   pragma warning(push)
#   pragma warning(disable : 4800)
#endif

// BEGIN INCLUDING 3RD PARTY HEADERS
#include <boost/geometry/index/rtree.hpp>
#include <boost/tuple/tuple.hpp>

#if defined __SUNPRO_CC
#   pragma enable_warn
#elif defined _MSC_VER
#   pragma warning(pop)
#endif


#include "../Gramat/include/gramat.h"

namespace tbm{
	
using namespace std;
using namespace gmt;
	
	//typedef		dvar	r_var;
	//typedef		dmat	r_mat;
	//typedef		zvar	x_var;
	//typedef		zmat	x_mat;

	typedef		svar	r_var;
	typedef		smat	r_mat;
	typedef		cvar	x_var;
	typedef		cmat	x_mat;
	
void ErrorMessage(string errorStr) {
	cout<<errorStr<<endl<<endl;
	exit(0);
}

#include "string_func.hpp"
#include "math_func.hpp"
#include "read_lat.hpp"
#include "atom.hpp"
#include "lattice.hpp"
#include "order_parameter.hpp"
#include "TBWannierHandler.hpp"
#include "TBDataSource.hpp"
#include "TBModelBase.hpp"
#include "TBClassicalSpinBase.hpp"

}

#endif
