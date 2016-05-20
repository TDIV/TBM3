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
//  string_func.h
//  TBM^3
//
//  Created by Yuan Yen Tai on 7/15/15.
//

#ifndef _string_func_h
#define _string_func_h

#include <vector>
#include <string>
using namespace std;


// The string replace function.
void			replaceAll(std::string& str, const std::string& from, const std::string& to) {
	if(from.empty())
		return;
	size_t start_pos = 0;
	while((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
	}
}

// The string split function.
vector<string>	split(string line, string delimiter) {
    std::vector<std::string> content;
    
    size_t pos = 0;
    std::string token;
    while ((pos = line.find(delimiter)) != std::string::npos) {
        token = line.substr(0, pos);
        if (token.size()>0) content.push_back(token);
        line.erase(0, pos + delimiter.length());
    }
    if (line.size()>0) content.push_back(line);
    
    return content;
}

int				StrToInt(string Str) {
    istringstream iss(Str);
    int var;
    iss>>var;
    return var;
}

string			IntToStr(int val) {
    ostringstream oss;
    oss<<val;
    return oss.str();
}

double			StrToDouble(string Str) {
	istringstream iss(Str);
	double var;
	iss>>var;
	return var;
}

x_var			StrToComplex(string Str) {
	replaceAll(Str, " ", "");
	
	x_var ret_var;
	if (Str[0] == '(' and Str[Str.size()-1] == ')') {
		Str[0] = ' ';
		Str[Str.size()-1] = ' ';
		auto var = split(Str, ",");
		ret_var = x_var( StrToDouble(var[0]), StrToDouble(var[1]));
	} else {
		ret_var = StrToDouble(Str);
	}
    return ret_var;
}

int				WordCount(string token, string line){
	
	int count=0;
	
    size_t pos = 0;
    while ((pos = line.find(token)) != std::string::npos) {
        line.erase(0, pos + token.length());
		count++;
    }
	
	
	return count;
}

// This function will extract the atom informations, for:
// 1. Number of orbitals.
// 2. If the atom has spin degree of freedom.
vector<string>	analyze_atom_label(string label){
	vector<string> lab_list;
	string L_num  = "";
	string L_spin = "";
	for (unsigned i=0; i<label.size(); i++) {
		if ( label[i] >= '0' and label[i] <= '9') {
			L_num.push_back(label[i]);
		}
		else{
			L_spin.push_back(label[i]);
		}
	}
	lab_list.push_back(L_num);
	lab_list.push_back(L_spin);
	return lab_list;
}

#endif











