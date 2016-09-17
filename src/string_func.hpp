/*-----------------------------------------------------------|
| Copyright (C) 2016 Yuan-Yen Tai, Hongchul Choi,            |
|                    Jian-Xin Zhu                            |
|                                                            |
| This file is distributed under the terms of the BSD        |
| Berkeley Software Distribution. See the file `LICENSE' in  |
| the root directory of the present distribution, or         |
| https://en.wikipedia.org/wiki/BSD_licenses.                |
|                                                            |
|-----------------------------------------------------------*/
//
//  string_func.hpp
//  TBM^3
//
//  Created by Yuan Yen Tai on 7/15/15.
//

#ifndef _string_func_h
#define _string_func_h

void			removeSpace			(std::string & str){
	if (str.size() > 0) {
		string tmpStr = "";
		for( unsigned i=0 ; i<str.size() ; i++){
			if (str[i] != ' ' and str[i] != '\t') {
				tmpStr += str[i];
			}
		}
		str = tmpStr;
	}
}
void			removeSpaceTopToe	(std::string & str){
	while( (str[0] == ' ' or str[0] == '\t') and str.size() > 0){
		str.erase(str.begin());
	}
	while( (str[str.size()-1] == ' ' or str[str.size()-1] == '\t') and str.size() > 0){
		str.pop_back();
	}
}

void			replaceAll		(std::string& str, const std::string& from, const std::string& to) {
	// The string replace function.
	if(from.empty())
		return;
	size_t start_pos = 0;
	while((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
	}
}
void			deleteComment	(std::string& str, char comment = '%'){
	string tmpStr = "";
	for(auto & c: str){
		if( c!= comment) tmpStr.push_back(c);
		else{
			str = tmpStr;
			return;
		}
	}
}
deque<string>	split			(string line, string delimiter){
	// The string split function.
	replaceAll(line, "\t", " ");
    std::deque<std::string> content;
    
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
deque<string>	splitByIntNumber(string line)		{
	// The string split function.
	replaceAll(line, "\t", " ");
    std::deque<std::string> content;
	line += " ";
	
	string word = "";
	for( unsigned i=0 ; i<line.size()-1 ; i++){
		char c = line[i];
		char cnext = line[i+1];
		
		if		( (c>= '0' and c<='9') or c=='.'){
			word.push_back(c);
			if( (cnext < '0' or cnext > '9' or cnext == ' ') and cnext != '.'){
				content.push_back(word);
				word = "";
			}
		}
		else if	( (c < '0' or c > '9') and c!='.'){
			word.push_back(c);
			if( (cnext >= '0' and cnext <='9') or cnext == ' '){
				content.push_back(word);
				word = "";
			}
		}
	}
	
    return content;
}


/* String-Variable convertors. */
int				StrToInt		(string Str)	{
    istringstream iss(Str);
    int var;
    iss>>var;
    return var;
}
string			IntToStr		(int val)		{
    ostringstream oss;
    oss<<val;
    return oss.str();
}
double			StrToDouble		(string Str)	{
	istringstream iss(Str);
	double var;
	iss>>var;
	return var;
}
string			DoubleToStr		(double val)	{
    ostringstream oss;
    oss<<val;
    return oss.str();
}
x_var			StrToComplex	(string Str)	{
	replaceAll(Str, " ", "");
	
	x_var ret_var;
	if (Str[0] == '(' and Str[Str.size()-1] == ')') {
		Str[0] = ' ';
		Str[Str.size()-1] = ' ';
		auto var = split(Str, ",");
		if		( var.size() == 1) ret_var = StrToDouble(var[0]);
		else if	( var.size() == 2) ret_var = x_var( StrToDouble(var[0]), StrToDouble(var[1]));
		else{
			ErrorMessage("Error, string:"+Str+" is not a valid from of complex variable.");
		}
		
	} else {
		ret_var = StrToDouble(Str);
	}
    return ret_var;
}
x_mat	&		StrToXVec		(string Str)	{
	/* Convert string type: "(1,2), 2,3, (3,4)" into a 1xN x_mat.*/
	bool xflag = false;
	vector<pair<bool, string> > subList;
	
	string tmpStr = "";
	for( auto & c : Str ){
		
		if( c == '(' ){ xflag = true;	tmpStr.push_back(c);	continue;}
		if( c == ')' ){ xflag = false;	tmpStr.push_back(c);	continue;}
		
		if( c!= ' ' and c!= '\t') {
			// Complex number goes here.
			if(xflag)	{ tmpStr.push_back(c); }
			else		{
				if( c!=',') tmpStr.push_back(c);
				else		tmpStr.push_back('@');
			}
		}
	}
	auto xList = split(tmpStr, "@");
	
	x_mat *ret = new x_mat(1,xList.size());
	for( unsigned ii=0 ; ii<xList.size() ; ii++){
		auto & xStr = xList[ii];
		ret->index(ii) = StrToComplex(xStr);
	}
	return *ret;
}

bool			IsIntStr(string Str)			{
	bool retVal = true;
	for(auto &c: Str) retVal = retVal and ( (c>='0' and c<='9') or c=='-' or c=='+');
	return retVal;
}
bool			IsFloatStr(string Str)			{
	bool retVal = true;
	for(auto &c: Str) retVal = retVal and ( (c>='0' and c<='9') or c=='-' or c=='+' or c=='.');
	return retVal;
}

int				WordCount(string token, string line){
	/* Counting how many words in one string. */
	
	int count=0;
	
    size_t pos = 0;
    while ((pos = line.find(token)) != std::string::npos) {
        line.erase(0, pos + token.length());
		count++;
    }
	
	
	return count;
}

/* Translate "+" or "-" into +1 or -1. */
double			PMStr(string str){ // return +1 or -1 for a given +
	if(str == "+") return +1;
	if(str == "-") return -1;
	return 0;
}

pair<bool,vector<double> > makeBoxFromStr(string pointA, string pointB){
	replaceAll(pointA, "[", "");
	replaceAll(pointA, "]", "");
	replaceAll(pointB, "[", "");
	replaceAll(pointB, "]", "");
	cout<<endl;
	vector<double> pA, pB;
	for( auto ii: split(pointA, ",")){
		pA.push_back(StrToDouble(ii));
	}
	
	for( auto ii: split(pointB, ",")){
		pB.push_back(StrToDouble(ii));
	}
	
	vector<double> returnSixBoxValue;
	if( pA.size() == 3 and pB.size() == 3){
		double min_x = pA[0] < pB[0] ? pA[0] : pB[0];
		double max_x = pA[0] > pB[0] ? pA[0] : pB[0];
		
		double min_y = pA[1] < pB[1] ? pA[1] : pB[1];
		double max_y = pA[1] > pB[1] ? pA[1] : pB[1];
		
		double min_z = pA[2] < pB[2] ? pA[2] : pB[2];
		double max_z = pA[2] > pB[2] ? pA[2] : pB[2];
		
		returnSixBoxValue.push_back(min_x);
		returnSixBoxValue.push_back(min_y);
		returnSixBoxValue.push_back(min_z);
		returnSixBoxValue.push_back(max_x);
		returnSixBoxValue.push_back(max_y);
		returnSixBoxValue.push_back(max_z);
		
		return make_pair(true,returnSixBoxValue);
	}
	
	return make_pair(false,returnSixBoxValue);
}

#endif











