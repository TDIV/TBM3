/* ------------------------------------------------------- */
/* -Output format of the matrix--------------------------- */
/* ------------------------------------------------------- */

string tostr(unsigned val){
	ostringstream oss;
	oss.precision(PRINT_PRECISION);
	oss<<val;
	return oss.str();
}
string tostr(float val){
	ostringstream oss;
	oss.precision(PRINT_PRECISION);
	oss<<val;
	return oss.str();
}
string tostr(double val){
	ostringstream oss;
	oss.precision(PRINT_PRECISION);
	oss<<val;
	return oss.str();
}
string tostr(cvar val){
	ostringstream oss;
	oss<<val.tostr();
	return oss.str();
}
string tostr(zvar val){
	ostringstream oss;
	oss<<val.tostr();
	return oss.str();
}
/* ------------------------------------------------------- */
/* ------------------------------------------------------- */



string fformat(unsigned val,size_t len=10){
	stringstream ss(stringstream::in | stringstream::out);
	ss << val;
	string str;
	str = ss.str();

	while(str.size() < len)
		str.push_back(' ');

	return str;
}
/* ------------------------------------------------------- */
string fformat(float val,size_t len=10){
	stringstream ss(stringstream::in | stringstream::out);
	ss << val;
	string str;
	if (val >= 0)	str = " "+ss.str();
	else			str = ss.str();
	
	if (abs(val)<0.000000000001) str = " 0";

	while(str.size() < len)
		str.push_back(' ');

	return str;
}
/* ------------------------------------------------------- */
string fformat(double val,size_t len=10){
	stringstream ss(stringstream::in | stringstream::out);
	ss << val;
	string str;
	if (val >= 0)	str = " "+ss.str();
	else			str = ss.str();

	if (abs(val)<0.000000000001) str = " 0";
	
	while(str.size() < len)
		str.push_back(' ');

	return str;
}
/* ------------------------------------------------------- */
string fformat(cvar val,size_t len=10){
	
	string str = val.tostr();
	while(str.size() < len)
		str.push_back(' ');

	return str;
}
/* ------------------------------------------------------- */
string fformat(zvar val,size_t len=10){
	
	string str = val.tostr();
	while(str.size() < len)
		str.push_back(' ');

	return str;
}
/* ------------------------------------------------------- */
string fformat(string val,size_t len=10){
	stringstream ss(stringstream::in | stringstream::out);
	ss << val;
	string str;
	str = ss.str();

	while(str.size() < len)
		str.push_back(' ');

	return str;
}

/* ------------------------------------------------------- */
/* ------------------------------------------------------- */

template<class T>
string fmt(T val, size_t len=10){ return fformat(val,len); }

