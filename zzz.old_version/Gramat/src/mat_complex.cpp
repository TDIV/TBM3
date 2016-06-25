/* The purpose of these class is to access 
	cuFloatComplex or cuDoubleComplex */


//--------------- --------------- --------------- --------------- --------------- --------------- ---------------

class zvar;
class cvar;

typedef float			svar ; // using svar for single variable
typedef double			dvar ; // using dvar for double variable

typedef complex<float>	cmplx; // C++ type single complex variable
typedef complex<double> zmplx; // C++ type double complex variable


const dvar pi=3.14159265359;

//--------------- --------------- --------------- --------------- --------------- --------------- ---------------

class zvar: public cuDoubleComplex{
public:
	zvar(dvar re=0, dvar im=0){ x=re; y=im;}
	dvar real()const{return x;}
	dvar imag()const{return y;}
	zvar& operator=(const dvar & val)			{ x=val; y=0; return *this; }
	zvar& operator=(const cuDoubleComplex & val){ x=val.x; y=val.y; return *this; }

	zvar  operator+ ()					const { return zvar(+this->x,+this->y); }
	zvar  operator+ (const zvar	& val)	const { return zvar(x+val.x, y+val.y);	}
	zvar& operator+=(const zvar	& val)	{ x+=val.x; y+=val.y; return *this; }
	zvar& operator+=(const dvar	& val)	{ x+=val;				return *this; }

	zvar  operator- ()					const { return zvar(-this->x,-this->y); }
	zvar  operator- (const zvar	& val)	const { return zvar(x-val.x, y-val.y);	}
	zvar& operator-=(const zvar	& val)	{ x-=val.x; y-=val.y;	return *this; }
	zvar& operator-=(const dvar	& val)	{ x-=val;				return *this; }

	zvar  operator* (const zvar	& val)	const { return zvar(x*val.x-y*val.y, x*val.y+y*val.x); }
	zvar& operator*=(const zvar	& val)	{ dvar ax=x, ay=y; x=(ax*val.x-ay*val.y); y=(ax*val.y+ay*val.x); return *this; }
	zvar& operator*=(const dvar	& val)	{ x*=val; y*=val; return *this; }

	zvar  operator/ (const zvar	& val)	const { return zvar((x*val.x+y*val.y)/(val.x*val.x+val.y*val.y), (y*val.x-x*val.y)/(val.x*val.x+val.y*val.y)); }
	zvar  operator/ (const dvar & val)	const { return zvar(x/val,y/val);}
	zvar& operator/=(const zvar	& val)	{ dvar ax=x, ay=y; x=(ax*val.x+ay*val.y)/(val.x*val.x+val.y*val.y); y=(ay*val.x-ax*val.y)/(val.x*val.x+val.y*val.y); return *this; }
	zvar& operator/=(const dvar	& val)	{ x/=val; y/=val; return *this; }
	string tostr(){
		ostringstream oss;
		oss.precision(PRINT_PRECISION);
		
		double re=real();
		double im=imag();
		if (abs(re) < 0.0000000001) { re=0; }
		if (abs(im) < 0.0000000001) { im=0; }
		oss<<"("<< real() << "," <<imag() <<")";
		return oss.str();
	}
	
};

zvar operator+ (const dvar & val, const zvar & me  ) { return zvar(me.x+val, me.y);}
zvar operator+ (const zvar & me , const dvar & val ) { return zvar(me.x+val, me.y);}
zvar operator- (const dvar & val, const zvar & me  ) { return zvar(val-me.x,-me.y);}
zvar operator- (const zvar & me , const dvar & val ) { return zvar(me.x-val, me.y);}
zvar operator* (const dvar & val, const zvar & me  ) { return zvar(me.x*val,me.y*val);}
zvar operator* (const zvar & me , const dvar & val ) { return zvar(me.x*val,me.y*val);}
zvar operator/ (const dvar & val, const zvar & me  ) { return zvar((val*me.x)/(me.x*me.x+me.y*me.y), (-val*me.y)/(me.x*me.x+me.y*me.y)); } //Need to be rewrite


ostream & operator<< (ostream & out, zvar & m){ out<<"("<<m.x<<","<<m.y<<")"; return out; }

istream & operator>> (istream & in , zvar & m){
	complex<dvar> tmp;
	in>>tmp;
	m.x=tmp.real();
	m.y=tmp.imag();
	return in;
}

//--------------- --------------- --------------- --------------- --------------- --------------- ---------------

class cvar: public cuFloatComplex{
public:
	cvar(zvar val){ x=(float)val.x; y=(float)val.y;}
	cvar(dvar re=0, dvar im=0){ x=re; y=im;}
	dvar real()const{return x;}
	dvar imag()const{return y;}
	cvar& operator=(const dvar & val)			{ x=val; y=0; return *this; }
	cvar& operator=(const cuFloatComplex & val)	{ x=val.x; y=val.y; return *this; }
	cvar& operator=(const zvar & val)			{ x=(float)val.x; y=(float)val.y; return *this;}

	cvar  operator+ ()					const { return cvar(+this->x,+this->y); }
	cvar  operator+ (const cvar	& val)	const { return cvar(x+val.x, y+val.y);	}
	cvar& operator+=(const cvar	& val)	{ x+=val.x; y+=val.y; return *this; }
	cvar& operator+=(const dvar	& val)	{ x+=val;				return *this; }

	cvar  operator- ()					const { return cvar(-this->x,-this->y); }
	cvar  operator- (const cvar	& val)	const { return cvar(x-val.x, y-val.y);	}
	cvar& operator-=(const cvar	& val)	{ x-=val.x; y-=val.y;	return *this; }
	cvar& operator-=(const dvar	& val)	{ x-=val;				return *this; }

	cvar  operator* (const cvar	& val)	const { return cvar(x*val.x-y*val.y, x*val.y+y*val.x); }
	cvar& operator*=(const cvar	& val)	{ dvar ax=x, ay=y; x=(ax*val.x-ay*val.y); y=(ax*val.y+ay*val.x); return *this; }
	cvar& operator*=(const dvar	& val)	{ x*=val; y*=val; return *this; }

	cvar  operator/ (const cvar	& val)	const { return cvar((x*val.x+y*val.y)/(val.x*val.x+val.y*val.y), (y*val.x-x*val.y)/(val.x*val.x+val.y*val.y)); }
	cvar  operator/ (const dvar & val)	const { return cvar(x/val,y/val);}
	cvar& operator/=(const cvar	& val)	{ dvar ax=x, ay=y; x=(ax*val.x+ay*val.y)/(val.x*val.x+val.y*val.y); y=(ay*val.x-ax*val.y)/(val.x*val.x+val.y*val.y); return *this; }
	cvar& operator/=(const dvar	& val)	{ x/=val; y/=val; return *this; }
	string tostr(){
		ostringstream oss;
		oss.precision(PRINT_PRECISION);
		double re=real();
		double im=imag();
		
		if (abs(re) < 0.0000000001) { re=0; }
		if (abs(im) < 0.0000000001) { im=0; }
		oss<<"("<< re << "," <<im <<")";
		return oss.str();
	}
	
	operator zvar()	{return zvar(x,y);}
};

cvar operator+ (const dvar & val, const cvar & me  ) { return cvar(me.x+val, me.y);}
cvar operator+ (const cvar & me , const dvar & val ) { return cvar(me.x+val, me.y);}
cvar operator- (const dvar & val, const cvar & me  ) { return cvar(val-me.x,-me.y);}
cvar operator- (const cvar & me , const dvar & val ) { return cvar(me.x-val, me.y);}
cvar operator* (const dvar & val, const cvar & me  ) { return cvar(me.x*val,me.y*val);}
cvar operator* (const cvar & me , const dvar & val ) { return cvar(me.x*val,me.y*val);}
cvar operator/ (const dvar & val, const cvar & me  ) { return cvar((val*me.x)/(me.x*me.x+me.y*me.y), (-val*me.y)/(me.x*me.x+me.y*me.y)); } //Need to be rewrite


ostream & operator<< (ostream & out, cvar & m){ out<<"("<<m.x<<","<<m.y<<")"; return out; }

istream & operator>> (istream & in , cvar & m){
	complex<dvar> tmp;
	in>>tmp;
	m.x=tmp.real();
	m.y=tmp.imag();
	return in;
}


// Define a pure imaginary variable for convenient operations
const zvar Im(0,1);

	
	

	
	
