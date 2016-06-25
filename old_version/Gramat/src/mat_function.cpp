#ifndef __MATRIX_FUNCTION_OPERATION__
#define __MATRIX_FUNCTION_OPERATION__

///* ------------------------------------------------------- */
///* -Functional operations of Array like Matrix------------ */
///* ------------------------------------------------------- */
cvar conj(cvar c){ return cvar(c.real(),-c.imag()); }
zvar conj(zvar c){ return zvar(c.real(),-c.imag()); }

//--------------------------------------------------------------
//---Calculate the abs of the complex variable
//--------------------------------------------------------------
double abs(cvar val){ return sqrt(val.real()*val.real()+val.imag()*val.imag()); }
double abs(zvar val){ return sqrt(val.real()*val.real()+val.imag()*val.imag()); }

//--- Vector operations ---
template<class vv>
matrix<vv> curl(matrix<vv> a, matrix<vv> b){
	matrix<vv> z(1,3);
	if (a.size()==3 and b.size()==3 ) {
		z[0]=a[1]*b[2]-a[2]*b[1];
		z[1]=a[2]*b[0]-a[0]*b[2];
		z[2]=a[0]*b[1]-a[1]*b[0];
	}
	return z;
}

template<class vv>
vv cdot(matrix<vv> a, matrix<vv> b){
	vv z=0;
	if (a.size()==3 and b.size()==3 ) {
		z+=a[0]*b[0];
		z+=a[1]*b[1];
		z+=a[2]*b[2];
	}
	return z;
}


//--------------------------------------------------------------
//---Calculate the square root of the variable------------------
//--------------------------------------------------------------
zvar sqrt(zvar _X){
	complex<double> X(_X.x,_X.y);
	X=sqrt(X);
	return zvar(X.real(),X.imag());
}

cvar sqrt(cvar _X){
	complex<float> X(_X.x,_X.y);
	X=sqrt(X);
	return cvar(X.real(),X.imag());
}

zvar	exp(zvar val){
	zmplx vv(val.x,val.y);
	vv=exp(vv);
	return zvar(vv.real(), vv.imag());
}


//--------------------------------------------------------------
//---Perform Gram-Schmidt orthogonalization to the matrix-------
//--------------------------------------------------------------

template<class vvar>
void Orthonormalize(matrix<vvar> & a, matrix<vvar> & b)
{
	b=matrix<vvar>(1,a.cols()).get();

	for(int i=0 ; i<a.cols() ; i++)
	{
		if(i>0)
		for( int j=0 ; j<i ; j++)
		{
			vvar c1=0;
			for( int k=0 ; k<a.cols() ; k++) { c1=c1+conj(a(k,j))*a(k,i); }
			for( int k=0 ; k<a.cols() ; k++) { a(k,i)=a(k,i)-c1*a(k,j);	  }
		}

		b[i]=0;
		for( int k=0 ; k<a.cols() ; k++){
			b[i]=b[i]+conj(a(k,i))*a(k,i);
		}
		b[i]=sqrt(b[i]);
		for( int k=0 ; k<a.cols() ; k++){
			a(k,i)=a(k,i)/b[i];
		}
	}
}

	

#endif
