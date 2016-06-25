#ifndef __MAT_MUL_H__
#define __MAT_MUL_H__

//namespace gmt

// Calling magma_sgemm
double magma_smul(smat & A, smat & B, smat & C){

	if( A.size() ==0 or B.size() ==0 )	return -1;
	if( A.rows() != B.cols() )			return -1;

	C = smat(A.cols(), B.rows()).get();

	magma_init();
	magma_int_t m = A.cols();
	magma_int_t k = A.rows();
	magma_int_t n = B.rows();
	magma_int_t mk = m*k;
	magma_int_t kn = k*n;
	magma_int_t mn = m*n;
	float	*a=A.get_ptr(),
			*b=B.get_ptr(),
			*c=C.get_ptr(),
			*d_a, *d_b, *d_c;

	float alpha = MAGMA_S_MAKE(1.0, 0.0);
	float beta	= MAGMA_S_MAKE(0.0, 0.0);

	// device mem. for a, b, c
	magma_smalloc( &d_a, mk );
	magma_smalloc( &d_b, kn );
	magma_smalloc( &d_c, mn );

	// copy data from host to device
	magma_ssetmatrix( m, k, a, m, d_a, m );
	magma_ssetmatrix( k, n, b, k, d_b, k );
	magma_ssetmatrix( m, n, c, m, d_c, m );

	real_Double_t calc_time = magma_wtime();
	magma_sgemm(MagmaNoTrans,MagmaNoTrans,m,n,k,alpha,d_a,m,d_b,k, beta,d_c,m);
	calc_time = calc_time-magma_wtime();

	magma_sgetmatrix( m, n, d_c, m, c, m );

	magma_free(d_a);
	magma_free(d_b);
	magma_free(d_c);
	magma_finalize ();

	return calc_time;
}

// Calling magma_dgemm 
double magma_dmul(dmat & A, dmat & B, dmat & C){

	if( A.size() ==0 or B.size() ==0 )	return -1;
	if( A.rows() != B.cols() )			return -1;

	C = dmat(A.cols(), B.rows()).get();

	magma_init();
	magma_int_t m = A.cols();
	magma_int_t k = A.rows();
	magma_int_t n = B.rows();
	magma_int_t mk = m*k;
	magma_int_t kn = k*n;
	magma_int_t mn = m*n;
	double	*a=A.get_ptr(),
			*b=B.get_ptr(),
			*c=C.get_ptr(),
			*d_a, *d_b, *d_c;

	double	alpha	= MAGMA_D_MAKE(1.0, 0.0);
	double	beta	= MAGMA_D_MAKE(0.0, 0.0);

	// device mem. for a, b, c
	magma_dmalloc( &d_a, mk );
	magma_dmalloc( &d_b, kn );
	magma_dmalloc( &d_c, mn );

	// copy data from host to device
	magma_dsetmatrix( m, k, a, m, d_a, m );
	magma_dsetmatrix( k, n, b, k, d_b, k );
	magma_dsetmatrix( m, n, c, m, d_c, m );

	real_Double_t calc_time = magma_wtime();
	magma_dgemm(MagmaNoTrans,MagmaNoTrans,m,n,k,alpha,d_a,m,d_b,k, beta,d_c,m);
	calc_time = calc_time-magma_wtime();

	magma_dgetmatrix( m, n, d_c, m, c, m );

	magma_free(d_a);
	magma_free(d_b);
	magma_free(d_c);
	magma_finalize ();

	return calc_time;
}

// Calling magma_cgemm 
double magma_cmul(cmat & A, cmat & B, cmat & C){

	if( A.size() ==0 or B.size() ==0 )	return -1;
	if( A.rows() != B.cols() )			return -1;

	C = cmat(A.cols(), B.rows()).get();

	magma_init();
	magma_int_t m = A.cols();
	magma_int_t k = A.rows();
	magma_int_t n = B.rows();
	magma_int_t mk = m*k;
	magma_int_t kn = k*n;
	magma_int_t mn = m*n;
	cuFloatComplex	*a=A.get_ptr(),
					*b=B.get_ptr(),
					*c=C.get_ptr(),
					*d_a, *d_b, *d_c;

	cuFloatComplex 	alpha	= MAGMA_C_MAKE(1.0, 0.0);
	cuFloatComplex	beta	= MAGMA_C_MAKE(0.0, 0.0);

	// device mem. for a, b, c
	magma_cmalloc( &d_a, mk );
	magma_cmalloc( &d_b, kn );
	magma_cmalloc( &d_c, mn );

	// copy data from host to device
	magma_csetmatrix( m, k, a, m, d_a, m );
	magma_csetmatrix( k, n, b, k, d_b, k );
	magma_csetmatrix( m, n, c, m, d_c, m );

	real_Double_t calc_time = magma_wtime();
	magma_cgemm(MagmaNoTrans,MagmaNoTrans,m,n,k,alpha,d_a,m,d_b,k, beta,d_c,m);
	calc_time = calc_time-magma_wtime();

	magma_cgetmatrix( m, n, d_c, m, c, m );

	magma_free(d_a);
	magma_free(d_b);
	magma_free(d_c);
	magma_finalize ();

	return calc_time;
}

// Calling magma_zgemm 
double magma_zmul(zmat & A, zmat & B, zmat & C){

	if( A.size() ==0 or B.size() ==0 )	return -1;
	if( A.rows() != B.cols() )			return -1;

	C = zmat(A.cols(), B.rows()).get();

	magma_init();
	magma_int_t m = A.cols();
	magma_int_t k = A.rows();
	magma_int_t n = B.rows();
	magma_int_t mk = m*k;
	magma_int_t kn = k*n;
	magma_int_t mn = m*n;
	cuDoubleComplex	*a=A.get_ptr(),
						*b=B.get_ptr(),
						*c=C.get_ptr(),
						*d_a, *d_b, *d_c;

	cuDoubleComplex 	alpha=MAGMA_Z_MAKE(1.0, 0.0);
	cuDoubleComplex		beta= MAGMA_Z_MAKE(0.0, 0.0);

	// device mem. for a, b, c
	magma_zmalloc( &d_a, mk );
	magma_zmalloc( &d_b, kn );
	magma_zmalloc( &d_c, mn );

	// copy data from host to device
	magma_zsetmatrix( m, k, a, m, d_a, m );
	magma_zsetmatrix( k, n, b, k, d_b, k );
	magma_zsetmatrix( m, n, c, m, d_c, m );

	real_Double_t calc_time = magma_wtime();
	magma_zgemm(MagmaNoTrans,MagmaNoTrans,m,n,k,alpha,d_a,m,d_b,k, beta,d_c,m);
	calc_time = calc_time-magma_wtime();

	magma_zgetmatrix( m, n, d_c, m, c, m );

	magma_free(d_a);
	magma_free(d_b);
	magma_free(d_c);
	magma_finalize ();

	return calc_time;
}



#endif
