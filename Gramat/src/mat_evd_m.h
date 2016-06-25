//namespace gmt

//-----------------------------------
// N  : dimension of the matrix,
// A  : The matrix to be diagonalized
// w1 : The eigenvalues
// h_A: The eigenvectors
//-----------------------------------

//-----------------------------------
//magma_ssyevd
//------------------------------------
float gevdm(unsigned ngpu, unsigned int N, float *A, float *w1, float *h_A, string &Status) {

	// Start counting the calculation time
	real_Double_t calc_time = magma_wtime();


	magma_init ();
	magma_int_t n=N;
	float  *h_work;	// single precision
	magma_int_t lwork, *iwork, liwork, info;


	lapackf77_slacpy( MagmaUpperLowerStr, &n, &n, A, &n, h_A, &n ); // single precision, copy A->h_A

	// Query for workspace sizes
	float  aux_work [1]; // single precision
	magma_int_t aux_iwork[1];

	magma_ssyevd_m(ngpu,MagmaVec,MagmaLower,n,h_A,n,w1,aux_work ,-1, aux_iwork ,-1,&info ); //Single prcision query for dimension

	lwork = (magma_int_t) aux_work[0];
	liwork = aux_iwork[0];


	iwork=(magma_int_t*)malloc(liwork*sizeof(magma_int_t));
	magma_smalloc_cpu(&h_work,lwork); //memory query


	// Perform eigen-value problem
	magma_ssyevd_m(ngpu,MagmaVec,MagmaLower,n,h_A,n,w1,h_work,lwork, iwork,liwork,&info);

	free(h_work);
	magma_finalize();


	if (info != 0)  Status=magma_strerror( info );
	else			Status="Success";

	// End counting the calculation time
	calc_time = magma_wtime() - calc_time;
	return calc_time;
}

//-----------------------------------
//magma_dsyevd
//-----------------------------------
float gevdm(unsigned ngpu, unsigned int N, double *A, double *w1, double *h_A, string &Status) {

	// Start counting the calculation time
	real_Double_t calc_time = magma_wtime();

	magma_init ();
	magma_int_t n=N;
	double *h_work;	// double precision
	magma_int_t lwork, *iwork, liwork, info;


	lapackf77_dlacpy( MagmaUpperLowerStr, &n, &n, A, &n, h_A, &n ); // double precision, copy A->h_A

	// Query for workspace sizes
	double aux_work [1]; // double precision
	magma_int_t aux_iwork[1];

	magma_dsyevd_m(ngpu, MagmaVec,MagmaLower,n,h_A,n,w1,aux_work ,-1, aux_iwork ,-1,&info ); //Double prcision query for dimension

	lwork = (magma_int_t) aux_work[0];
	liwork = aux_iwork[0];


	iwork=(magma_int_t*)malloc(liwork*sizeof(magma_int_t));
	magma_dmalloc_cpu(&h_work,lwork); //memory query


	// Perform eigen-value problem
	magma_dsyevd_m(ngpu,MagmaVec,MagmaLower,n,h_A,n,w1,h_work,lwork, iwork,liwork,&info);

	free(h_work);
	magma_finalize();


	if (info != 0)  Status=magma_strerror( info );
	else			Status="Success";

	// End counting the calculation time
	calc_time = magma_wtime() - calc_time;
	return calc_time;
}

//-----------------------------------
//magma_cheevd
//-----------------------------------
float gevdm(unsigned ngpu,unsigned int N, cuFloatComplex *A, float *w1, cuFloatComplex *h_A, string &Status) {

	// Start counting the calculation time
	real_Double_t calc_time = magma_wtime();

	magma_init ();
	magma_int_t n=N;
	cuFloatComplex *h_work;	// single complex precision
	magma_int_t lwork, *iwork, liwork, lrwork, info;
	float aux_rwork[1],*rwork;

	lapackf77_clacpy( MagmaUpperLowerStr, &n, &n, A, &n, h_A, &n ); // double precision, copy A->h_A

	// Query for workspace sizes
	cuFloatComplex aux_work [1]; // double precision
	magma_int_t aux_iwork[1];

	magma_cheevd_m(ngpu,MagmaVec, MagmaLower, n,h_A,n,w1	,aux_work ,-1 ,aux_rwork,-1 ,aux_iwork,-1 ,&info ); //single complex prcision query for dimension

    lwork  = (magma_int_t) MAGMA_C_REAL( aux_work[0] );
    lrwork = (magma_int_t) aux_rwork[0];
    liwork = aux_iwork[0];

	iwork=(magma_int_t*)malloc(liwork*sizeof(magma_int_t));
	magma_cmalloc_cpu(&h_work,lwork); //memory query
	magma_smalloc_cpu(&rwork,lrwork); //memory query

	// Perform eigen-value problem
	magma_cheevd_m(ngpu,MagmaVec, MagmaLower, n, h_A, n, w1, h_work, lwork, rwork, lrwork, iwork, liwork, &info );

	free(h_work);
	free(rwork);
	magma_finalize();

	if (info != 0)  Status=magma_strerror( info );
	else			Status="Success";

	//// End counting the calculation time
	calc_time = magma_wtime() - calc_time;
	return calc_time;
}

//-----------------------------------
//magma_zheevd
//-----------------------------------
float gevdm(unsigned ngpu, unsigned int N, cuDoubleComplex *A, double *w1, cuDoubleComplex *h_A, string &Status) {

	// Start counting the calculation time
	real_Double_t calc_time = magma_wtime();

	magma_init ();
	magma_int_t n=N;
	cuDoubleComplex *h_work;	// single complex precision
	magma_int_t lwork, *iwork, liwork, lrwork, info;
	double	aux_rwork[1],*rwork;

	lapackf77_zlacpy( MagmaUpperLowerStr, &n, &n, A, &n, h_A, &n ); // double precision, copy A->h_A

	// Query for workspace sizes
	cuDoubleComplex aux_work [1]; // double precision
	magma_int_t aux_iwork[1];

	magma_zheevd_m(ngpu,MagmaVec, MagmaLower, n,h_A,n,w1	,aux_work ,-1 ,aux_rwork,-1 ,aux_iwork,-1 ,&info ); //single complex prcision query for dimension

	lwork  = (magma_int_t) MAGMA_C_REAL( aux_work[0] );
	lrwork = (magma_int_t) aux_rwork[0];
	liwork = aux_iwork[0];

	iwork=(magma_int_t*)malloc(liwork*sizeof(magma_int_t));
	magma_zmalloc_cpu(&h_work,lwork); //memory query
	magma_dmalloc_cpu(&rwork,lrwork); //memory query

	// Perform eigen-value problem
	magma_zheevd_m(ngpu,MagmaVec, MagmaLower, n, h_A, n, w1, h_work, lwork, rwork, lrwork, iwork, liwork, &info );

	free(h_work);
	free(rwork);
	magma_finalize();

	if (info != 0)  Status=magma_strerror( info );
	else			Status="Success";

	//// End counting the calculation time
	calc_time = magma_wtime() - calc_time;
	return calc_time;
}
