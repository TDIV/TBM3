#ifndef __MATRIX_OPERATOR_OPERATION__
#define __MATRIX_OPERATOR_OPERATION__

void err(string s){
	cout<<s<<endl;
	throw s;
}
template<class var>
matrix<var> & matrix<var>::operator=(const matrix &rhs) {
	
	{	//// Copy the entire data base to another instance.

        auto_ptr<Arr<var> > _mat(new Arr<var>(rhs.cols()*rhs.rows()));
        _mat->_cols = rhs.cols();
        _mat->_rows = rhs.rows();

        for (unsigned i=0 ; i < rhs.size(); i++){ _mat->at(i) = rhs.const_index(i); }
        mat = _mat;

	}
	return *this;
}

template<class var>
matrix<var> & matrix<var>::operator=(auto_ptr<Arr<var> > rhs){
	set_col(rhs->_cols);
	set_row(rhs->_rows);
	mat = rhs;
	return *this;
}
/* -------------------------------------------------------
-- Adding operation --------------------------------------
---------------------------------------------------------*/
template<class var>
matrix<var>& operator+(matrix<var> & lhs, matrix<var> & rhs){
	if (lhs.cols() != rhs.cols()) err( "Adding in different _cols");
	if (lhs.rows() != rhs.rows()) err( "Adding in different _rows");

	matrix<var>  *ret = new matrix<var>(lhs.cols(),lhs.rows());
	ret->set_cpy_type(PTR);

	for (int i=0 ; i< lhs.cols()*lhs.rows() ; i++) ret->index(i) = lhs[i]+rhs[i];

	return *ret;
}

template<class var, class T>
matrix<var>& operator+(matrix<var> & lhs, T rhs){
	matrix<var>  *ret = new matrix<var>(lhs.cols(),lhs.rows());
	ret->set_cpy_type(PTR);

	for (int i=0 ; i< lhs.cols()*lhs.rows() ; i++) ret->index(i) = lhs[i]+rhs;

	return *ret;
}

template<class var, class T>
matrix<var>& operator+(T lhs, matrix<var> & rhs){
	matrix<var>  *ret = new matrix<var>(rhs.cols(),rhs.rows());
	ret->set_cpy_type(PTR);

	for (int i=0 ; i< rhs.cols()*rhs.rows() ; i++) ret->index(i) = lhs+rhs[i];

	return *ret;
}

/* ------------------------------------------------------
-- Minus operation --------------------------------------
---------------------------------------------------------*/
template<class var>
matrix<var>& operator-(matrix<var> & lhs, matrix<var> & rhs){
	if (lhs.cols() != rhs.cols()) err( "Adding in different _cols");
	if (lhs.rows() != rhs.rows()) err( "Adding in different _rows");

	matrix<var>  *ret = new matrix<var>(lhs.cols(),lhs.rows());
	ret->set_cpy_type(PTR);

	for (int i=0 ; i< lhs.cols()*lhs.rows() ; i++) ret->index(i) = lhs[i]-rhs[i];

	return *ret;
}

template<class var, class T>
matrix<var>& operator-(matrix<var> & lhs, T rhs){
	matrix<var>  *ret = new matrix<var>(lhs.cols(),lhs.rows());
	ret->set_cpy_type(PTR);

	for (int i=0 ; i< lhs.cols()*lhs.rows() ; i++) ret->index(i) = lhs[i]+rhs;

	return *ret;
}

template<class var, class T>
matrix<var>& operator-(T lhs, matrix<var> & rhs){
	matrix<var>  *ret = new matrix<var>(rhs.cols(),rhs.rows());
	ret->set_cpy_type(PTR);

	for (int i=0 ; i< rhs.cols()*rhs.rows() ; i++) ret->index(i) = lhs+rhs[i];

	return *ret;
}

/* -------------------------------------------------------
-- Multiply operation ------------------------------------
---------------------------------------------------------*/
template<class var>
matrix<var>& operator*(matrix<var> & lhs, matrix<var> & rhs){

	matrix<var>  *ret = new matrix<var>(lhs.cols(),rhs.rows());
	ret->set_cpy_type(PTR);

	if		(lhs.TYPE() == ARRAY	and rhs.TYPE() == ARRAY){
		if (lhs.cols() != rhs.cols() and lhs.rows() != rhs.cols()) err( "Array dimension does not match!");
		for (unsigned int i=0 ; i< lhs.cols()*lhs.rows() ; i++){ ret->index(i) = lhs[i]*rhs[i]; }
	}
	else if	(lhs.TYPE() == MATRIX	and rhs.TYPE() == MATRIX){
		if (lhs.cols() != rhs.rows()) err( "Matrix does not match in k-index!");
		magma_mul(lhs, rhs, *ret); // Call Magma_multipy function
	}
	else if	(lhs.TYPE() == ARRAY	and rhs.TYPE() == MATRIX)	err( "Data type does not match: lhs->Array & rhs->Matrix");
	else if	(lhs.TYPE() == MATRIX	and rhs.TYPE() == ARRAY)	err( "Data type does not match: lhs->Matrix & rhs->ARRAY");
	
	ret->set_type(lhs.TYPE());

	return *ret;
}


template<class var>
matrix<var>& operator*(matrix<var> & lhs, double rhs){
	matrix<var>  *ret = new matrix<var>(lhs.cols(),lhs.rows());
	ret->set_cpy_type(PTR);

	for (unsigned int i=0 ; i< lhs.cols()*lhs.rows() ; i++) ret->index(i) = lhs[i]*rhs;

	ret->set_type(lhs.TYPE());

	return *ret;
}

template<class var>
matrix<var>& operator*(double lhs, matrix<var> & rhs){
	matrix<var>  *ret = new matrix<var>(rhs.cols(),rhs.rows());
	ret->set_cpy_type(PTR);

	for (int i=0 ; i< rhs.cols()*rhs.rows() ; i++) ret->index(i) = lhs*rhs[i];

	ret->set_type(rhs.TYPE());

	return *ret;
}


//
///* ------------------------------------------------------- */
///* -Index accessing of the operator ---------------------- */
///* ------------------------------------------------------- */
template<class var>
var & matrix<var>::at(unsigned int i, unsigned int j){

	/* Indexing operation for the 2D matrix */
	if (i>=cols()) {throw "Index i out of range.";}
	if (j>=rows()) {throw "Index j out of range.";}
	return index(j*cols()+i);
}
//
///* ------------------------------------------------------- */
///* -Column access----------------------------------------- */
///* ------------------------------------------------------- */
template<class var>
matrix<var> & 	matrix<var>::col(unsigned int index)
{	if (index>=cols()) {throw "Index out of range.";}

	matrix<var> * ret = new matrix<var>(cols(),1);
	ret->set_cpy_type(PTR);

	for( unsigned int i=0 ; i<cols(); i++) ret->at(i,0) = at(i,index);

	return *ret;
}

///* ------------------------------------------------------- */
///* -Row access-------------------------------------------- */
///* ------------------------------------------------------- */
template<class var>
matrix<var> & 	matrix<var>::row(unsigned int index)
{	if (index>=rows()) {throw "Index out of range.";}

	matrix<var> * ret = new matrix<var>(1,rows());
	ret->set_cpy_type(PTR);

	for( unsigned int i=0 ; i<rows(); i++) ret->at(0,i) = at(index,i);

	return *ret;
}

///* ------------------------------------------------------- */
///* -Swap two columns-------------------------------------- */
///* ------------------------------------------------------- */
template<class var>
void	matrix<var>::col_swap	(unsigned  i,unsigned j){
	if (i >=rows()) {throw "Index out of range.";}
	if (j >=rows()) {throw "Index out of range.";}

	for( unsigned index=0 ; index<cols(); index++){
		var buffer=at(i,index);
		at(i,index) = at(j,index);
		at(j,index) = buffer;
	}
}
//
///* ------------------------------------------------------- */
///* -Swap two rows----------------------------------------- */
///* ------------------------------------------------------- */
template<class var>
void	matrix<var>::row_swap	(unsigned  i,unsigned j){
	if (i >=cols()) {throw "Index out of range.";}
	if (j >=cols()) {throw "Index out of range.";}

	for( unsigned index=0 ; index<rows(); index++){
		var buffer=at(index,i);
		at(index,i) = at(index,j);
		at(index,j) = buffer;
	}
}
//
///* ------------------------------------------------------- */
///* -Output of the matrix---------------------------------- */
///* ------------------------------------------------------- */
template<class var>
ostream & operator<< (ostream & out, matrix<var> & m){
	
	// Determine the max length of each component
	int max_len=m.PrintLength();
	for (unsigned int i=0 ; i< m.cols() ; i++)
	for (unsigned int j=0 ; j< m.rows() ; j++){
		string ss = tostr(m(i,j));
		if (max_len < ss.length() )
			max_len = ss.length();
	}
	

	for (unsigned int i=0 ; i< m.cols() ; i++){
		if ( i>0 )out<<endl;
		for (unsigned int j=0 ; j< m.rows() ; j++){
			if (i==0 and j==0) out<<"[[ ";
			else if (j==0) out<<" [ ";
			else out<<" ";
			
			out<<fformat(m(i,j), max_len+1);
			
			if (i==m.cols()-1 and j==m.rows()-1) out<<"]]";
			else if (j==m.rows()-1) out<<"] ";
		}
	}
	return out;
}

	

#endif
