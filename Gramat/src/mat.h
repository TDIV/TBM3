#ifndef __MATRIX_BASE_OPERATION__
#define __MATRIX_BASE_OPERATION__

// includes, system
#include <stdio.h>
#include <stdlib.h>

// STD library
#include <string>
#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <complex>
#include <cmath>
#include <memory>

// includes, magma / cuda library
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cublas_v2.h>
#include "magma.h"
#include "magma_lapack.h"


using namespace std;


namespace gmt {
	
	unsigned PRINT_PRECISION=8;
	unsigned NGPU = 1;
	
	#include "mat_complex.h"
	
	float gevd(unsigned int NN, float			*A,	float	*w1, float				*h_A, string &Status);
	float gevd(unsigned int NN, double			*A,	double	*w1, double				*h_A, string &Status);
	float gevd(unsigned int NN, cuFloatComplex	*A,	float	*w1, cuFloatComplex		*h_A, string &Status);
	float gevd(unsigned int NN, cuDoubleComplex	*A,	double	*w1, cuDoubleComplex	*h_A, string &Status);
	#include "mat_evd.h"
	
	float gevdm(unsigned ngpu, unsigned int NN, float			*A,	float	*w1, float				*h_A, string &Status);
	float gevdm(unsigned ngpu, unsigned int NN, double			*A,	double	*w1, double				*h_A, string &Status);
	float gevdm(unsigned ngpu, unsigned int NN, cuFloatComplex	*A,	float	*w1, cuFloatComplex		*h_A, string &Status);
	float gevdm(unsigned ngpu, unsigned int NN, cuDoubleComplex	*A,	double	*w1, cuDoubleComplex	*h_A, string &Status);
	#include "mat_evd_m.h"
	
	#include "mat_fmt.h"


	string fformat(float	val, size_t len);
	string fformat(double	val, size_t len);
	string fformat(cmplx	val, size_t len);
	string fformat(zmplx	val, size_t len);
	
	/* ---------------------------------------------
	--- Class for data storage and transfer. -------
	-----------------------------------------------*/
	template<class T>
	class Arr {
		T * arr;
	public:
		unsigned _cols,_rows;
		 Arr( unsigned N=0): arr(new T[N]){}
		~Arr(){ delete [] arr;}
		T& at(unsigned i){return arr[i];}
		T * get_ptr(){return arr;}
	
	};
	
	/* ---------------------------------------------
	--- Class for operations of data. --------------
	-----------------------------------------------*/
	
	template<class var> class matrix;
	
	typedef matrix<svar>	smat;
	typedef matrix<dvar>	dmat;
	typedef matrix<cvar>	cmat;
	typedef matrix<zvar>	zmat;
	
	enum type{MATRIX, ARRAY};
	enum CPY_TYPE{VAR,PTR};
	
	
	template<class var>
	class matrix {
	private:
		type			_TYPE;
		CPY_TYPE		cpy_type;
		unsigned int	print_len;
		unsigned int	print_flag;
		void			set_col		(unsigned  c)	{mat->_cols=c;}
		void			set_row		(unsigned  r)	{mat->_rows=r;}
	
		auto_ptr<Arr<var> >		mat;
	public:
		//matrix(): mat(new Arr<var>(1)){
		//	mat->_cols=1; mat->_rows=1; zerolize();
		//	print_len = 10;
		//	print_flag=0;
		//	_TYPE=MATRIX;
		//	cpy_type=VAR;
		//}
	
		matrix(unsigned C=1, type tp=MATRIX): mat(new Arr<var>(C*C)){
			mat->_cols=C; mat->_rows=C; zerolize();
			print_len = 10;
			print_flag=0;
			_TYPE=tp;
			cpy_type=VAR;
		}
	
		matrix(unsigned C, unsigned R, type tp=MATRIX): mat(new Arr<var>(C*R)){
			mat->_cols=C; mat->_rows=R; zerolize();
			print_len = 12;
			print_flag=0;
			_TYPE=tp;
			cpy_type=VAR;
		}
	
		matrix(const matrix &  cm):mat(new Arr<var>(cm.cols()*cm.rows())) {
			mat->_cols = cm.cols();
			mat->_rows = cm.rows();
			_TYPE = cm._TYPE;
			cpy_type = VAR;
			print_len = cm.print_len;
			print_flag= cm.print_flag;
			for( unsigned i=0 ; i<size() ; i++) index(i)=cm.const_index(i);
		}
	
		/* ---------------------------------------------
		--- Setting the Type of the data base. ---------
		-----------------------------------------------*/
		type			TYPE	()				{return _TYPE;}
		matrix &		toA		()				{_TYPE=ARRAY;	return *this;}
		matrix &		toM		()				{_TYPE=MATRIX;	return *this;}
		void			set_type(type T)		{_TYPE = T;}
		void			zerolize ()				{ for (unsigned int i=0 ; i<cols()*rows() ; i++) index(i)=0; }
	
		/* ---------------------------------------------
		--- Copy operators -----------------------------
		-----------------------------------------------*/
		matrix&			operator=	( const matrix& );
		matrix&			operator=	(auto_ptr<Arr<var> > );
		CPY_TYPE		get_cpy_type() const	 {return cpy_type;}
		void			set_cpy_type(CPY_TYPE tp){ cpy_type = tp;}
	
		/* ---------------------------------------------
		--- Row(col) access ----------------------------
		-----------------------------------------------*/
		unsigned	cols		()		const	{return mat->_cols;}
		unsigned	rows		()		const	{return mat->_rows;}
		unsigned	size		()		const	{return cols()*rows();}
	
		matrix &			 	col			(unsigned  i);
		matrix &			 	row			(unsigned  i);
		void			col_swap	(unsigned  i,unsigned j);
		void			row_swap	(unsigned  i,unsigned j);
	
		/* ---------------------------------------------
		--- Pointer operation for data transfer. -------
		-----------------------------------------------*/
		var *			get_ptr		(){return mat->get_ptr();}
		auto_ptr<Arr<var> > get     (){return mat;}
	
		/* ---------------------------------------------
		--- Indexing operation -------------------------
		-----------------------------------------------*/
		var &			operator[]	(unsigned i)		{
			if (i>=cols()*rows())	{
				string ErrorStr = "Error, index out of range of matrix element";
				throw ErrorStr;
			}
		return mat->at(i);
		}
	    var &			index		(unsigned i)		{
			if (i>=cols()*rows())			{
				string ErrorStr = "Error, index out of range of matrix element";
				throw ErrorStr;
			}
		return mat->at(i);
		}
	    var &			const_index (unsigned i) const	{
			if (i>=cols()*rows())	{
				string ErrorStr = "Error, index out of range of matrix element";
				throw ErrorStr;
			}
			return mat->at(i);
		}
	
		var &			operator()	(unsigned i, unsigned j){return at(i,j);}
		var &			at			(unsigned , unsigned );
		
		matrix&			operator-()	{
			matrix  *ret = new matrix(cols(),rows());
			ret->set_cpy_type(PTR);

			for (int i=0 ; i< size() ; i++) ret->index(i) = -index(i);

			return *ret;
		}
	
		void			setPrintFlag(unsigned flag){print_flag=flag;}
		void			setPrintLength(unsigned len){print_len = len;}
		unsigned		PrintLength(){return print_len;}
		unsigned		PrintFlag(){return print_flag;}
	
		///* --------------------------------------------- */
		///* - The entry of the gevd --------------------- */
		///* - Which diagonalize symmetric matrix--------- */
		///* --------------------------------------------- */
		///*	By calling one of the following function.  */
		///*	Returns the status of calculations		   */
		string			evd			(smat & E, matrix & V);
		string			evd			(dmat & E, matrix & V);
	};
	
	template<class var>
	string matrix<var>::evd(smat & E, matrix & V) {
		if (cols() != rows()) throw "Object is not a square matrix!";
		E = smat(1,cols()).get();
		V = matrix(cols(),cols()).get();
	
		string Status;
		//float time=gevd(cols(), mat->get_ptr(), E.get_ptr(), V.get_ptr(), Status);
		if( NGPU == 1){
			gevd(cols(), mat->get_ptr(), E.get_ptr(), V.get_ptr(), Status);
		} else {
			gevdm(NGPU,cols(), mat->get_ptr(), E.get_ptr(), V.get_ptr(), Status);
		}
		return Status;
	}
	
	template<class var>
	string matrix<var>::evd(dmat & E, matrix & V) {
		if (cols() != rows()) throw "Object is not a square matrix!";
		E = dmat(1,cols()).get();
		V = matrix(cols(),cols()).get();
	
		string Status;
		//float time=gevd(cols(), mat->get_ptr(), E.get_ptr(), V.get_ptr(), Status);
		if( NGPU == 1){
			gevd(cols(), mat->get_ptr(), E.get_ptr(), V.get_ptr(), Status);
		} else {
			gevdm(NGPU,cols(), mat->get_ptr(), E.get_ptr(), V.get_ptr(), Status);
		}
		return Status;
	}
	
	#include "mat_multiply.h"
	#include "mat_function.h"
	#include "mat_operator.h"
	
} // End of "gramat" namespace


#endif
