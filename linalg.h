/*
	Linear algebra functions and MKL/OpenBLAS Matrix Operation Wrapper by Pak Ki Henry Tsang 
  
  Jan-2021
	
	Compile instructions:
	
	1. Please compile with c++17 onwards enabled (-std=c++17)
	2. for intel compiler with mkl, use compiler flag : -mkl
	3. for intel compiler with openblas, use compiler flags : -I/usr/local/include/openblas/ -lopenblas -lpthread -lgomp -lgfortran
	4. for gcc mkl has to be linked using compiler flags : -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
	
	Compile-time functions are usually much faster than their run-time counterparts
	
	Note that these functions are made to be as simple as possible, which means that they do not check for errors and it is the user's responsibility to check the eligibility of the operation and other possible numerical errors.
	
	TODO :
	1. Conjugate Transpose
	2. non-symmetric eigenvalue problems
	3. eigenvalue-only routines
	4. Matrix Exponential
	5. Permanent
*/


#ifndef _LINALG_H_
#define _LINALG_H_
#endif

#if defined(_MKL_H_) || defined(_MKL_LAPACKE_H_)
#define SCOMPLEX MKL_Complex8
#define DCOMPLEX MKL_Complex16
#define SLAMBDA slamch("S")
#define DLAMBDA dlamch("S")
#elif defined(CBLAS_H) || defined(_LAPACKE_H_)
#define SLAMBDA LAPACKE_slamch( 'S' )
#define DLAMBDA LAPACKE_dlamch( 'S' )
#define SCOMPLEX __complex__ float
#define DCOMPLEX __complex__ double
#endif

#include<complex>
#include <type_traits>

namespace linalg{
	
	//Zero and Identity matrices
	template<class T, size_t m, size_t n> void zero(T* A);
	template<class T> void zero(T* A, const size_t m, const size_t n);
	template<class T, size_t m> void identity(T* A);
	template<class T> void identity(T* A, const size_t m);
	
	//Matrix Addition
	template<class T, size_t m, size_t n> void matadd(const T alpha,T* A, const T beta, T* B,  T* C);
	template<class T> void matadd(const T alpha,T* A, const T beta, T* B, const size_t m, const size_t n,  T* C);
	
	//Matrix Multiplication
	template<class T,size_t m,size_t n,size_t k> void matmul(T* A, T* B, T* C);
	template<class T> void matmul(T* A, T* B, const size_t m, const size_t n, const size_t k, T* C); 
	
	//Matrix Power
	template<class T, size_t n> void matpow(T* A, const size_t m, T* B);
	template<class T> void matpow(T* A, const size_t m, const size_t n, T* B);

	//Matrix Trace
	template<class T, size_t m> T trace(T* A);
	template<class T> T trace(T* A, const size_t m);
	
	//Matrix Transpose
	template<class T> void transpose(T* A,const size_t m, const size_t n, T* B);
	template<class T, size_t m, size_t n> void transpose(T* A, T* B);
	
	void conjugatetranspose(std::complex<float>* A,const size_t m, const size_t n, std::complex<float>* B);
	void conjugatetranspose(std::complex<double>* A,const size_t m, const size_t n, std::complex<double>* B);
	
	template<size_t m, size_t n> void conjugatetranspose(std::complex<float>* A, std::complex<float>* B);
	template<size_t m, size_t n> void conjugatetranspose(std::complex<double>* A, std::complex<double>* B);
	
	//Matrix Minor
	template<class T, size_t m, size_t i, size_t j> void minor(T* A, T* B);
	template<class T> void minor(T* A, const size_t m, const size_t i, const size_t j, T* B);
	
	//Matrix Determinant
	template<class T> T det(T* A, const size_t m);
	template<class T, size_t m> T det(T* A);
	template<class T> T det_2x2(T* A);
	template<class T> T det_2x2(const T a, const T b, const T c,const T d);
	template<class T> T det_3x3(T* A);
	template<class T> T det_4x4(T* A);
	template<class T> T det_leibniz(T* A, const size_t m);
	template<class T> T det_LU(T* A, const size_t m);
	
	//Matrix Inverse
	template<class T> void inv_2x2(T* A, T* B);
	template<class T> void inv_3x3(T* A, T* B);
	template<class T> void inv_4x4(T* A, T* B);
	void inv_LU(double* A, const size_t m, double* B);
	void inv_LU(float* A, const size_t m, float* B);
	void inv_LU(std::complex<double>* A, const size_t m, std::complex<double>* B);
	void inv_LU(std::complex<float>* A, const size_t m, std::complex<float>* B);
	template<class T> void inv(T* A, const size_t m, T* B);
	
	//Eigenvalue problem
	void eigh(std::complex<double>* A , const size_t m, double* eigval, std::complex<double>* eigvec);
	void eigh(std::complex<float>* A , const size_t m, float* eigval, std::complex<float>* eigvec);
	void eigh(double* A, const size_t m, double* eigval, double* eigvec);
	void eigh(float* A, const size_t m, float* eigval, float* eigvec);
	
	//Misc
	constexpr size_t factorial(size_t n);

};




/*******************
	Matrix summation
*******************/

template<class T> void linalg::matadd(const T alpha,T* A, const T beta, T* B, const size_t m, const size_t n,  T* C){
	/*
		General matrix addition for 2 mxn matrices
	*/
	for (size_t i=0;i<m;i++){	for (size_t j=0;j<n;j++){
		C[m*i+j]+=alpha*A[m*i+j]+beta*B[m*i+j];
	}}
}

template<class T, size_t m, size_t n> void linalg::matadd(const T alpha, T* A, const T beta, T* B,  T* C){
	/*
		General matrix addition for 2 mxn matrices
	*/
	for (size_t i=0;i<m;i++){	for (size_t j=0;j<n;j++){
		C[m*i+j]+=alpha*A[m*i+j]+beta*B[m*i+j];
	}}
}



/************************
	Matrix Multiplication
************************/

template<class T> void linalg::matmul(T* A, T* B, const size_t m, const size_t n, const size_t k, T* C){
	/*
		General matrix multiplication for A(mxk) and B(kxn) matrix, output to C(mxn) matrix
	*/
	for (size_t i=0;i<m;i++){	for (size_t j=0;j<n;j++){
		for (size_t l=0;l<k;l++){
			C[m*i+j]+=A[m*i+l]*B[n*l+j];
		}
	}}
}

template<class T,size_t m,size_t n,size_t k> void linalg::matmul(T* A, T* B, T* C){
	/*
		General matrix multiplication for A(mxk) and B(kxn) matrix, output to C(mxn) matrix (Compile time)
	*/
	for (size_t i=0;i<m;i++){	for (size_t j=0;j<n;j++){
		for (size_t l=0;l<k;l++){
			C[m*i+j]+=A[m*i+l]*B[n*l+j];
		}
	}}
}

/***************
	Matrix Power
***************/

template<class T, size_t m, size_t n> void linalg::matpow(T* A, T* B){
	/*
		mxm square matrix to the power n : A^n (Compile time)
	*/
	if constexpr(n==0){
		identity(B,m);
	}
	else if constexpr(n%2==0){ //even n
		T* C = new T[m*m];
		matpow<T,m,n/2>(A, C);
		zero<T, m, m>(B);
		matmul<T, m, m, m>(C, C, B);
		delete [] C;
	}
	else{ //odd n
		T* C = new T[m*m];
		T* D = new T[m*m];
		matpow<T, m,n/2>(A, C);
		zero<T, m, m>(B);
		zero<T, m, m>(D);
		matmul<T, m, m, m>(C, C, D); //D=C^2
		matmul<T, m, m, m>(A, D, B); //B=A*(C^2)
		delete [] C;
		delete [] D;
	}
}

template<class T, size_t n> void linalg::matpow(T* A, const size_t m, T* B){
	/*
		mxm square matrix to the power n : A^n (Compile time)
	*/
	if constexpr(n==0){
		identity(B,m);
	}
	else if constexpr(n%2==0){ //even n
		T* C = new T[m*m];
		matpow<T,n/2>(A, m, C);
		zero<T>(B, m, m);
		matmul<T>(C, C, m, m, m, B);
		delete [] C;
	}
	else{ //odd n
		T* C = new T[m*m];
		T* D = new T[m*m];
		matpow<T,n/2>(A, m, C);
		zero<T>(B, m, m);
		zero<T>(D, m, m);
		matmul<T>(C, C, m, m, m, D); //D=C^2
		matmul<T>(A, D, m, m, m, B); //B=A*(C^2)
		delete [] C;
		delete [] D;
	}
}

template<class T> void linalg::matpow(T* A, const size_t m, const size_t n, T* B){
	/*
		mxm square matrix to the power n : A^n
	*/
	if (n==0){
		identity(B,m);
	}
	else if (n%2==0){ //even n
		T* C = new T[m*m];
		matpow<T>(A, m, n/2, C);
		zero<T>(B, m, m);
		matmul<T>(C, C, m, m, m, B);
		delete [] C;
	}
	else{ //odd n
		T* C = new T[m*m];
		T* D = new T[m*m];
		matpow<T>(A, m, n/2, C);
		zero<T>(B, m, m);
		zero<T>(D, m, m);
		matmul<T>(C, C, m, m, m, D); //D=C^2
		matmul<T>(A, D, m, m, m, B); //B=A*(C^2)
		delete [] C;
		delete [] D;
	}
}

template<class T> void linalg::zero(T* A, const size_t m, const size_t n){
	for (size_t i=0;i<m;i++){ for(size_t j=0;j<n;j++) {
		A[i*m+j]=0;
	}}
}

template<class T, size_t m, size_t n> void linalg::zero(T* A){
	for (size_t i=0;i<m;i++){ for(size_t j=0;j<n;j++) {
		A[i*m+j]=0;
	}}
}

template<class T> void linalg::identity(T* A, const size_t m){

	zero<T>(A, m,m);
	
	for (size_t i=0;i<m;i++){
		A[i*m+i]=1;
	}
}

template<class T, size_t m> void linalg::identity(T* A){

	zero<T,m,m>(A);
	
	for (size_t i=0;i<m;i++){
		A[i*m+i]=1;
	}
}


/***************
	Matrix trace
***************/

template<class T> T linalg::trace(T* A, const size_t m){
	/*
		Computes the trace of a square matrix and return
	*/

	T result=0;
	
	for (size_t i=0;i<m;i++) result+=A[m*i+i];
	
	return result;
}

template<class T, const size_t m> T linalg::trace(T* A){
	/*
		Computes the trace of a square matrix and return
	*/

	T result=0;
	
	for (size_t i=0;i<m;i++) result+=A[m*i+i];
	
	return result;
}
/*******************
	Matrix Transpose
*******************/

template<class T> void linalg::transpose(T* A,const size_t m, const size_t n, T* B){
	/*
		Transpose mxn matrix A, and write into nxm matrix B
	*/
	
	for (size_t i=0;i<m;i++){	for (size_t j=0;j<n;j++){
		B[n*j+i]=A[m*i+j];
	}}	
}

template<class T, size_t m, size_t n> void linalg::transpose(T* A, T* B){
	/*
		Transpose mxn matrix A, and write into nxm matrix B (Compile-time)
	*/
	
	for (size_t i=0;i<m;i++){	for (size_t j=0;j<n;j++){
		B[n*j+i]=A[m*i+j];
	}}	
}

void linalg::conjugatetranspose(std::complex<float>* A,const size_t m, const size_t n, std::complex<float>* B){
	/*
		Transpose mxn matrix A, and write into nxm matrix B
	*/
	
	for (size_t i=0;i<m;i++){	for (size_t j=0;j<n;j++){
		B[n*j+i]=std::conj(A[m*i+j]);
	}}	
}

void linalg::conjugatetranspose(std::complex<double>* A,const size_t m, const size_t n, std::complex<double>* B){
	/*
		Transpose mxn matrix A, and write into nxm matrix B
	*/
	
	for (size_t i=0;i<m;i++){	for (size_t j=0;j<n;j++){
		B[n*j+i]=std::conj(A[m*i+j]);
	}}	
}

template<size_t m, size_t n> void linalg::conjugatetranspose(std::complex<float>* A, std::complex<float>* B){
	/*
		Transpose mxn matrix A, and write into nxm matrix B
	*/
	
	for (size_t i=0;i<m;i++){	for (size_t j=0;j<n;j++){
		B[n*j+i]=std::conj(A[m*i+j]);
	}}	
}

template<size_t m, size_t n> void linalg::conjugatetranspose(std::complex<double>* A, std::complex<double>* B){
	/*
		Transpose mxn matrix A, and write into nxm matrix B
	*/
	
	for (size_t i=0;i<m;i++){	for (size_t j=0;j<n;j++){
		B[n*j+i]=std::conj(A[m*i+j]);
	}}	
}

/***************************************************
	Matrix determinant
	
	for 3x3 and 4x4 the cofactor expansions are used
	for higher dimensions LU factorization is faster
***************************************************/

template<class T> T linalg::det(T* A, const size_t m){
	/*
		Computes matrix determinant of mxm square matrix
	*/
	switch(m){
		case 1:
			return 0;
			break;
		case 2:
			return det_2x2<T>(A);
			break;
		case 3:
			return det_3x3<T>(A);
			break;
		case 4:
			return det_4x4<T>(A);
			break;
		default:
			return det_LU<T>(A,m);
			break;
	}
}

template<class T, size_t m> T linalg::det(T* A){
	/*
		Computes matrix determinant of mxm square matrix
	*/
	if constexpr(m==1)
		return 0;
	else if constexpr(m==2)
		return det_2x2<T>(A);
	else if constexpr(m==3)
		return det_3x3<T>(A);
	else if constexpr(m==4)
		return det_4x4<T>(A);
	else
		return det_LU<T>(A,m);
}

template<class T> T linalg::det_2x2(T* A){
	return A[0]*A[2]-A[1]*A[3];
}

template<class T> T linalg::det_2x2(const T a, const T b, const T c,const T d){
	return a*d-b*c;
}

template<class T> T linalg::det_3x3(T* A){ //Cofactor expansion
	return A[0]*(A[4]*A[8]-A[5]*A[7])+A[1]*(A[5]*A[6]-A[3]*A[8])+A[2]*(A[3]*A[7]-A[4]*A[6]);
}

template<class T> T linalg::det_4x4(T* A){ //Cofactor expansion

	T result=0;
	
	T* M11 = new T[9]{A[5],A[6],A[7],A[9],A[10],A[11],A[13],A[14],A[15]};
	T* M12 = new T[9]{A[4],A[6],A[7],A[8],A[10],A[11],A[12],A[14],A[15]};
	T* M13 = new T[9]{A[4],A[5],A[7],A[8],A[9],A[11],A[12],A[13],A[15]};
	T* M14 = new T[9]{A[4],A[5],A[6],A[8],A[9],A[10],A[12],A[13],A[14]};
	
	result=A[0]*det_3x3<T>(M11)-A[1]*det_3x3<T>(M12)+A[2]*det_3x3<T>(M13)-A[3]*det_3x3<T>(M14);
	
	delete [] M11;
	delete [] M12;
	delete [] M13;
	delete [] M14;
	
	return result;
}

template<class T> T linalg::det_LU(T* A, const size_t m){
	T* B = new T[m*m];
	for (size_t i=0;i<m*m;i++) B[i]=A[i]; //copy A to B
	
	int* ipiv = new int[m]; 
	
	if constexpr(std::is_same<T,float>::value)
		LAPACKE_sgetrf(LAPACK_ROW_MAJOR, m , m , B , m , ipiv );
	else if constexpr(std::is_same<T,double>::value)
		LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m , m , B , m , ipiv );
	else if constexpr(std::is_same<T,std::complex<float>>::value)
		//LAPACKE_cgetrf(LAPACK_ROW_MAJOR, m , m , reinterpret_cast <__complex__ float*>(B) , m , ipiv );
		LAPACKE_cgetrf(LAPACK_ROW_MAJOR, m , m , B , m , ipiv );
	else if constexpr(std::is_same<T,std::complex<double>>::value)
		//LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m , m , reinterpret_cast <__complex__ double*>(B) , m , ipiv );
		LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m , m , B , m , ipiv );
	else //error
		abort();
	
	T result=1;
	for (size_t i=0;i<m;i++){
		if ( i!=(ipiv[i]-1) )
			result*=-B[i*m+i];
		else
			result*=B[i*m+i];
	}
	delete [] ipiv;
	return result;
}
	

template<class T> T linalg::det_leibniz(T* A, const size_t m){
	/*
		Not yet implemented
	*/
	T result=0;
	return result;
}


//Matrix Inversion
template<class T> void linalg::inv_2x2(T* A, T* B){
	const T det=det_2x2<T>(A);
	B[0]=A[3]/det;
	B[1]=-A[1]/det;
	B[2]=-A[2]/det;
	B[3]=A[0]/det;
}


template<class T> void linalg::inv_3x3(T* A, T* B){
	const T det=det_3x3<T>(A);
	B[0]=det_2x2<T>(A[4],A[5],A[7],A[8])/det;
	B[1]=det_2x2<T>(A[2],A[1],A[8],A[7])/det;
	B[2]=det_2x2<T>(A[1],A[2],A[4],A[5])/det;
	B[3]=det_2x2<T>(A[5],A[3],A[8],A[6])/det;
	B[4]=det_2x2<T>(A[0],A[2],A[6],A[8])/det;
	B[5]=det_2x2<T>(A[2],A[0],A[5],A[3])/det;
	B[6]=det_2x2<T>(A[3],A[4],A[6],A[7])/det;
	B[7]=det_2x2<T>(A[1],A[0],A[7],A[6])/det;
	B[8]=det_2x2<T>(A[0],A[1],A[3],A[4])/det;
}

template<class T> void linalg::inv_4x4(T* A, T* B){
	const T detA = det<T,4>(A);
	printf("detA=%f\n\n",detA);
	T* Identity = new T[16];
	T* A2 = new T[16]();
	T* A3 = new T[16]();
	
	zero<T,4,4>(B);//set B to zero
	//zero<T,4,4>(A2);//set A2 to zero
	//zero<T,4,4>(A3);//set A3 to zero
	
	matmul<T,4,4,4>(A,A,A2); //A^2
	matmul<T,4,4,4>(A,A2,A3); //A^3
	
	const T traceA = trace<T,4>(A);
	printf("trA=%f\n\n",traceA);
	const T traceA2 = trace<T,4>(A2);
	printf("trA2=%f\n\n",traceA2);
	const T traceA3 = trace<T,4>(A3);
	printf("trA3=%f\n\n",traceA3);
	
	identity<T,4>(Identity);//identity matrix
	
	matadd<T,4,4>(traceA/detA,A2,-1/detA,A3,B); //B += tr(A)*A^2 -A^3
	matadd<T,4,4>((traceA*traceA*traceA-3*traceA*traceA2+2*traceA3)/6/detA,Identity,-(traceA*traceA-traceA2)/2/detA,A,B);
	
	
	delete [] Identity;
	delete [] A2;
	delete [] A3;
}

template<class T> void linalg::inv_LU(T* A, const size_t m, T* B){
	/*
		Matrix Inversion using LU decomposition (Lapack)
	*/
	
	//Allocate pivot memory
	int* ipiv = new int[m]; 
	
	//Copy the matrix to be inverted, as it has to go through LU factorization and then inverse
	for (size_t i=0;i<m*m;i++) B[i]=A[i];
	
	if constexpr(std::is_same<T,float>::value){
		int LU_FACTOR = LAPACKE_sgetrf (LAPACK_ROW_MAJOR , m, m , B , m , ipiv ); //Use LAPACK routines to get LU Factorization of the input matrix
		int INVERSE = LAPACKE_sgetri (LAPACK_ROW_MAJOR , m, B , m , ipiv ); //Feed the above results to get the inverse
	}
	else if constexpr(std::is_same<T,double>::value){
		int LU_FACTOR = LAPACKE_dgetrf (LAPACK_ROW_MAJOR , m, m , B , m , ipiv ); 
		int INVERSE = LAPACKE_dgetri (LAPACK_ROW_MAJOR , m, B , m , ipiv );
	} 
	else if constexpr(std::is_same<T,std::complex<float>>::value){
		int LU_FACTOR = LAPACKE_cgetrf (LAPACK_ROW_MAJOR , m, m , reinterpret_cast <SCOMPLEX*>(B) , m , ipiv ); 
		int INVERSE = LAPACKE_cgetri (LAPACK_ROW_MAJOR , m, reinterpret_cast <SCOMPLEX*>(B) , m , ipiv ); 
		//int LU_FACTOR = LAPACKE_cgetrf (LAPACK_ROW_MAJOR , m, m , B , m , ipiv ); 
		//int INVERSE = LAPACKE_cgetri (LAPACK_ROW_MAJOR , m, B , m , ipiv ); 
	}
	else if constexpr(std::is_same<T,std::complex<double>>::value){
		int LU_FACTOR = LAPACKE_zgetrf (LAPACK_ROW_MAJOR , m, m , reinterpret_cast <DCOMPLEX*>(B) , m , ipiv ); 
		int INVERSE = LAPACKE_zgetri (LAPACK_ROW_MAJOR , m, reinterpret_cast <DCOMPLEX*>(B) , m , ipiv ); 
		//int LU_FACTOR = LAPACKE_zgetrf (LAPACK_ROW_MAJOR , m, m , B , m , ipiv ); 
		//int INVERSE = LAPACKE_zgetri (LAPACK_ROW_MAJOR , m, B , m , ipiv ); 
	}
	else //error
		abort();
		
	delete [] ipiv; 
}


template<class T> void linalg::inv(T* A, const size_t m, T* B){
	switch(m){
		case 1:
			B[0]=1.0/A[0];
			break;
		case 2:
			inv_2x2<T>(A,B);
			break;
		case 3:
			inv_3x3<T>(A,B);
			break;
		case 4:
			inv_4x4<T>(A,B);
			break;
		default:
			inv_LU<T>(A, m, B);
			break;
	}
}

template<class T, size_t m> void linalg::inv(T* A, T* B){
	if constexpr(m==1)
		B[0]=1.0/A[0];
	else if constexpr(m==2)
		inv_2x2<T>(A,B);
	else if constexpr(m==3)
		inv_3x3<T>(A,B);
	else if constexpr(m==4)
		inv_4x4<T>(A,B);
	else
		inv_LU<T>(A, m, B);
}

/*
	Eigenvalue problem
*/

void linalg::eigh(std::complex<double>* A , const size_t m, double* eigval, std::complex<double>* eigvec){
	/*
		Computes eigenvalue and eigenvectors of a Hermitian mxm matrix A
		output is mxm eigvec and m eigval
		
		eigenvectors are placed in eigvec in the column form (v1 v2 ...) and this corresponse to the matrix U such that U^dagger * A * U is diagonal
	*/
	//int* isuppz = new int[2*m];
	//int eigval_found;
	//int SYEVR = LAPACKE_zheevr (LAPACK_ROW_MAJOR,'V', 'A' , 'U', m, reinterpret_cast <DCOMPLEX*>(A), m, 0, 0, 0, 0, DLAMBDA, &eigval_found, eigval,reinterpret_cast <DCOMPLEX*>(eigvec), m, isuppz);
	//delete [] isuppz;
	for (size_t i=0;i<m*m;i++) eigvec[i]=A[i];//Copy A to eigvec
	int ZHEEVD = LAPACKE_zheevd(LAPACK_ROW_MAJOR,'V',  'L', m, reinterpret_cast <DCOMPLEX*>(eigvec), m, eigval );
}


void linalg::eigh(std::complex<float>* A , const size_t m, float* eigval, std::complex<float>* eigvec){
	//int* isuppz = new int[2*m];
	//int eigval_found;
	//int SYEVR = LAPACKE_cheevr (LAPACK_ROW_MAJOR,'V', 'A' , 'U', m, reinterpret_cast <SCOMPLEX*>(A), m, 0, 0, 0, 0, SLAMBDA, &eigval_found, eigval, reinterpret_cast <SCOMPLEX*>(eigvec), m, isuppz);
	//delete [] isuppz;
	for (size_t i=0;i<m*m;i++) eigvec[i]=A[i];//Copy A to eigvec
	int CHEEVD = LAPACKE_cheevd(LAPACK_ROW_MAJOR,'V',  'L', m, reinterpret_cast <SCOMPLEX*>(eigvec), m, eigval );
}

void linalg::eigh(double* A, const size_t m, double* eigval, double* eigvec){
	//int* isuppz = new int[2*m];
	//int eigval_found;
	//int SYEVR = LAPACKE_dsyevr (LAPACK_ROW_MAJOR,'V', 'A' , 'U', m, A, m, 0, 0, 0, 0, DLAMBDA, &eigval_found, eigval, eigvec, m, isuppz);
	//delete [] isuppz;
	for (size_t i=0;i<m*m;i++) eigvec[i]=A[i];//Copy A to eigvec
	int DSYEVD = LAPACKE_dsyevd(LAPACK_ROW_MAJOR,'V',  'L', m, eigvec, m, eigval );
}


void linalg::eigh(float* A, const size_t m, float* eigval, float* eigvec){
	//int* isuppz = new int[2*m];
	//int eigval_found;
	//int SYEVR = LAPACKE_ssyevr (LAPACK_ROW_MAJOR,'V', 'A' , 'U', m, A, m, 0, 0, 0, 0, SLAMBDA, &eigval_found, eigval, eigvec, m, isuppz);
	//delete [] isuppz;
	for (size_t i=0;i<m*m;i++) eigvec[i]=A[i];//Copy A to eigvec
	int SSYEVD = LAPACKE_ssyevd(LAPACK_ROW_MAJOR,'V',  'L', m, eigvec, m, eigval );
}




/*
	Misc Functions
*/
constexpr size_t linalg::factorial(size_t n) {
	//Taken from https://en.cppreference.com/w/cpp/language/constexpr
  return n==0 ? 1 : (n * factorial(n - 1));
}


