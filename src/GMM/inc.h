#ifndef _NR3_H_
#define _NR3_H_

//#define _CHECKBOUNDS_ 1
//#define _USESTDVECTOR_ 1
//#define _USENRERRORCLASS_ 1
//#define _TURNONFPES_ 1

// all the system #include's we'll ever need
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

using namespace std;

// macro-like inline functions

template<class T>
inline T SQR(const T a) {return a*a;}

template<class T>
inline const T &MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

inline float MAX(const double &a, const float &b)
        {return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
        {return b > a ? float(b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}

inline float MIN(const double &a, const float &b)
        {return b < a ? (b) : float(a);}

inline float MIN(const float &a, const double &b)
        {return b < a ? float(b) : (a);}

template<class T>
inline T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
	{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}

// exception handling

#ifndef _USENRERRORCLASS_
#define throw(message) \
{printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); throw(1);}
#else
struct NRerror {
	char *message;
	char *file;
	int line;
	NRerror(char *m, char *f, int l) : message(m), file(f), line(l) {}
};
#define throw(message) throw(NRerror(message,__FILE__,__LINE__));
void NRcatch(NRerror err) {
	printf("ERROR: %s\n     in file %s at line %d\n",
		err.message, err.file, err.line);
	exit(1);
}
#endif

// usage example:
//
//	try {
//		somebadroutine();
//	}
//	catch(NRerror s) {NRcatch(s);}
//
// (You can of course substitute any other catch body for NRcatch(s).)


// Vector and Matrix Classes

#ifdef _USESTDVECTOR_
#define NRvector vector
#else

template <class T>
class NRvector {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	NRvector();
	explicit NRvector(int n);		// Zero-based array
	NRvector(int n, const T &a);	//initialize to constant value
	NRvector(int n, const T *a);	// Initialize to array
	NRvector(const NRvector &rhs);	// Copy constructor
	NRvector & operator=(const NRvector &rhs);	//assignment
	typedef T value_type; // make T available externally
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	void resize(int newn); // resize (contents not preserved)
	void assign(int newn, const T &a); // resize and assign a constant value
	~NRvector();
};

// NRvector definitions

template <class T>
NRvector<T>::NRvector() : nn(0), v(NULL) {}

template <class T>
NRvector<T>::NRvector(int n) : nn(n), v(n>0 ? new T[n] : NULL) {}

template <class T>
NRvector<T>::NRvector(int n, const T& a) : nn(n), v(n>0 ? new T[n] : NULL)
{
	for(int i=0; i<n; i++) v[i] = a;
}

template <class T>
NRvector<T>::NRvector(int n, const T *a) : nn(n), v(n>0 ? new T[n] : NULL)
{
	for(int i=0; i<n; i++) v[i] = *a++;
}

template <class T>
NRvector<T>::NRvector(const NRvector<T> &rhs) : nn(rhs.nn), v(nn>0 ? new T[nn] : NULL)
{
	for(int i=0; i<nn; i++) v[i] = rhs[i];
}

template <class T>
NRvector<T> & NRvector<T>::operator=(const NRvector<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			if (v != NULL) delete [] (v);
			nn=rhs.nn;
			v= nn>0 ? new T[nn] : NULL;
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

template <class T>
inline T & NRvector<T>::operator[](const int i)	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRvector subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline const T & NRvector<T>::operator[](const int i) const	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRvector subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline int NRvector<T>::size() const
{
	return nn;
}

template <class T>
void NRvector<T>::resize(int newn)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
}

template <class T>
void NRvector<T>::assign(int newn, const T& a)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
	for (int i=0;i<nn;i++) v[i] = a;
}

template <class T>
NRvector<T>::~NRvector()
{
	if (v != NULL) delete[] (v);
}

// end of NRvector definitions

#endif //ifdef _USESTDVECTOR_

template <class T>
class NRmatrix {
private:
	int nn;
	int mm;
	T **v;
public:
	NRmatrix();
	NRmatrix(int n, int m);			// Zero-based array
	NRmatrix(int n, int m, const T &a);	//Initialize to constant
	NRmatrix(int n, int m, const T *a);	// Initialize to array
	NRmatrix(const NRmatrix &rhs);		// Copy constructor
	NRmatrix & operator=(const NRmatrix &rhs);	//assignment
	typedef T value_type; // make T available externally
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline int nrows() const;
	inline int ncols() const;
	void resize(int newn, int newm); // resize (contents not preserved)
	void assign(int newn, int newm, const T &a); // resize and assign a constant value
	~NRmatrix();
};

template <class T>
NRmatrix<T>::NRmatrix() : nn(0), mm(0), v(NULL) {}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1;i<n;i++) v[i] = v[i-1] + m;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

template <class T>
NRmatrix<T>::NRmatrix(const NRmatrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL)
{
	int i,j,nel=mm*nn;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
NRmatrix<T> & NRmatrix<T>::operator=(const NRmatrix<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i,j,nel;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != NULL) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			v = nn>0 ? new T*[nn] : NULL;
			nel = mm*nn;
			if (v) v[0] = nel>0 ? new T[nel] : NULL;
			for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
		}
		for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
inline T* NRmatrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRmatrix subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline const T* NRmatrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRmatrix subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline int NRmatrix<T>::nrows() const
{
	return nn;
}

template <class T>
inline int NRmatrix<T>::ncols() const
{
	return mm;
}

template <class T>
void NRmatrix<T>::resize(int newn, int newm)
{
	int i,nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
}

template <class T>
void NRmatrix<T>::assign(int newn, int newm, const T& a)
{
	int i,j,nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::~NRmatrix()
{
	if (v != NULL) {
		delete[] (v[0]);
		delete[] (v);
	}
}

template <class T>
class NRMat3d {
private:
	int nn;
	int mm;
	int kk;
	T ***v;
public:
	NRMat3d();
	NRMat3d(int n, int m, int k);
	inline T** operator[](const int i);	//subscripting: pointer to row i
	inline const T* const * operator[](const int i) const;
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	~NRMat3d();
};

template <class T>
NRMat3d<T>::NRMat3d(): nn(0), mm(0), kk(0), v(NULL) {}

template <class T>
NRMat3d<T>::NRMat3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
	int i,j;
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
	for(i=1; i<n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
	}
}

template <class T>
inline T** NRMat3d<T>::operator[](const int i) //subscripting: pointer to row i
{
	return v[i];
}

template <class T>
inline const T* const * NRMat3d<T>::operator[](const int i) const
{
	return v[i];
}

template <class T>
inline int NRMat3d<T>::dim1() const
{
	return nn;
}

template <class T>
inline int NRMat3d<T>::dim2() const
{
	return mm;
}

template <class T>
inline int NRMat3d<T>::dim3() const
{
	return kk;
}

template <class T>
NRMat3d<T>::~NRMat3d()
{
	if (v != NULL) {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}


// basic type names (redefine if your bit lengths don't match)

typedef int Int; // 32 bit integer
typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

typedef char Char; // 8 bit integer
typedef unsigned char Uchar;

typedef double Doub; // default floating type
typedef long double Ldoub;

typedef complex<double> Complex; // default complex type

typedef bool Bool;

// NaN: uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

static const Doub NaN = numeric_limits<Doub>::quiet_NaN();

//Uint proto_nan[2]={0xffffffff, 0x7fffffff};
//double NaN = *( double* )proto_nan;

//Doub NaN = sqrt(-1.);

// vector types

typedef const NRvector<Int> VecInt_I;
typedef NRvector<Int> VecInt, VecInt_O, VecInt_IO;

typedef const NRvector<Uint> VecUint_I;
typedef NRvector<Uint> VecUint, VecUint_O, VecUint_IO;

typedef const NRvector<Llong> VecLlong_I;
typedef NRvector<Llong> VecLlong, VecLlong_O, VecLlong_IO;

typedef const NRvector<Ullong> VecUllong_I;
typedef NRvector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;

typedef const NRvector<Char> VecChar_I;
typedef NRvector<Char> VecChar, VecChar_O, VecChar_IO;

typedef const NRvector<Char*> VecCharp_I;
typedef NRvector<Char*> VecCharp, VecCharp_O, VecCharp_IO;

typedef const NRvector<Uchar> VecUchar_I;
typedef NRvector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;

typedef const NRvector<Doub> VecDoub_I;
typedef NRvector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef const NRvector<Doub*> VecDoubp_I;
typedef NRvector<Doub*> VecDoubp, VecDoubp_O, VecDoubp_IO;

typedef const NRvector<Complex> VecComplex_I;
typedef NRvector<Complex> VecComplex, VecComplex_O, VecComplex_IO;

typedef const NRvector<Bool> VecBool_I;
typedef NRvector<Bool> VecBool, VecBool_O, VecBool_IO;

// matrix types

typedef const NRmatrix<Int> MatInt_I;
typedef NRmatrix<Int> MatInt, MatInt_O, MatInt_IO;

typedef const NRmatrix<Uint> MatUint_I;
typedef NRmatrix<Uint> MatUint, MatUint_O, MatUint_IO;

typedef const NRmatrix<Llong> MatLlong_I;
typedef NRmatrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;

typedef const NRmatrix<Ullong> MatUllong_I;
typedef NRmatrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;

typedef const NRmatrix<Char> MatChar_I;
typedef NRmatrix<Char> MatChar, MatChar_O, MatChar_IO;

typedef const NRmatrix<Uchar> MatUchar_I;
typedef NRmatrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;

typedef const NRmatrix<Doub> MatDoub_I;
typedef NRmatrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;

typedef const NRmatrix<Bool> MatBool_I;
typedef NRmatrix<Bool> MatBool, MatBool_O, MatBool_IO;

// 3D matrix types

typedef const NRMat3d<Doub> Mat3DDoub_I;
typedef NRMat3d<Doub> Mat3DDoub, Mat3DDoub_O, Mat3DDoub_IO;

// Floating Point Exceptions for Microsoft compilers

#ifdef _TURNONFPES_
#ifdef _MSC_VER
struct turn_on_floating_exceptions {
	turn_on_floating_exceptions() {
		int cw = _controlfp( 0, 0 );
		cw &=~(EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE );
		_controlfp( cw, MCW_EM );
	}
};
turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif /* _MSC_VER */
#endif /* _TURNONFPES */

#endif /* _NR3_H_ */

struct Ran {
	Ullong u,v,w;
	Ran(Ullong j) : v(4101842887655102017LL), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}
	inline Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};
struct Ranq1 {
	Ullong v;
	Ranq1(Ullong j) : v(4101842887655102017LL) {
		v ^= j;
		v = int64();
	}
	inline Ullong int64() {
		v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
		return v * 2685821657736338717LL;
	}
	inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};
struct Ranq2 {
	Ullong v,w;
	Ranq2(Ullong j) : v(4101842887655102017LL), w(1) {
		v ^= j;
		w = int64();
		v = int64();
	}
	inline Ullong int64() {
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		return v ^ w;
	}
	inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};
struct Ranhash {
	inline Ullong int64(Ullong u) {
		Ullong v = u * 3935559000370003845LL + 2691343689449507681LL;
		v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
		v *= 4768777513237032717LL;
		v ^= v << 20; v ^= v >> 41; v ^= v << 5;
		return  v;
	}
	inline Uint int32(Ullong u)
		{ return (Uint)(int64(u) & 0xffffffff) ; }
	inline Doub doub(Ullong u)
		{ return 5.42101086242752217E-20 * int64(u); }
};
struct Ranbyte {
	Int s[256],i,j,ss;
	Uint v;
	Ranbyte(Int u) {
		v = 2244614371U ^ u;
		for (i=0; i<256; i++) {s[i] = i;}
		for (j=0, i=0; i<256; i++) {
			ss = s[i];
			j = (j + ss + (v >> 24)) & 0xff;
			s[i] = s[j]; s[j] = ss;
			v = (v << 24) | (v >> 8);
		}
		i = j = 0;
		for (Int k=0; k<256; k++) int8();
	}
	inline unsigned char int8() {
		i = (i+1) & 0xff;
		ss = s[i];
		j = (j+ss) & 0xff;
		s[i] = s[j]; s[j] = ss;
		return (unsigned char)(s[(s[i]+s[j]) & 0xff]);
	}
	Uint int32() {
		v = 0;
		for (int k=0; k<4; k++) {
			i = (i+1) & 0xff;
			ss = s[i];
			j = (j+ss) & 0xff;
			s[i] = s[j]; s[j] = ss;
			v = (v << 8) | s[(s[i]+s[j]) & 0xff];
		}
		return v;
	}
	Doub doub() {
		return 2.32830643653869629E-10 * ( int32() +
			   2.32830643653869629E-10 * int32() );
	}
};
struct Ranfib {
	Doub dtab[55], dd;
	Int inext, inextp;
	Ranfib(Ullong j) : inext(0), inextp(31) {
		Ranq1 init(j);
		for (int k=0; k<55; k++) dtab[k] = init.doub();
	}
	Doub doub() {
		if (++inext == 55) inext = 0;
		if (++inextp == 55) inextp = 0;
		dd = dtab[inext] - dtab[inextp];
		if (dd < 0) dd += 1.0;
		return (dtab[inext] = dd);
	}
	inline unsigned long int32()
		{ return (unsigned long)(doub() * 4294967295.0);}
};
struct Ranlim32 {
	Uint u,v,w1,w2;
	Ranlim32(Uint j) : v(2244614371U), w1(521288629U), w2(362436069U) {
		u = j ^ v; int32();
		v = u; int32();
	}
	inline Uint int32() {
		u = u * 2891336453U + 1640531513U;
		v ^= v >> 13; v ^= v << 17; v ^= v >> 5;
		w1 = 33378 * (w1 & 0xffff) + (w1 >> 16);
		w2 = 57225 * (w2 & 0xffff) + (w2 >> 16);
		Uint x = u ^ (u << 9); x ^= x >> 17; x ^= x << 6;
		Uint y = w1 ^ (w1 << 17); y ^= y >> 15; y ^= y << 5;
		return (x + v) ^ (y + w2);
	}
	inline Doub doub() { return 2.32830643653869629E-10 * int32(); }
	inline Doub truedoub() {
		return 2.32830643653869629E-10 * ( int32() +
		2.32830643653869629E-10 * int32() );
	}
};

struct Cholesky{
	Int n;
	MatDoub el;
	Cholesky(MatDoub_I &a) : n(a.nrows()), el(a) {
		Int i,j,k;
		VecDoub tmp;
		Doub sum;
		if (el.ncols() != n) throw("need square matrix");
		for (i=0;i<n;i++) {
			for (j=i;j<n;j++) {
				for (sum=el[i][j],k=i-1;k>=0;k--) sum -= el[i][k]*el[j][k];
				if (i == j) {
					if (sum <= 0.0)
						throw("Cholesky failed");
					el[i][i]=sqrt(sum);
				} else el[j][i]=sum/el[i][i];
			}
		}
		for (i=0;i<n;i++) for (j=0;j<i;j++) el[j][i] = 0.;
	}
	void solve(VecDoub_I &b, VecDoub_O &x) {
		Int i,k;
		Doub sum;
		if (b.size() != n || x.size() != n) throw("bad lengths in Cholesky");
		for (i=0;i<n;i++) {
			for (sum=b[i],k=i-1;k>=0;k--) sum -= el[i][k]*x[k];
			x[i]=sum/el[i][i];
		}
		for (i=n-1;i>=0;i--) {
			for (sum=x[i],k=i+1;k<n;k++) sum -= el[k][i]*x[k];
			x[i]=sum/el[i][i];
		}		
	}
	void elmult(VecDoub_I &y, VecDoub_O &b) {
		Int i,j;
		if (b.size() != n || y.size() != n) throw("bad lengths");
		for (i=0;i<n;i++) {
			b[i] = 0.;
			for (j=0;j<=i;j++) b[i] += el[i][j]*y[j];
		}
	}
	void elsolve(VecDoub_I &b, VecDoub_O &y) {
		Int i,j;
		Doub sum;
		if (b.size() != n || y.size() != n) throw("bad lengths");
		for (i=0;i<n;i++) {
			for (sum=b[i],j=0; j<i; j++) sum -= el[i][j]*y[j];
			y[i] = sum/el[i][i];
		}
	}
	void inverse(MatDoub_O &ainv) {
		Int i,j,k;
		Doub sum;
		ainv.resize(n,n);
		for (i=0;i<n;i++) for (j=0;j<=i;j++){
			sum = (i==j? 1. : 0.);
			for (k=i-1;k>=j;k--) sum -= el[i][k]*ainv[j][k];
			ainv[j][i]= sum/el[i][i];
		}
		for (i=n-1;i>=0;i--) for (j=0;j<=i;j++){
			sum = (i<j? 0. : ainv[j][i]);
			for (k=i+1;k<n;k++) sum -= el[k][i]*ainv[j][k];
			ainv[i][j] = ainv[j][i] = sum/el[i][i];
		}				
	}
	Doub logdet() {
		Doub sum = 0.;
		for (Int i=0; i<n; i++) sum += log(el[i][i]);
		return 2.*sum;
	}
};
Doub gammln(const Doub xx) {
	Int j;
	Doub x,tmp,y,ser;
	static const Doub cof[14]={57.1562356658629235,-59.5979603554754912,
	14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
	.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
	-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
	.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
	if (xx <= 0) throw("bad arg in gammln");
	y=x=xx;
	tmp = x+5.24218750000000000;
	tmp = (x+0.5)*log(tmp)-tmp;
	ser = 0.999999999999997092;
	for (j=0;j<14;j++) ser += cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}
Doub factrl(const Int n) {
	static VecDoub a(171);
	static Bool init=true;
	if (init) {
		init = false;
		a[0] = 1.;
		for (Int i=1;i<171;i++) a[i] = i*a[i-1];
	}
	if (n < 0 || n > 170) throw("factrl out of range");
	return a[n];
}
Doub factln(const Int n) {
	static const Int NTOP=2000;
	static VecDoub a(NTOP);
	static Bool init=true;
	if (init) {
		init = false;
		for (Int i=0;i<NTOP;i++) a[i] = gammln(i+1.);
	}
	if (n < 0) throw("negative arg in factln");
	if (n < NTOP) return a[n];
	return gammln(n+1.);
}
Doub bico(const Int n, const Int k) {
	if (n<0 || k<0 || k>n) throw("bad args in bico");
	if (n<171) return floor(0.5+factrl(n)/(factrl(k)*factrl(n-k)));
	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}
Doub beta(const Doub z, const Doub w) {
	return exp(gammln(z)+gammln(w)-gammln(z+w));
}
struct Gauleg18 {
	static const Int ngau = 18;
	static const Doub y[18];
	static const Doub w[18];
};
const Doub Gauleg18::y[18] = {0.0021695375159141994,
0.011413521097787704,0.027972308950302116,0.051727015600492421,
0.082502225484340941, 0.12007019910960293,0.16415283300752470,
0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
0.87126389619061517, 0.95698180152629142};
const Doub Gauleg18::w[18] = {0.0055657196642445571,
0.012915947284065419,0.020181515297735382,0.027298621498568734,
0.034213810770299537,0.040875750923643261,0.047235083490265582,
0.053244713977759692,0.058860144245324798,0.064039797355015485,
0.068745323835736408,0.072941885005653087,0.076598410645870640,
0.079687828912071670,0.082187266704339706,0.084078218979661945,
0.085346685739338721,0.085983275670394821};
struct Gamma : Gauleg18 {
	static const Int ASWITCH=100;
	static const Doub EPS;
	static const Doub FPMIN;
	Doub gln;

	Doub gammp(const Doub a, const Doub x) {
		if (x < 0.0 || a <= 0.0) throw("bad args in gammp");
		if (x == 0.0) return 0.0;
		else if ((Int)a >= ASWITCH) return gammpapprox(a,x,1);
		else if (x < a+1.0) return gser(a,x);
		else return 1.0-gcf(a,x);
	}

	Doub gammq(const Doub a, const Doub x) {
		if (x < 0.0 || a <= 0.0) throw("bad args in gammq");
		if (x == 0.0) return 1.0;
		else if ((Int)a >= ASWITCH) return gammpapprox(a,x,0);
		else if (x < a+1.0) return 1.0-gser(a,x);
		else return gcf(a,x);
	}

	Doub gser(const Doub a, const Doub x) {
		Doub sum,del,ap;
		gln=gammln(a);
		ap=a;
		del=sum=1.0/a;
		for (;;) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				return sum*exp(-x+a*log(x)-gln);
			}
		}
	}

	Doub gcf(const Doub a, const Doub x) {
		Int i;
		Doub an,b,c,d,del,h;
		gln=gammln(a);
		b=x+1.0-a;
		c=1.0/FPMIN;
		d=1.0/b;
		h=d;
		for (i=1;;i++) {
			an = -i*(i-a);
			b += 2.0;
			d=an*d+b;
			if (fabs(d) < FPMIN) d=FPMIN;
			c=b+an/c;
			if (fabs(c) < FPMIN) c=FPMIN;
			d=1.0/d;
			del=d*c;
			h *= del;
			if (fabs(del-1.0) <= EPS) break;
		}
		return exp(-x+a*log(x)-gln)*h;
	}

	Doub gammpapprox(Doub a, Doub x, Int psig) {
		Int j;
		Doub xu,t,sum,ans;
		Doub a1 = a-1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
		gln = gammln(a);
		if (x > a1) xu = MAX(a1 + 11.5*sqrta1, x + 6.0*sqrta1);
		else xu = MAX(0.,MIN(a1 - 7.5*sqrta1, x - 5.0*sqrta1));
		sum = 0;
		for (j=0;j<ngau;j++) {
			t = x + (xu-x)*y[j];
			sum += w[j]*exp(-(t-a1)+a1*(log(t)-lna1));
		}
		ans = sum*(xu-x)*exp(a1*(lna1-1.)-gln);
		return (psig?(ans>0.0? 1.0-ans:-ans):(ans>=0.0? ans:1.0+ans));
	}

	Doub invgammp(Doub p, Doub a);

};
const Doub Gamma::EPS = numeric_limits<Doub>::epsilon();
const Doub Gamma::FPMIN = numeric_limits<Doub>::min()/EPS;
Doub Gamma::invgammp(Doub p, Doub a) {
	Int j;
	Doub x,err,t,u,pp,lna1,afac,a1=a-1;
	const Doub EPS=1.e-8;
	gln=gammln(a);
	if (a <= 0.) throw("a must be pos in invgammap");
	if (p >= 1.) return MAX(100.,a + 100.*sqrt(a));
	if (p <= 0.) return 0.0;
	if (a > 1.) {
		lna1=log(a1);
		afac = exp(a1*(lna1-1.)-gln);
		pp = (p < 0.5)? p : 1. - p;
		t = sqrt(-2.*log(pp));
		x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
		if (p < 0.5) x = -x;
		x = MAX(1.e-3,a*pow(1.-1./(9.*a)-x/(3.*sqrt(a)),3));
	} else {
		t = 1.0 - a*(0.253+a*0.12);
		if (p < t) x = pow(p/t,1./a);
		else x = 1.-log(1.-(p-t)/(1.-t));
	}
	for (j=0;j<12;j++) {
		if (x <= 0.0) return 0.0;
		err = gammp(a,x) - p;
		if (a > 1.) t = afac*exp(-(x-a1)+a1*(log(x)-lna1));
		else t = exp(-x+a1*log(x)-gln);
		u = err/t;
		x -= (t = u/(1.-0.5*MIN(1.,u*((a-1.)/x - 1))));
		if (x <= 0.) x = 0.5*(x + t);
		if (fabs(t) < EPS*x ) break;
	}
	return x;
}
struct Beta : Gauleg18 {
	static const Int SWITCH=3000;
	static const Doub EPS, FPMIN;

	Doub betai(const Doub a, const Doub b, const Doub x) {
		Doub bt;
		if (a <= 0.0 || b <= 0.0) throw("Bad a or b in routine betai");
		if (x < 0.0 || x > 1.0) throw("Bad x in routine betai");
		if (x == 0.0 || x == 1.0) return x;
		if (a > SWITCH && b > SWITCH) return betaiapprox(a,b,x);
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
		if (x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
		else return 1.0-bt*betacf(b,a,1.0-x)/b;
	}

	Doub betacf(const Doub a, const Doub b, const Doub x) {
		Int m,m2;
		Doub aa,c,d,del,h,qab,qam,qap;
		qab=a+b;
		qap=a+1.0;
		qam=a-1.0;
		c=1.0;
		d=1.0-qab*x/qap;
		if (fabs(d) < FPMIN) d=FPMIN;
		d=1.0/d;
		h=d;
		for (m=1;m<10000;m++) {
			m2=2*m;
			aa=m*(b-m)*x/((qam+m2)*(a+m2));
			d=1.0+aa*d;
			if (fabs(d) < FPMIN) d=FPMIN;
			c=1.0+aa/c;
			if (fabs(c) < FPMIN) c=FPMIN;
			d=1.0/d;
			h *= d*c;
			aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
			d=1.0+aa*d;
			if (fabs(d) < FPMIN) d=FPMIN;
			c=1.0+aa/c;
			if (fabs(c) < FPMIN) c=FPMIN;
			d=1.0/d;
			del=d*c;
			h *= del;
			if (fabs(del-1.0) <= EPS) break;
		}
		return h;
	}

	Doub betaiapprox(Doub a, Doub b, Doub x) {
		Int j;
		Doub xu,t,sum,ans;
		Doub a1 = a-1.0, b1 = b-1.0, mu = a/(a+b);
		Doub lnmu=log(mu),lnmuc=log(1.-mu);
		t = sqrt(a*b/(SQR(a+b)*(a+b+1.0)));
		if (x > a/(a+b)) {
			if (x >= 1.0) return 1.0;
			xu = MIN(1.,MAX(mu + 10.*t, x + 5.0*t));
		} else {
			if (x <= 0.0) return 0.0;
			xu = MAX(0.,MIN(mu - 10.*t, x - 5.0*t));
		}
		sum = 0;
		for (j=0;j<18;j++) {
			t = x + (xu-x)*y[j];
			sum += w[j]*exp(a1*(log(t)-lnmu)+b1*(log(1-t)-lnmuc));
		}
		ans = sum*(xu-x)*exp(a1*lnmu-gammln(a)+b1*lnmuc-gammln(b)+gammln(a+b));
		return ans>0.0? 1.0-ans : -ans;
	}

	Doub invbetai(Doub p, Doub a, Doub b) {
		const Doub EPS = 1.e-8;
		Doub pp,t,u,err,x,al,h,w,afac,a1=a-1.,b1=b-1.;
		Int j;
		if (p <= 0.) return 0.;
		else if (p >= 1.) return 1.;
		else if (a >= 1. && b >= 1.) {
			pp = (p < 0.5)? p : 1. - p;
			t = sqrt(-2.*log(pp));
			x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
			if (p < 0.5) x = -x;
			al = (SQR(x)-3.)/6.;
			h = 2./(1./(2.*a-1.)+1./(2.*b-1.));
			w = (x*sqrt(al+h)/h)-(1./(2.*b-1)-1./(2.*a-1.))*(al+5./6.-2./(3.*h));
			x = a/(a+b*exp(2.*w));
		} else {
			Doub lna = log(a/(a+b)), lnb = log(b/(a+b));
			t = exp(a*lna)/a;
			u = exp(b*lnb)/b;
			w = t + u;
			if (p < t/w) x = pow(a*w*p,1./a);
			else x = 1. - pow(b*w*(1.-p),1./b);
		}
		afac = -gammln(a)-gammln(b)+gammln(a+b);
		for (j=0;j<10;j++) {
			if (x == 0. || x == 1.) return x;
			err = betai(a,b,x) - p;
			t = exp(a1*log(x)+b1*log(1.-x) + afac);
			u = err/t;
			x -= (t = u/(1.-0.5*MIN(1.,u*(a1/x - b1/(1.-x)))));
			if (x <= 0.) x = 0.5*(x + t);
			if (x >= 1.) x = 0.5*(x + t + 1.);
			if (fabs(t) < EPS*x && j > 0) break;
		}
		return x;
	}

};
const Doub Beta::EPS = numeric_limits<Doub>::epsilon();
const Doub Beta::FPMIN = numeric_limits<Doub>::min()/EPS;
struct Gammadist : Gamma {
	Doub alph, bet, fac;
	Gammadist(Doub aalph, Doub bbet = 1.) : alph(aalph), bet(bbet) {
		if (alph <= 0. || bet <= 0.) throw("bad alph,bet in Gammadist");
		fac = alph*log(bet)-gammln(alph);
	}
	Doub p(Doub x) {
		if (x <= 0.) throw("bad x in Gammadist");
		return exp(-bet*x+(alph-1.)*log(x)+fac);
	}
	Doub cdf(Doub x) {
		if (x < 0.) throw("bad x in Gammadist");
		return gammp(alph,bet*x);
	}
	Doub invcdf(Doub p) {
		if (p < 0. || p >= 1.) throw("bad p in Gammadist");
		return invgammp(p,alph)/bet;
	}
};
struct Betadist : Beta {
	Doub alph, bet, fac;
	Betadist(Doub aalph, Doub bbet) : alph(aalph), bet(bbet) {
		if (alph <= 0. || bet <= 0.) throw("bad alph,bet in Betadist");
		fac = gammln(alph+bet)-gammln(alph)-gammln(bet);
	}
	Doub p(Doub x) {
		if (x <= 0. || x >= 1.) throw("bad x in Betadist");
		return exp((alph-1.)*log(x)+(bet-1.)*log(1.-x)+fac);
	}
	Doub cdf(Doub x) {
		if (x < 0. || x > 1.) throw("bad x in Betadist");
		return betai(alph,bet,x);
	}
	Doub invcdf(Doub p) {
		if (p < 0. || p > 1.) throw("bad p in Betadist");
		return invbetai(p,alph,bet);
	}
};
struct Studenttdist : Beta {
	Doub nu, mu, sig, np, fac;
	Studenttdist(Doub nnu, Doub mmu = 0., Doub ssig = 1.)
	: nu(nnu), mu(mmu), sig(ssig) {
		if (sig <= 0. || nu <= 0.) throw("bad sig,nu in Studentdist");
		np = 0.5*(nu + 1.);
		fac = gammln(np)-gammln(0.5*nu);
	}
	Doub p(Doub t) {
		return exp(-np*log(1.+SQR((t-mu)/sig)/nu)+fac)
			/(sqrt(3.14159265358979324*nu)*sig);
	}
	Doub cdf(Doub t) {
		Doub p = 0.5*betai(0.5*nu, 0.5, nu/(nu+SQR((t-mu)/sig)));
		if (t >= mu) return 1. - p;
		else return p;
	}
	Doub invcdf(Doub p) {
		if (p <= 0. || p >= 1.) throw("bad p in Studentdist");
		Doub x = invbetai(2.*MIN(p,1.-p), 0.5*nu, 0.5);
		x = sig*sqrt(nu*(1.-x)/x);
		return (p >= 0.5? mu+x : mu-x);
	}
	Doub aa(Doub t) {
		if (t < 0.) throw("bad t in Studentdist");
		return 1.-betai(0.5*nu, 0.5, nu/(nu+SQR(t)));
	}
	Doub invaa(Doub p) {
		if (p < 0. || p >= 1.) throw("bad p in Studentdist");
		Doub x = invbetai(1.-p, 0.5*nu, 0.5);
		return sqrt(nu*(1.-x)/x);
	}
};
struct Poissondist : Gamma {
	Doub lam;
	Poissondist(Doub llam) : lam(llam) {
		if (lam <= 0.) throw("bad lam in Poissondist");	
	}
	Doub p(Int n) {
		if (n < 0) throw("bad n in Poissondist");
		return exp(-lam + n*log(lam) - gammln(n+1.));
	}
	Doub cdf(Int n) {
		if (n < 0) throw("bad n in Poissondist");
		if (n == 0) return 0.;
		return gammq((Doub)n,lam);
	}
	Int invcdf(Doub p) {
		Int n,nl,nu,inc=1;
		if (p <= 0. || p >= 1.) throw("bad p in Poissondist");
		if (p < exp(-lam)) return 0;
		n = (Int)MAX(sqrt(lam),5.);
		if (p < cdf(n)) {
			do {
				n = MAX(n-inc,0);
				inc *= 2;
			} while (p < cdf(n));
			nl = n; nu = n + inc/2;
		} else {
			do {
				n += inc;
				inc *= 2;
			} while (p > cdf(n));
			nu = n; nl = n - inc/2;
		}
		while (nu-nl>1) {
			n = (nl+nu)/2;
			if (p < cdf(n)) nu = n;
			else nl = n;
		}
		return nl;
	}
};
struct Binomialdist : Beta {
	Int n;
	Doub pe, fac;
	Binomialdist(Int nn, Doub ppe) : n(nn), pe(ppe) {
		if (n <= 0 || pe <= 0. || pe >= 1.) throw("bad args in Binomialdist");
		fac = gammln(n+1.);
	}
	Doub p(Int k) {
		if (k < 0) throw("bad k in Binomialdist");
		if (k > n) return 0.;
		return exp(k*log(pe)+(n-k)*log(1.-pe)
			+fac-gammln(k+1.)-gammln(n-k+1.));
	}
	Doub cdf(Int k) {
		if (k < 0) throw("bad k in Binomialdist");
		if (k == 0) return 0.;
		if (k > n) return 1.;
		return 1. - betai((Doub)k,n-k+1.,pe);
	}
	Int invcdf(Doub p) {
		Int k,kl,ku,inc=1;
		if (p <= 0. || p >= 1.) throw("bad p in Binomialdist");
		k = MAX(0,MIN(n,(Int)(n*pe)));
		if (p < cdf(k)) {
			do {
				k = MAX(k-inc,0);
				inc *= 2;
			} while (p < cdf(k));
			kl = k; ku = k + inc/2;
		} else {
			do {
				k = MIN(k+inc,n+1);
				inc *= 2;
			} while (p > cdf(k));
			ku = k; kl = k - inc/2;
		}
		while (ku-kl>1) {
			k = (kl+ku)/2;
			if (p < cdf(k)) ku = k;
			else kl = k;
		}
		return kl;
	}
};
struct Chisqdist : Gamma {
	Doub nu,fac;
	Chisqdist(Doub nnu) : nu(nnu) {
		if (nu <= 0.) throw("bad nu in Chisqdist");
		fac = 0.693147180559945309*(0.5*nu)+gammln(0.5*nu);
	}
	Doub p(Doub x2) {
		if (x2 <= 0.) throw("bad x2 in Chisqdist");
		return exp(-0.5*(x2-(nu-2.)*log(x2))-fac);
	}
	Doub cdf(Doub x2) {
		if (x2 < 0.) throw("bad x2 in Chisqdist");
		return gammp(0.5*nu,0.5*x2);
	}
	Doub invcdf(Doub p) {
		if (p < 0. || p >= 1.) throw("bad p in Chisqdist");
		return 2.*invgammp(p,0.5*nu);
	}
};
struct Fdist : Beta {
	Doub nu1,nu2;
	Doub fac;
	Fdist(Doub nnu1, Doub nnu2) : nu1(nnu1), nu2(nnu2) {
		if (nu1 <= 0. || nu2 <= 0.) throw("bad nu1,nu2 in Fdist");
		fac = 0.5*(nu1*log(nu1)+nu2*log(nu2))+gammln(0.5*(nu1+nu2))
			-gammln(0.5*nu1)-gammln(0.5*nu2);
	}
	Doub p(Doub f) {
		if (f <= 0.) throw("bad f in Fdist");
		return exp((0.5*nu1-1.)*log(f)-0.5*(nu1+nu2)*log(nu2+nu1*f)+fac);
	}
	Doub cdf(Doub f) {
		if (f < 0.) throw("bad f in Fdist");
		return betai(0.5*nu1,0.5*nu2,nu1*f/(nu2+nu1*f));
	}
	Doub invcdf(Doub p) {
		if (p <= 0. || p >= 1.) throw("bad p in Fdist");
		Doub x = invbetai(p,0.5*nu1,0.5*nu2);
		return nu2*x/(nu1*(1.-x));
	}
};
struct preGaumixmod {
	static Int mmstat;
	struct Mat_mm : MatDoub {Mat_mm() : MatDoub(mmstat,mmstat) {} };
	preGaumixmod(Int mm) {mmstat = mm;}
};
Int preGaumixmod::mmstat = -1;

struct Gaumixmod : preGaumixmod {
	Int nn, kk, mm;
	MatDoub data, means, resp;
	VecDoub frac, lndets;
	vector<Mat_mm> sig;
	Doub loglike;
	Gaumixmod(MatDoub &ddata, MatDoub &mmeans) : preGaumixmod(ddata.ncols()),
	nn(ddata.nrows()), kk(mmeans.nrows()), mm(mmstat), data(ddata), means(mmeans),
	resp(nn,kk), frac(kk), lndets(kk), sig(kk) {
		Int i,j,k;
		for (k=0;k<kk;k++) {
			frac[k] = 1./kk;
			for (i=0;i<mm;i++) {
				for (j=0;j<mm;j++) sig[k][i][j] = 0.;
				sig[k][i][i] = 1.0e-5; /*1.0e-10;*/
			}
		}
		estep();
		mstep();
	}
	Doub estep() {
		Int k,m,n;
		Doub tmp,sum,max,oldloglike;
		VecDoub u(mm),v(mm);
		oldloglike = loglike;
		for (k=0;k<kk;k++) {
			Cholesky choltmp(sig[k]);
			lndets[k] = choltmp.logdet();
			for (n=0;n<nn;n++) {
				for (m=0;m<mm;m++) u[m] = data[n][m]-means[k][m];
				choltmp.elsolve(u,v);
				for (sum=0.,m=0; m<mm; m++) sum += SQR(v[m]);
				resp[n][k] = -0.5*(sum + lndets[k]) + log(frac[k]);
			}
		}
		loglike = 0;
		for (n=0;n<nn;n++) {
			max = -99.9e99;
			for (k=0;k<kk;k++) if (resp[n][k] > max) max = resp[n][k];
			for (sum=0.,k=0; k<kk; k++) sum += exp(resp[n][k]-max);
			tmp = max + log(sum);
			for (k=0;k<kk;k++) resp[n][k] = exp(resp[n][k] - tmp);
			loglike +=tmp;
		}
		return loglike - oldloglike;
	}
	void mstep() {
		Int j,n,k,m;
		Doub wgt,sum;
		for (k=0;k<kk;k++) {
			wgt=0.;
			for (n=0;n<nn;n++) wgt += resp[n][k];
			frac[k] = wgt/nn;
			for (m=0;m<mm;m++) {
				for (sum=0.,n=0; n<nn; n++) sum += resp[n][k]*data[n][m];
				means[k][m] = sum/wgt;
				for (j=0;j<mm;j++) {
					for (sum=0.,n=0; n<nn; n++) {
						sum += resp[n][k]*
							(data[n][m]-means[k][m])*(data[n][j]-means[k][j]);
					}
					sig[k][m][j] = sum/wgt;
				}
			}
		}
	}
};
struct Expondev : Ran {
	Doub beta;
	Expondev(Doub bbeta, Ullong i) : Ran(i), beta(bbeta) {}
	Doub dev() {
		Doub u;
		do u = doub(); while (u == 0.);
		return -log(u)/beta;
	}
};
struct Logisticdev : Ran {
	Doub mu,sig;
	Logisticdev(Doub mmu, Doub ssig, Ullong i) : Ran(i), mu(mmu), sig(ssig) {}
	Doub dev() {
		Doub u;
		do u = doub(); while (u*(1.-u) == 0.);
		return mu + 0.551328895421792050*sig*log(u/(1.-u));
	}
};
struct Normaldev_BM : Ran {
	Doub mu,sig;
	Doub storedval;
	Normaldev_BM(Doub mmu, Doub ssig, Ullong i)
	: Ran(i), mu(mmu), sig(ssig), storedval(0.) {}
	Doub dev() {
		Doub v1,v2,rsq,fac;
		if (storedval == 0.) {
			do {
				v1=2.0*doub()-1.0;
				v2=2.0*doub()-1.0;
				rsq=v1*v1+v2*v2;
			} while (rsq >= 1.0 || rsq == 0.0);
			fac=sqrt(-2.0*log(rsq)/rsq);
			storedval = v1*fac;
			return mu + sig*v2*fac;
		} else {
			fac = storedval;
			storedval = 0.;
			return mu + sig*fac;
		}
	}
};
struct Cauchydev : Ran {
	Doub mu,sig;
	Cauchydev(Doub mmu, Doub ssig, Ullong i) : Ran(i), mu(mmu), sig(ssig) {}
	Doub dev() {
		Doub v1,v2;
		do {
			v1=2.0*doub()-1.0;
			v2=doub();
		} while (SQR(v1)+SQR(v2) >= 1. || v2 == 0.);
		return mu + sig*v1/v2;
	}
};
struct Normaldev : Ran {
	Doub mu,sig;
	Normaldev(Doub mmu, Doub ssig, Ullong i)
	: Ran(i), mu(mmu), sig(ssig){}
	Doub dev() {
		Doub u,v,x,y,q;
		do {
			u = doub();
			v = 1.7156*(doub()-0.5);
			x = u - 0.449871;
			y = abs(v) + 0.386595;
			q = SQR(x) + y*(0.19600*y-0.25472*x);
		} while (q > 0.27597
			&& (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));
		return mu + sig*v/u;
	}
};
struct Gammadev : Normaldev {
	Doub alph, oalph, bet;
	Doub a1,a2;
	Gammadev(Doub aalph, Doub bbet, Ullong i)
	: Normaldev(0.,1.,i), alph(aalph), oalph(aalph), bet(bbet) {
		if (alph <= 0.) throw("bad alph in Gammadev");
		if (alph < 1.) alph += 1.;
		a1 = alph-1./3.;
		a2 = 1./sqrt(9.*a1);
	}
	Doub dev() {
		Doub u,v,x;
		do {
			do {
				x = Normaldev::dev();
				v = 1. + a2*x;
			} while (v <= 0.);
			v = v*v*v;
			u = doub();
		} while (u > 1. - 0.331*SQR(SQR(x)) &&
			log(u) > 0.5*SQR(x) + a1*(1.-v+log(v)));
		if (alph == oalph) return a1*v/bet;
		else {
			do u=doub(); while (u == 0.);
			return pow(u,1./oalph)*a1*v/bet;
		}
	}
};
struct Poissondev : Ran {
	Doub lambda, sqlam, loglam, lamexp, lambold;
	VecDoub logfact;
	Int swch;
	Poissondev(Doub llambda, Ullong i) : Ran(i), lambda(llambda),
		logfact(1024,-1.), lambold(-1.) {}
	Int dev() {
		Doub u,u2,v,v2,p,t,lfac;
		Int k;
		if (lambda < 5.) {
			if (lambda != lambold) lamexp=exp(-lambda);
			k = -1;
			t=1.;
			do {
				++k;
				t *= doub();
			} while (t > lamexp);
		} else {
			if (lambda != lambold) {
				sqlam = sqrt(lambda);
				loglam = log(lambda);
			}
			for (;;) {
				u = 0.64*doub();
				v = -0.68 + 1.28*doub();
				if (lambda > 13.5) {
					v2 = SQR(v);
					if (v >= 0.) {if (v2 > 6.5*u*(0.64-u)*(u+0.2)) continue;}
					else {if (v2 > 9.6*u*(0.66-u)*(u+0.07)) continue;}
				}
				k = Int(floor(sqlam*(v/u)+lambda+0.5));
				if (k < 0) continue;
				u2 = SQR(u);
				if (lambda > 13.5) {
					if (v >= 0.) {if (v2 < 15.2*u2*(0.61-u)*(0.8-u)) break;}
					else {if (v2 < 6.76*u2*(0.62-u)*(1.4-u)) break;}
				}
				if (k < 1024) {
					if (logfact[k] < 0.) logfact[k] = gammln(k+1.);
					lfac = logfact[k];
				} else lfac = gammln(k+1.);
				p = sqlam*exp(-lambda + k*loglam - lfac);
				if (u2 < p) break;
			}
		}
		lambold = lambda;
		return k;
	}
	Int dev(Doub llambda) {
		lambda = llambda;
		return dev();
	}
};
struct Binomialdev : Ran {
	Doub pp,p,pb,expnp,np,glnp,plog,pclog,sq;
	Int n,swch;
	Ullong uz,uo,unfin,diff,rltp;
	Int pbits[5];
	Doub cdf[64];
	Doub logfact[1024];
	Binomialdev(Int nn, Doub ppp, Ullong i) : Ran(i), pp(ppp), n(nn) {
		Int j;
		pb = p = (pp <= 0.5 ? pp : 1.0-pp);
		if (n <= 64) {
			uz=0;
			uo=0xffffffffffffffffLL;
			rltp = 0;
			for (j=0;j<5;j++) pbits[j] = 1 & ((Int)(pb *= 2.));
			pb -= floor(pb);
			swch = 0;
		} else if (n*p < 30.) {
			cdf[0] = exp(n*log(1-p));
			for (j=1;j<64;j++) cdf[j] =  cdf[j-1] + exp(gammln(n+1.)
				-gammln(j+1.)-gammln(n-j+1.)+j*log(p)+(n-j)*log(1.-p));
			swch = 1;
		} else {
			np = n*p;
			glnp=gammln(n+1.);
			plog=log(p);
			pclog=log(1.-p);
			sq=sqrt(np*(1.-p));
			if (n < 1024) for (j=0;j<=n;j++) logfact[j] = gammln(j+1.);
			swch = 2;
		}
	}	
	Int dev() {
		Int j,k,kl,km;
		Doub y,u,v,u2,v2,b;
		if (swch == 0) {
			unfin = uo;
			for (j=0;j<5;j++) {
				diff = unfin & (int64()^(pbits[j]? uo : uz));
				if (pbits[j]) rltp |= diff;
				else rltp = rltp & ~diff;
				unfin = unfin & ~diff;
			}
			k=0;
			for (j=0;j<n;j++) {
				if (unfin & 1) {if (doub() < pb) ++k;}
				else {if (rltp & 1) ++k;}
				unfin >>= 1;
				rltp >>= 1;
			}
		} else if (swch == 1) {
			y = doub();
			kl = -1;
			k = 64;
			while (k-kl>1) {
				km = (kl+k)/2;
				if (y < cdf[km]) k = km;
				else kl = km;
			}
		} else {
			for (;;) {
				u = 0.645*doub();
				v = -0.63 + 1.25*doub();
				v2 = SQR(v);
				if (v >= 0.) {if (v2 > 6.5*u*(0.645-u)*(u+0.2)) continue;}
				else {if (v2 > 8.4*u*(0.645-u)*(u+0.1)) continue;}
				k = Int(floor(sq*(v/u)+np+0.5));
				if (k < 0) continue;
				u2 = SQR(u);
				if (v >= 0.) {if (v2 < 12.25*u2*(0.615-u)*(0.92-u)) break;}
				else {if (v2 < 7.84*u2*(0.615-u)*(1.2-u)) break;}
				b = sq*exp(glnp+k*plog+(n-k)*pclog
					- (n < 1024 ? logfact[k]+logfact[n-k]
						: gammln(k+1.)+gammln(n-k+1.)));
				if (u2 < b) break;
			}
		}
		if (p != pp) k = n - k;
		return k;
	}
};
