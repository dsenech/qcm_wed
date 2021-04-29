#ifndef Q_matrix_h
#define Q_matrix_h

#include <iostream>
#include <iomanip>
#include <cstdio>

#include "parser.hpp"
#include "block_matrix.hpp"

//! Used to store the Lehman representation of a part of the Green function (a given symmetry block)
template<typename HilbertField>
struct Q_matrix
{

	size_t L; //!< number of orbital indices (number of rows)
	size_t M; //!< number of eigenvalues (number of columns)
	vector<double> e; //!< eigenvalues
	matrix<HilbertField> v; //!< matrix itself

  Q_matrix();
  Q_matrix(size_t _L, size_t _M);
  Q_matrix(const Q_matrix<HilbertField> &q);
  Q_matrix(const vector<double>& _e, const matrix<HilbertField>& _v);
  Q_matrix(istream& fin);
  Q_matrix<HilbertField>& operator=(const Q_matrix<HilbertField> &q);
  void append(Q_matrix &q);
  void check_norm(double threshold, double norm = 1.0);
  void Green_function(const Complex &z, matrix<Complex> &G);
  void integrated_Green_function(matrix<Complex> &G);
  void Spectral_function(const Complex &z, matrix<Complex> &G);
  void streamline();
};


//==============================================================================
// implementation


/**
 default constructor
 */
template<typename HilbertField>
Q_matrix<HilbertField>::Q_matrix(): L(0), M(0), v() {}

/**
 constructor from sizes
 */
template<typename HilbertField>
Q_matrix<HilbertField>::Q_matrix(size_t _L, size_t _M): L(_L), M(_M) {
  e.resize(M);
  v.set_size(_L,_M);
}

/**
 copy constructor
 */
template<typename HilbertField>
Q_matrix<HilbertField>::Q_matrix(const Q_matrix<HilbertField> &q)
{
  L = q.L;
  M = q.M;
  e = q.e;
  v = q.v;
}

/**
 constructor from arrays
 */
template<typename HilbertField>
Q_matrix<HilbertField>::Q_matrix(const vector<double>& _e, const matrix<HilbertField>& _v) : e(_e), v(_v)
{
  assert(e.size() == v.r);
  L = v.c;
  M = e.size();
}

/**
 Constructor from input stream (ASCII file)
 */
template<typename HilbertField>
Q_matrix<HilbertField>::Q_matrix(std::istream &flux){
  size_t i,j;
  
  parser::find_next(flux, "w");
  flux >> L >> M;
  
  assert(M > 0 and L > 0);
  v.set_size(L,M);
  e.resize(M);
  
  for(i=0; i<M; ++i){
    flux >> e[i];
    for(j=0; j<L; ++j) flux >> v(j,i);
  }
}



template<typename HilbertField>
Q_matrix<HilbertField>& Q_matrix<HilbertField>::operator=(const Q_matrix<HilbertField> &q){
  L = q.L;
  M = q.M;
  e = q.e;
  v = q.v;
  return *this;
}

/**
 Appends to the qmatrix another one (increasing the number of columns for the same number of rows)
 @param q q_matrix to append
 */
template<typename HilbertField>
void Q_matrix<HilbertField>::append(Q_matrix &q)
{
  if(q.M==0) return;
  M += q.M;
  if(e.size()>0){
    assert(L == q.L);
    vector<double> e_temp = e;
    matrix<HilbertField> v_temp(v);
    v.concatenate(v_temp, q.v);
    e.resize(e_temp.size()+q.e.size());
    copy(&e_temp[0],&e_temp[0]+e_temp.size(),&e[0]);
    copy(&q.e[0],&q.e[0]+q.e.size(),&e[e_temp.size()]);
  }
  else{
    L = q.L;
    v = q.v;
    e = q.e;
  }
}


/**
 partial Green function evaluation
 @param z complex frequency
 @param [out] G Green function (adds to previous value)
 */
template<typename HilbertField>
void Q_matrix<HilbertField>::Green_function(const Complex &z, matrix<Complex> &G)
{  
  for(size_t i=0; i<M; ++i){
    Complex u = (1.0/(z-e[i]));
    for(size_t a=0; a<L; ++a){
      for(size_t b=0; b<L; ++b){
        G(b,a) += v(a,i)*conjugate(v(b,i))*u; // original was G(a,b) but was wrong with complex operators
      }
    }
  }
}




/**
 integrated Green function evaluation
 @param [out] G integrated Green function (adds to previous value)
 */
template<typename HilbertField>
void Q_matrix<HilbertField>::integrated_Green_function(matrix<Complex> &G)
{
  for(size_t i=0; i<M; ++i){
    if(e[i] >= 0.0) continue;
    for(size_t a=0; a<L; ++a){
      for(size_t b=0; b<L; ++b){
        G(b,a) += v(a,i)*conjugate(v(b,i)); // original was G(a,b) but was wrong with complex operators
      }
    }
  }
}




/**
 partial Spectral function evaluation
 clears G before filling
 @param z complex frequency
 @param [out] G Green function
 */
template<typename HilbertField>
void Q_matrix<HilbertField>::Spectral_function(const Complex &z, matrix<Complex> &G)
{
  G.zero();
  for(size_t i=0; i<M; ++i){
    double u = -(1.0/(z-e[i])).imag();
    for(size_t a=0; a<L; ++a){
      for(size_t b=0; b<L; ++b){
        G(a,b) += v(a,i)*conjugate(v(b,i))*u;
      }
    }
  }
}







/**
 Eliminates the small contributions from the Q matrix
 */
template<typename HilbertField>
void Q_matrix<HilbertField>::streamline()
{
  vector<int> f(M);
  double Qmatrix_tolerance = global_double("Qmatrix_tolerance");
  
  size_t k=0;
  for(size_t i=0; i<M; ++i){
    size_t j;
    if(f[i]) continue;
    for(j=0; j<L; ++j) if(abs(v(j,i)) > Qmatrix_tolerance) break;
    if(j==L) f[i]=1;
    else k++;
  }
  
  if(k==M) return;
  else{
    size_t M_old=M;
    M = k;
    vector<double> e_temp = e;
    matrix<HilbertField> v_temp = v;
    v.set_size(L,M);
    e.resize(M);
    for(size_t i=0, k=0; i<M_old; ++i){
      if(f[i]) continue;
      e[k] = e_temp[i];
      for(size_t j=0; j<L; ++j) v(j,k) = v_temp(j,i);
      k++;
    }
  }
}






/**
 Checks the normalization
 @param threshold limit under which a normalization violation is signaled
 @param norm expected norm of the Q matrix (1.0 by default)
 */
template<typename HilbertField>
void Q_matrix<HilbertField>::check_norm(double threshold, double norm)
{
  matrix<HilbertField> tmp(L);
  HilbertField z, ztot=0.0;
  for(size_t i=0; i<L; ++i){
    z = 0.0;
    for(size_t k=0; k<M; ++k) z += v(i,k)*conjugate(v(i,k));
    tmp(i,i) = chop(z-1.0);
    ztot += (z-norm)*conjugate(z-norm);
    for(size_t j=0; j<i; ++j){
      HilbertField zz = 0.0;
      for(int k=0; k<M; ++k) zz += v(i,k)*conjugate(v(j,k));
      tmp(i,j) = chop(zz);
      tmp(j,i) = chop(conjugate(zz));
      ztot += 2.0*zz*conjugate(zz); // correction (remarquée par A. Foley)
    }
  }
  double ztot2 = sqrt(realpart(ztot));
  if(ztot2 > threshold){
    qcm_ED_throw("Q-matrix does not satisfy unitary condition by "+to_string(ztot2));
  }
}




/**
 Printing to ASCII file
 */
template<typename HilbertField>
std::ostream & operator<<(std::ostream &flux, const Q_matrix<HilbertField> &Q){
  size_t i,j;
  
  flux << "w" << std::setprecision(LONG_DISPLAY);
  flux << '\t' << Q.L << '\t' << Q.M << std::endl;
  for(i=0; i<Q.M; ++i){
    flux << Q.e[i];
    for(j=0; j<Q.L; ++j) flux << '\t' << std::setprecision(LONG_DISPLAY) << chop(Q.v(j,i));
    flux << std::endl;
  }
  return flux;
}




#endif
