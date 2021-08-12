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
  void streamline();
  matrix<HilbertField> streamlineQ(const matrix<HilbertField> &Q, double threshold);
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
        G(a,b) += v(a,i)*conjugate(v(b,i))*u; 
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
        G(a,b) += v(a,i)*conjugate(v(b,i));
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
  
  vector<pair<double, matrix<HilbertField>>> LR;
  LR.reserve(M);

  matrix<HilbertField> v2(L,M);
  vector<double> e2(M);
  int ii=0;

  for(int i=0; i<M; i++){
    int j;
    for(j=i+1; j<M; j++){
      if(fabs(e[j]-e[i]) > Qmatrix_tolerance) break;
    }
    int n = j-i;
    matrix<HilbertField> Q(L, n);
    v.move_sub_matrix(L, n, 0, i, 0, 0, Q);
    matrix<HilbertField> Us = streamlineQ(Q, 1e-8);
    if(Us.c > 0){
      LR.push_back({e[i], Us});
      // cout << "e = " << e[i] << ", Us = " << Us << endl;
      memcpy(&e2[ii], &e[i], Us.c*sizeof(double));
      memcpy(&v2.v[L*ii], &Us.v[0], Us.c*L*sizeof(HilbertField));
      ii += Us.c;
    }
    i += n-1;
  }
  M = ii;
  e = e2;
  e.resize(M);
  v.v = v2.v;
  v.v.resize(M*L);
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
    tmp(i,i) = z-norm;
    ztot += (z-norm)*conjugate(z-norm);
    for(size_t j=0; j<i; ++j){
      HilbertField zz = 0.0;
      for(int k=0; k<M; ++k) zz += v(i,k)*conjugate(v(j,k));
      tmp(i,j) = zz;
      tmp(j,i) = conjugate(zz);
      ztot += 2.0*zz*conjugate(zz); // correction (remarquée par A. Foley)
    }
  }
  double ztot2 = sqrt(realpart(ztot));
  if(console::level > 5) cout << "Q-matrix norm error : " << ztot2 << endl;
  if(ztot2 > threshold){
    cout << "faulty Q-matrix:\n" << *this << endl;
    cout << "Q-matrix norm is off by :\n" << tmp << endl;
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

/**
From Q^dagger Q, eliminates the small contributions and returns a new Q with fewer columns
*/
template<typename HilbertField>
matrix<HilbertField> Q_matrix<HilbertField>::streamlineQ(const matrix<HilbertField> &Q, double threshold)
{
  // computing QQ = Q*Q^\dagger
  matrix<HilbertField> QQ(L);
  for(int i=0; i<L; i++){
    QQ(i,i) = 0;
    for(int k=0; k<Q.c; k++){
      QQ(i,i) += Q(i,k)*conjugate(Q(i,k));
    }
    for(int j=0; j<i; j++){
      QQ(i,j) = 0;
      for(int k=0; k<Q.c; k++){
        QQ(i,j) += Q(i,k)*conjugate(Q(j,k));
      }
      QQ(j,i) = conjugate(QQ(i,j));
    }
  }
  // eigendecomposition
  vector<double> d(L);
  matrix<HilbertField> U(L);
  QQ.eigensystem(d, U);
  // cout << "streamlining eigenvalues from " << Q.c << " terms : " << d << endl;
  vector<bool> I(0);
  I.reserve(L);
  for(int i=0; i<L; i++) if(fabs(d[i]) > threshold) I.push_back(i);
  matrix<HilbertField> Us(L,I.size());
  for(int i=0; i<I.size(); i++){
    //! extracts a rectangular submatrix of size (R,C), starting at (i,j) of this, and copies it to position (I, J) of A
      U.move_sub_matrix(L, 1, 0, I[i], 0, i, Us, sqrt(d[I[i]]));
  }
  return Us;
}

#endif
