#ifndef CSR_hermitian_h
#define CSR_hermitian_h

#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

//! sparse matrix in compressed-sparse-row (CSR) format
/**
 template (type is expected to be double or complex<double>)
 if sym_store is false, only the lower triangle and the diagonal are stored
 otherwise the full matrix is stored, which takes almost twice as much memory
 but allows for faster application on a vector via openMP
 */
template <typename T>
struct CSR_hermitian {
	bool sym_store;  //!< True if the matrix is stored symmetrically, instead of only the upper diagonal
	vector<int32_t> J; //!< array of column indices
	vector<int32_t> Iptr; //!< array of indices of J where the row index changes
	vector<T> v; //!< array of values of the matrix elements
	vector<double> diag; //!< diagonal elements (stored separately)
	
  
  /**
   Applies the matrix on a vector x and adds the result to vector y
   */
	void apply(const vector<T> &x, vector<T> &y){
		// first multiply by the diagonal elements
		for(uint32_t i=0; i<diag.size(); i++) y[i] += diag[i]*x[i];
		
		uint32_t imax = Iptr.size()-1;
    // parallelize the action on a vector via openMP (sym_store = true)
		if(sym_store){
#ifdef _OPENMP
			int level = omp_get_level();
			int num=1;
			for(int i=0; i<level; i++) num *= omp_get_team_size(i+1);
			num = omp_get_max_threads()/num;
#endif
#pragma omp parallel for schedule(guided) num_threads(num)
			for(uint32_t i=0; i<imax; i++){
				uint32_t jmax = Iptr[i+1];
				for(uint32_t j=Iptr[i]; j<jmax; j++) y[i] += x[J[j]]*v[j];
			}
		}
		
		else{
			// then multiply by the upper triangle
			for(uint32_t i=0; i<imax; i++){
				uint32_t jmax = Iptr[i+1];
				for(uint32_t j=Iptr[i]; j<jmax; j++) y[i] += x[J[j]]*v[j];
			}
			
			// then by the lower triangle
			for(uint32_t i=0; i<imax; i++){
				uint32_t jmax = Iptr[i+1];
				for(uint32_t j=Iptr[i]; j<jmax; j++) y[J[j]] += x[i]*conjugate(v[j]);
			}
		}
	}
	
  /**
   Applies the matrix on a vector x and adds z times the result to vector y
   */
	void apply(const vector<T> &x, vector<T> &y, T z){
		// first multiply by the diagonal elements
		for(uint32_t i=0; i<diag.size(); i++) y[i] += z*diag[i]*x[i];
		
		uint32_t imax = Iptr.size()-1;
		if(sym_store){
#ifdef _OPENMP
			int level = omp_get_level();
			int num=1;
			for(int i=0; i<level; i++) num *= omp_get_team_size(i+1);
			num = omp_get_max_threads()/num;
#endif
#pragma omp parallel for schedule(guided) num_threads(num) // copyin(console::tab_prefix,console::level)
			for(uint32_t i=0; i<imax; i++){
				uint32_t jmax = Iptr[i+1];
				for(uint32_t j=Iptr[i]; j<jmax; j++) y[i] += z*x[J[j]]*v[j];
			}
		}
		
		else{
			// then multiply by the upper triangle
			for(uint32_t i=0; i<imax; i++){
				uint32_t jmax = Iptr[i+1];
				for(uint32_t j=Iptr[i]; j<jmax; j++) y[i] += z*x[J[j]]*v[j];
			}
			
			// then by the lower triangle
			for(uint32_t i=0; i<imax; i++){
				uint32_t jmax = Iptr[i+1];
				for(uint32_t j=Iptr[i]; j<jmax; j++) y[J[j]] += z*x[i]*conjugate(v[j]);
			}
		}
	}
	
	void clear(){
		erase(J);
		erase(Iptr);
		erase(v);
		erase(diag);
	}
	
	
  /**
   Prints the matrix if of modest size
   */
	friend std::ostream & operator<<(std::ostream &s, const CSR_hermitian &m)
	{
		if(m.Iptr.size()>100) return s;
		s << "diagonal elements:\n";
		for(size_t i=0; i<m.diag.size(); i++) s << i << '\t' << i << '\t' << m.diag[i] << endl;
		
		s << "upper triangle:\n";
		size_t imax = m.Iptr.size()-1;
		for(size_t i=0; i<imax; i++){
			for(size_t j=m.Iptr[i]; j<m.Iptr[i+1]; j++){
				s << i << '\t' << m.J[j] << '\t' << m.v[j] << endl;
			}
		}
		
		return s;
	}
	
};


#endif
