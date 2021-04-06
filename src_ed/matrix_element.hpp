#ifndef matrix_element_h
#define matrix_element_h

#include <iostream>
#include <string>

using namespace std;


//! used to store a matrix element, location and value
template<typename T>
struct matrix_element
{
	size_t r;
	size_t c;
  T v;
	
	matrix_element() : r(0), c(0), v(0) {}
	
	matrix_element(size_t _r, size_t _c, T _v) : r(_r), c(_c), v(_v) {}
  
	friend std::ostream & operator<<(std::ostream &s, const matrix_element &e)
	{
    s << '(' << e.r+1 << ',' << e.c+1 << ") : " << e.v;
		return s;
	}
	
  string str() const{
		ostringstream sout;
		sout << *this;
		return sout.str();
	}

	friend inline bool operator<(const matrix_element &x, const matrix_element &y){
		if(x.r < y.r) return true;
		else if(x.r > y.r) return false;
		else if(x.c < y.c) return true;
		else return false;
	}
	
	
};




struct diagonal_matrix_element{
  uint32_t r;
  double v;
  
  diagonal_matrix_element(uint32_t _r, double _v) : r(_r), v(_v) {}

  friend inline bool operator<(const diagonal_matrix_element &x, const diagonal_matrix_element &y){
    if(x.r < y.r) return true;
    else return false;
  }
} ;


#define MAXORB 64
template<typename op_field>
struct interaction_matrix_element
{
  size_t i;
  size_t j;
  size_t k;
  size_t l;
  op_field v;
  
  interaction_matrix_element() : i(0), j(0), k(0), l(0), v(0) {}
  
  interaction_matrix_element(size_t _i, size_t _j, size_t _k, size_t _l, op_field _v) : i(_i), j(_j), k(_k), l(_l), v(_v) {}
  
  inline size_t rank() const
  {
    return i + MAXORB*(j + MAXORB*(k + MAXORB*l));
  }
  
  friend std::ostream & operator<<(std::ostream &s, const interaction_matrix_element &e)
  {
    s << '(' << e.i+1 << ',' << e.j+1 << "|" << e.k+1 << ',' << e.l+1 << ") : " << e.v;
    return s;
  }
  
  string str() const{
    ostringstream sout;
    sout << *this;
    return sout.str();
  }
};

template<typename op_field>
bool operator<(const interaction_matrix_element<op_field> &x, const interaction_matrix_element<op_field> &y){
  return x.rank() < y.rank();
}

#endif
