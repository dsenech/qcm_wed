#ifndef general_interaction_operator_h
#define general_interaction_operator_h

#include "Hermitian_operator.hpp"

template<typename op_field>
struct general_interaction_operator : Hermitian_operator
{
  vector<interaction_matrix_element<op_field>> elements; //!< matrix elements

  general_interaction_operator(const string &_name, shared_ptr<model> the_model, const vector<interaction_matrix_element<op_field>>& _elements);
  void check_spin_symmetry();
  void set_target(vector<bool> &in_bath);
  shared_ptr<HS_Hermitian_operator> build_HS_operator(shared_ptr<ED_basis> B, bool complex_Hilbert_space);
  void set_hopping_matrix(double value, matrix<double>& tc, bool spin_down, int sys_mixing){}
  void set_hopping_matrix(double value, matrix<Complex>& tc, bool spin_down, int sys_mixing){}
  double average_from_GF(matrix<Complex>& Gave, bool spin_down){return 0.0;}
  void diag(vector<double> &d, double z);
  void print(ostream& fout);
  vector<matrix_element<Complex>> matrix_elements();
  string type() {return string("general_interaction");}
  void multiply_add_OTF(const vector<double> &x, vector<double> &y, double z, shared_ptr<ED_basis> B);
  void multiply_add_OTF(const vector<Complex> &x, vector<Complex> &y, double z, shared_ptr<ED_basis> B);
};


/**
 returns a list of complexified matrix elements
 */
template<typename op_field>
vector<matrix_element<Complex>> general_interaction_operator<op_field>::matrix_elements()
{
  return vector<matrix_element<Complex>>(0);
}


template<typename op_field_HS, typename op_field>
struct HS_general_interaction_operator : HS_nondiagonal_operator<op_field_HS>
{
  // general_interaction_operator<op_field>& op;
  HS_general_interaction_operator(const general_interaction_operator<op_field>& _op, shared_ptr<ED_basis> _B);
};

/**
 Specialization of the above for a real Hilbert space
 */
template<typename op_field>
struct HS_general_interaction_operator<double, op_field> : HS_nondiagonal_operator<double>
{
  // general_interaction_operator<op_field>& op;
  
  HS_general_interaction_operator(const general_interaction_operator<op_field>& _op, shared_ptr<ED_basis> _B)
  : HS_nondiagonal_operator(_B, _op.name)//, op(_op)
  {
    // to be completed
  }
};


/**
 Specialization of the above for a complex Hilbert space
 */
template<typename op_field>
struct HS_general_interaction_operator<Complex, op_field> : HS_nondiagonal_operator<Complex>
{
  // general_interaction_operator<op_field>& op;
  
  HS_general_interaction_operator(const general_interaction_operator<op_field>& _op, shared_ptr<ED_basis> _B)
  : HS_nondiagonal_operator(_B, _op.name)//, op(_op)
  {
    // to be completed
  }
};

//==============================================================================
// implementation of general_interaction_operator




/**
 Constructor from name and matrix elements
 @param _name   name of the operator
 @param _the_model   model
 @param _elements   nonzero one-body matrix elements
 */
template<typename op_field>
general_interaction_operator<op_field>::general_interaction_operator(const string &_name, shared_ptr<model> _the_model, const vector<interaction_matrix_element<op_field>>& _elements)
: Hermitian_operator(_name, _the_model)
{
  // checks not done here.
}




/**
 determines whether the operator is symmetric under the exchange of up and down spins
 */
template<typename op_field>
void general_interaction_operator<op_field>::check_spin_symmetry()
{
}





/**
 set the target of an operator
 1 : cluster
 2 : bath only
 3 : hybridization
 @param in_bath vector of bool defining the status of each site
 */
template<typename op_field>
void general_interaction_operator<op_field>::set_target(vector<bool> &in_bath){
  this->target = 1;
}




/**
 returns a pointer to, and constructs the associated HS operator in the sector with basis B.
 */
template<typename op_field>
shared_ptr<HS_Hermitian_operator>  general_interaction_operator<op_field>::build_HS_operator(shared_ptr<ED_basis> B, bool complex_Hilbert_space)
{
  if(complex_Hilbert_space) return make_shared<HS_general_interaction_operator<Complex, op_field>>(*this, B);
  else return make_shared<HS_general_interaction_operator<double, op_field>>(*this, B);
}




/**
 prints definition to a file
 @param fout output stream
 */
template<typename op_field>
void general_interaction_operator<op_field>::print(ostream& fout)
{
  fout << "\ninteraction operator " << name << "\t (target " << target << ")" << endl;
  for(auto& x : elements) fout << x << endl;
}

template<typename op_field>
void general_interaction_operator<op_field>::multiply_add_OTF(const vector<double> &x, vector<double> &y, double z, shared_ptr<ED_basis> B)
{
  qcm_ED_throw("on the fly computation imnpossible with a general interaction operator");
}

template<typename op_field>
void general_interaction_operator<op_field>::multiply_add_OTF(const vector<Complex> &x, vector<Complex> &y, double z, shared_ptr<ED_basis> B)
{
  qcm_ED_throw("on the fly computation imnpossible with a general interaction operator");
}

#endif
