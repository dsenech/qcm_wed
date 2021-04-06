/*
 Class for a Hamiltonian of a given model_instance in a given Hilbert space sector
*/

#ifndef Hamiltonian_h
#define Hamiltonian_h

#include "model.hpp"
#include "state.hpp"
#include "Hermitian_operator.hpp"
#include "continued_fraction.hpp"
#include "Q_matrix.hpp"
#include "Lanczos.hpp"
#include "Davidson.hpp"
#include "CSR_hermitian.hpp"

extern double max_gap;
extern std::normal_distribution<double> normal_dis;

#define DEGENERACY_LIMIT 0.001


//! Represents the Hamiltonian of the system, for a given value of the parameters and a given sector of the HS
template<typename HilbertField>
struct Hamiltonian
{
  CSR_hermitian<HilbertField> csr; //!< CSR matrix of the Hamiltonian
  H_FORMAT format; //!< format of the Hamiltonian (see H_FORMAT enumerator above)
  map<shared_ptr<Hermitian_operator>, double> ops; //!< correpondence between terms in H and their coefficients
  map<shared_ptr<HS_Hermitian_operator>, double> sparse_ops; //!< correpondence between terms in H and their coefficients
  matrix<HilbertField> H_dense; //!< dense form of the Hamiltonian (used when the dimension is small)
  sector sec; //!< sector of the HS on which the Hamiltonian acts
  shared_ptr<ED_mixed_basis> B; //!<  pointer to basis of the space on which the Hamiltonian is defined
  shared_ptr<model> the_model; //!< backtrace to the cluster model
  size_t dim; //!< dimension of the HS sector on which the Hamiltonian acts
  vector<double> alpha; //!< main diagonal of the projected Hamiltonian in the Lanczos basis
  vector<double> beta; //!< second diagonal of the projected Hamiltonian in the Lanczos basis

  Hamiltonian(shared_ptr<model> the_model, const map<string, double> &value, sector _sec);
  map<shared_ptr<HS_Hermitian_operator>, double> HS_ops_map(const map<string, double> &value);
  map<shared_ptr<Hermitian_operator>, double> ops_map(const map<string, double> &value);
  void dense_form();
  void mult_add(vector<HilbertField> &x, vector<HilbertField> &y);
  vector<shared_ptr<state<HilbertField>>> states(double& GS_energy);
  double GS_energy();
  void diag(vector<double> &d);
  Q_matrix<HilbertField> build_Q_matrix(vector<vector<HilbertField>> &phi);
  void print(ostream& fout);
};

//==============================================================================
// implementation

/**
 constructor
 @param _B : pointer to Hilbert space basis
 @param _sparse_ops : map providing values to the different terms of the Hamiltonian
 */
template<typename HilbertField>
Hamiltonian<HilbertField>::Hamiltonian(shared_ptr<model> _the_model, const map<string, double> &value, sector _sec)
: the_model(_the_model), sec(_sec), format(H_format_csr)
{
  if(the_model->is_factorized){
    dim = the_model->factorized_basis.at(_sec)->dim;
    format = H_format_factorized;
  }
  else{
    B = the_model->basis.at(_sec);
    dim = B->dim;
    int dim_max_full = global_int("max_dim_full");
    if(dim < dim_max_full) format = H_format_dense;
    else format = Hamiltonian_format;
    if(format == H_format_factorized) format = H_format_csr;
  }

  if(dim == 0) return;

  if(format != H_format_onthefly) sparse_ops = HS_ops_map(value);

  switch(format){
    case H_format_ops:
    case H_format_factorized:
      break;

    case H_format_dense:
      dense_form();
      if(global_bool("print_Hamiltonian")) print(cout);
      break;

    case H_format_onthefly:
      ops = ops_map(value);
      break;

    case H_format_csr:
      int num=1;
      #ifdef _OPENMP
        num = omp_get_max_threads()/omp_get_num_threads();
      #endif  
      format = H_format_csr;
      map<index_pair,HilbertField> E;
      if(global_bool("CSR_sym_store") and num>1) csr.sym_store = true; // set up the CSR format for openMP parallelization
      else csr.sym_store = false;
      csr.diag.assign(dim, 0.0);
      if(csr.sym_store)
        console::message(5, "constructing the CSR Hamiltonian (stored symmetrically for openMP with "+to_string(num)+" threads)...");
      else console::message(5, "constructing the CSR Hamiltonian...");
      for(auto& h : sparse_ops){
        h.first->CSR_map(E, csr.diag, h.second, csr.sym_store);
      }
      size_t row = 0;
      size_t count=0;
      csr.Iptr.reserve(dim/2);
      csr.J.reserve(E.size());
      csr.v.reserve(E.size());
      csr.Iptr.push_back(0);
      for(auto &x : E){
        if(x.first.r != row){
          for(size_t i=row; i<x.first.r; i++) csr.Iptr.push_back(count);
          row = x.first.r;
        }
        csr.J.push_back(x.first.c);
        csr.v.push_back(x.second);
        count++;
      }
      csr.Iptr.push_back(count);
      break;

  }
}



/**
 builds HS_operators as needed
 */
template<typename HilbertField>
map<shared_ptr<HS_Hermitian_operator>, double> Hamiltonian<HilbertField>::HS_ops_map(const map<string, double> &value)
{
  map<shared_ptr<HS_Hermitian_operator>, double> sparse_ops;
  for(auto& x : value){
    Hermitian_operator& op = *the_model->term.at(x.first);
    if(op.HS_operator.find(sec) == op.HS_operator.end()){
      qcm_ED_throw("HS operator "+op.name+" for sector "+sec.name()+" not found!");
    }
    sparse_ops[op.HS_operator.at(sec)] = value.at(x.first);
  }
  return sparse_ops;
}



/**
 returns a map of Hermitian operators to value
 */
template<typename HilbertField>
map<shared_ptr<Hermitian_operator>, double> Hamiltonian<HilbertField>::ops_map(const map<string, double> &value)
{
  map<shared_ptr<Hermitian_operator>, double> ops;
  for(auto& x : value){
    ops[the_model->term.at(x.first)] = value.at(x.first);
  }
  return ops;
}



/**
 Puts a dense form of the Hamiltonian matrix in the matrix H_dense
 Mostly for debugging purposes, but also for small dimensions.
 */
template<typename HilbertField>
void Hamiltonian<HilbertField>::dense_form()
{
  if(H_dense.v.size()) return;
  H_dense.set_size(dim);
  if(format == H_format_onthefly){
    for(auto& i : ops) i.first->dense_form(H_dense, i.second, sec);
  }
  else{
    for(auto& i : sparse_ops) i.first->dense_form(H_dense, i.second);
  }
  return;
}




/**
 Applies the Hamiltonian: y = y +H.x
 @param y vector to which H.x is added to
 @param x input vector
 */
template<typename HilbertField>
void Hamiltonian<HilbertField>::mult_add(vector<HilbertField> &x, vector<HilbertField> &y)
{
  switch(format){
    case H_format_ops:
    case H_format_factorized:
      // applies all of the terms in H in turn
      for(auto& h : sparse_ops){
        h.first->multiply_add(x, y, h.second);
      }
      break;

    case H_format_csr:
      csr.apply(x, y); // applies the CSR matrix 
      break;
      
    case H_format_dense:
      H_dense.apply_add(x, y);
      break;
      
    case H_format_onthefly:
      for(auto& h : ops){
        h.first->multiply_add_OTF(x, y, h.second, B);
      }
      break;
  }
}




/**
 Applies the Lanczos algorithm for the lowest energy.
 @param alpha : first diagonal of the tridiagonal representation (returned by reference)
 @param beta : second diagonal of the tridiagonal representation (returned by reference)
returns the GS energy.
 */
template<typename HilbertField>
double Hamiltonian<HilbertField>::GS_energy()
{
  if(format == H_format_dense){
    vector<double> evalues(dim);
    H_dense.eigenvalues(evalues);
    return evalues[0];
  }
  else{
    vector<double> energy;
    size_t niter = 0;
    vector<HilbertField> x(dim);
    random(x, normal_dis);
    LanczosEigenvalue(*this, x, alpha, beta, energy, niter);
	  return energy[0];
  }
}





/**
 Applies the Lanczos or Davidson-Liu algorithm for the lowest-energy states and energies.
 returns a vector of pointers to states, to be treated by the parent model_instance.
 @param GS_energy : current ground-state energy of the model_instance, to be updated
 */
template<typename HilbertField>
vector<shared_ptr<state<HilbertField>>> Hamiltonian<HilbertField>::states(double& GS_energy)
{
  vector<shared_ptr<state<HilbertField>>> low_energy_states;
  
  if(format == H_format_dense){
    vector<double> evalues(dim);
    matrix<HilbertField> U((int)dim);
    H_dense.eigensystem(evalues,U);
    if(evalues[0] < GS_energy) GS_energy = evalues[0];
    for(size_t i=0; i<evalues.size(); i++){
      if(evalues[i]-GS_energy > max_gap) break;
      auto gs = make_shared<state<HilbertField>>(sec,dim);
      U.extract_column(i,gs->psi);
      gs->energy = evalues[i];
      low_energy_states.push_back(gs);
    }
  }
  else{
    vector<double> evalues;
    vector<vector<HilbertField> > evectors;
    size_t Davidson_states = global_int("Davidson_states");
    if(Davidson_states > 1){
      Davidson(*this, dim, Davidson_states, evalues, evectors, global_double("accur_Davidson"));
      if(evalues[0] < GS_energy) GS_energy = evalues[0];
      if(evalues.back()-GS_energy < max_gap){
        console::message(1,"ED WARNING! : not enough Davidson states (" + to_string(Davidson_states) + ") in sector " + sec.name());
      }
    }
    else{
      evalues.resize(1);
      evectors.resize(1);
      evectors[0].resize(dim);
      Lanczos(*this, dim, evalues[0], evectors[0]);
      if(evalues[0] < GS_energy) GS_energy = evalues[0];
    }
    for(size_t i=0; i<evectors.size(); i++){
      if(evalues[i]-GS_energy > max_gap) continue;
      auto gs = make_shared<state<HilbertField>>(sec,dim);
      gs->energy = evalues[i];
      gs->psi = evectors[i];
      low_energy_states.push_back(gs);
    }
  }
  return low_energy_states;
}




/**
 provides the diagonal d of H
 Used by the Davidson method
 @param d the diagonal of H (pre-allocated)
 */
template<typename HilbertField>
void Hamiltonian<HilbertField>::diag(vector<double> &d){
  for(auto& h : sparse_ops) h.first->diag(d, h.second);
}




/**
 Constructs the Q_matrix (Lehmann representation) from the Band Lanczos method,
 or full diagonalization if the dimension is small enough.
 @param phi the initial vectors
 */
template<typename HilbertField>
Q_matrix<HilbertField> Hamiltonian<HilbertField>::build_Q_matrix(vector<vector<HilbertField>> &phi)
{
  if(dim == 0 or phi.size()==0){
    return Q_matrix<HilbertField>(0,0);
  }
  
  //-----------------------------------------------------------------------------
  // case of small dimensions : the Hamiltonian is already stored in dense form
  
  if(format == H_format_dense){
    console::message(5,"Q_matrix : full diagonalization");
    Q_matrix<HilbertField> Q(phi.size(), dim);
    vector<HilbertField> y(dim);
    matrix<HilbertField> U(H_dense);
    H_dense.eigensystem(Q.e, U);
    assert(U.is_unitary(1e-6));
    for(size_t i=0; i<phi.size(); ++i){
      to_zero(y);
      U.left_apply_add(phi[i],y); // does y = phi[i] . U
      Q.v.insert_row(i,y); // inserts y as the ith row of Q
    }
    return Q;
  }
  
  //-----------------------------------------------------------------------------
  // Setting the maximum number of iterations
  
  int max_banditer = (int)(14*phi.size()*log(1.0*dim));
  int M = max_banditer; // essai
  assert(M>1);
  
  vector<double> eval; // eigenvalues of the reduced Hamiltonian
  matrix<HilbertField> U;  // eigenvectors of the reduced Hamiltonian
  matrix<HilbertField> P; // matrix of inner products <b[i]|v[j]>
  
  if(BandLanczos(*this, phi, eval, U, P, M)){
    Q_matrix<HilbertField> Q(phi.size(),M);
    if(Q.M > 0){
      Q.e = eval; 
      Q.v.product(P,U,phi.size()); //  tempo
      // Q.v.product(P,U; 
    }
    return Q;
  }
  else return Q_matrix<HilbertField>(phi.size(),0);

}




/**
 Prints a dense form of the Hamiltonian
 Mostly for debugging purposes on very small systems
 */
template<typename HilbertField>
void Hamiltonian<HilbertField>::print(ostream& fout)
{
  if(H_dense.v.size() == 0) return;
  console::banner('~', "Hamiltonian", fout);
  fout << *B;
  fout << "Hamiltonian (dense form):\n";
  fout << H_dense;
}

#endif
