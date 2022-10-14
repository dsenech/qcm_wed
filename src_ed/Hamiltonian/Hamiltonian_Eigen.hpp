/*
 Implementation of the Hamiltonian in CSR format with Eigen
*/

#ifndef Hamiltonian_eigen
#define Hamiltonian_eigen

#include "Hamiltonian_base.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>


template<typename HilbertField>
class Hamiltonian_Eigen : public Hamiltonian<HilbertField>
{
    public:
    
        	
        Eigen::SparseMatrix<HilbertField,Eigen::RowMajor> H_eigen;
        
        Hamiltonian_Eigen(
            shared_ptr<model> the_model, 
            const map<string, double> &value,
            sector _sec
        );
        void mult_add(vector<HilbertField> &x, vector<HilbertField> &y);
        void diag(vector<double> &d);
        Q_matrix<HilbertField> build_Q_matrix(vector<vector<HilbertField>> &phi);
        //vector<shared_ptr<state<HilbertField>>> states(double& GS_energy);
        
    private:
        map<shared_ptr<HS_Hermitian_operator>, double> sparse_ops; //!< correpondence between terms in H and their coefficients
        void HS_ops_map(const map<string, double> &value);

};


/**
 constructor
 */
template<typename HilbertField>
Hamiltonian_Eigen<HilbertField>::Hamiltonian_Eigen(
    shared_ptr<model> _the_model,
    const map<string, double> &value, 
    sector _sec
) {
    this->the_model = _the_model;
    this->sec = _sec;
    this->B = _the_model->provide_basis(_sec);
    this->dim = this->B->dim;
    if(this->dim == 0) return;

    HS_ops_map(value);
    
    //create the element
    map<index_pair,HilbertField> E;
    bool sym_store = true;
    std::vector<double> diag(this->dim);
    for(auto& h : sparse_ops){
        h.first->CSR_map(E, diag, h.second, sym_store);
    }
    
    std::vector< Eigen::Triplet<HilbertField> > tripletList;
    tripletList.reserve(E.size()+this->dim);
    //set up diagonal elements
    for (size_t i=0; i<this->dim; i++) {
        Eigen::Triplet<HilbertField> T(i,i,diag[i]);
        tripletList.push_back(T);
    }
    diag.resize(0); //clear vector
    //set non-diag element
    for(auto &x : E) {
        Eigen::Triplet<HilbertField> T(x.first.r,x.first.c,x.second);
        tripletList.push_back(T);
    }
    //create matrix
    H_eigen.resize(this->dim,this->dim);
    H_eigen.setFromTriplets(tripletList.begin(), tripletList.end());
}


/**
 Applies the Hamiltonian: y = y +H.x
 @param y vector to which H.x is added to
 @param x input vector
 */
template<typename HilbertField>
void Hamiltonian_Eigen<HilbertField>::mult_add(
    vector<HilbertField> &x, 
    vector<HilbertField> &y
) {
     //TODO see some optimization
     Eigen::Map< Eigen::Matrix<HilbertField,Eigen::Dynamic,1> > xe(x.data(), x.size());
     Eigen::Map< Eigen::Matrix<HilbertField,Eigen::Dynamic,1> > ye(y.data(), y.size());
     ye += H_eigen*xe; //this change value of ye in place, so for y
}


/**
 provides the diagonal d of H
 Used by the Davidson method
 @param d the diagonal of H (pre-allocated)
 */
template<typename HilbertField>
void Hamiltonian_Eigen<HilbertField>::diag(vector<double> &d){
    //for (size_t i=0; i<d.size(); i++) d[i] = H_eigen.diagonal()[i];
    for(auto& h : sparse_ops) h.first->diag(d, h.second);
}

/**
 builds HS_operators as needed
 */
template<typename HilbertField>
void Hamiltonian_Eigen<HilbertField>::HS_ops_map(const map<string, double> &value)
{
    bool is_complex = false;
    if(typeid(HilbertField) == typeid(Complex)) is_complex = true;
    for(auto& x : value){
        Hermitian_operator& op = *this->the_model->term.at(x.first);
        if(op.HS_operator.find(this->sec) == op.HS_operator.end()){
            op.HS_operator[this->sec] = op.build_HS_operator(this->sec, is_complex); // ***TEMPO***
        }
        sparse_ops[op.HS_operator.at(this->sec)] = value.at(x.first);
    }
}


/**
 Constructs the Q_matrix (Lehmann representation) from the Band Lanczos method,
 or full diagonalization if the dimension is small enough.
 @param phi the initial vectors
 */
template<typename HilbertField>
Q_matrix<HilbertField> Hamiltonian_Eigen<HilbertField>::build_Q_matrix(
    vector<vector<HilbertField>> &phi
) {
    if(this->dim == 0 or phi.size()==0){
        return Q_matrix<HilbertField>(0,0);
    }
  
    //-----------------------------------------------------------------------------
    // Setting the maximum number of iterations
  
    int max_banditer = (int)(14*phi.size()*log(1.0*this->dim));
    int M = max_banditer; // essai
    assert(M>1);

    vector<double> eval; // eigenvalues of the reduced Hamiltonian
    matrix<HilbertField> U;  // eigenvectors of the reduced Hamiltonian
    matrix<HilbertField> P; // matrix of inner products <b[i]|v[j]>
  
    if(BandLanczos(*this, phi, eval, U, P, M,  global_bool("verb_ED"))){
        Q_matrix<HilbertField> Q(phi.size(),M);
        if(Q.M > 0) {
            Q.e = eval; 
            Q.v.product(P,U,phi.size()); //  tempo
        }
        return Q;
    }
    else return Q_matrix<HilbertField>(phi.size(),0);
}

#endif
