/*
Base abstract class for a Hamiltonian of a given model_instance in a given
Hilbert space sector that serve as a proxy for the multiple implementation of
the Hamiltonian (Dense, legacy CSR, PETSc, Eigen, ...)
Also provide default method that could be overwrite by implementation
*/

#ifndef Hamiltonian_base_h
#define Hamiltonian_base_h

#include "model.hpp"
#include "state.hpp"
#include "Hermitian_operator.hpp"
#include "Q_matrix.hpp"
#include "Lanczos.hpp"

extern double max_gap;
extern std::normal_distribution<double> normal_dis;

template<typename HilbertField>
class Hamiltonian
{
    //legacy, unused attribute
    //H_FORMAT format; //now implicit with the implementation
    //CSR_hermitian<HilbertField> csr; //private property of the implementation
    //matrix<HilbertField> H_dense; 
    
    public:
    
        size_t dim; //!< dimension of the HS sector on which the Hamiltonian acts
        //this should not be belonging to the Hamiltonian class
        sector sec; //!< sector of the HS on which the Hamiltonian acts
        shared_ptr<ED_mixed_basis> B; //!<  pointer to basis of the space on which the Hamiltonian is defined
        shared_ptr<model> the_model; //!< backtrace to the cluster model
        
        Hamiltonian() {};
    
        virtual void mult_add(vector<HilbertField> &x, vector<HilbertField> &y) {};
        virtual void diag(vector<double> &d) {};
        virtual double GS_energy();
        virtual vector<shared_ptr<state<HilbertField>>> states(double& GS_energy);
        virtual void print(ostream& fout) {};
        //virtual matrix<HilbertField> to_dense();
        
        virtual Q_matrix<HilbertField> build_Q_matrix(vector<vector<HilbertField>> &phi);
        
    
    private:
    
        vector<double> alpha; //!< main diagonal of the projected Hamiltonian in the Lanczos basis
        vector<double> beta; //!< second diagonal of the projected Hamiltonian in the Lanczos basis


    // ??? what's below ?
    
    map<shared_ptr<HS_Hermitian_operator>, double> HS_ops_map(const map<string, double> &value);
    map<shared_ptr<Hermitian_operator>, double> ops_map(const map<string, double> &value);
    //void dense_form();
    
};


// ######################################
// Common method for every implementation
// ######################################

/**
 Applies the Lanczos algorithm for the lowest energy.
 Common with Hamiltonian in CSR and Factorized format.
 @param alpha : first diagonal of the tridiagonal representation (returned by reference)
 @param beta : second diagonal of the tridiagonal representation (returned by reference)
returns the GS energy.
 */
template<typename HilbertField>
double Hamiltonian<HilbertField>::GS_energy()
{
    vector<double> energy;
    size_t niter = 0;
    vector<HilbertField> x(dim);
    random(x, normal_dis);
    LanczosEigenvalue(*this, x, alpha, beta, energy, niter,  global_bool("verb_ED"));
    return energy[0];
}


/**
 Applies the Lanczos or Davidson-Liu algorithm for the lowest-energy states and energies.
 Common with Hamiltonian in CSR and Factorized format.
 returns a vector of pointers to states, to be treated by the parent model_instance.
 @param GS_energy : current ground-state energy of the model_instance, to be updated
 */
template<typename HilbertField>
vector<shared_ptr<state<HilbertField>>> Hamiltonian<HilbertField>::states(double& GS_energy)
{
    vector<shared_ptr<state<HilbertField>>> low_energy_states;
  
    vector<double> evalues;
    vector<vector<HilbertField> > evectors;
    size_t Davidson_states = global_int("Davidson_states");
    if(Davidson_states > 1) {
        Davidson(*this, dim, Davidson_states, evalues, evectors, global_double("accur_Davidson"),  global_bool("verb_ED"));
        if(evalues[0] < GS_energy) GS_energy = evalues[0];
        if(evalues.back()-GS_energy < max_gap and global_bool("verb_warning")) {
            cout << "ED WARNING! : not enough Davidson states (" << Davidson_states << ") in sector " << sec.name() << endl;
        }
    }
    else {
        evalues.resize(1);
        evectors.resize(1);
        evectors[0].resize(dim);
        Lanczos(*this, dim, evalues[0], evectors[0],  global_bool("verb_ED"));
        if(evalues[0] < GS_energy) GS_energy = evalues[0];
    }
    for(size_t i=0; i<evectors.size(); i++){
        if(evalues[i]-GS_energy > max_gap) continue;
        auto gs = make_shared<state<HilbertField>>(sec,dim);
        gs->energy = evalues[i];
        gs->psi = evectors[i];
        low_energy_states.push_back(gs);
    }
    return low_energy_states;
}


/**
 Constructs the Q_matrix (Lehmann representation) from the Band Lanczos method,
 or full diagonalization if the dimension is small enough.
 @param phi the initial vectors
 */
template<typename HilbertField>
Q_matrix<HilbertField> Hamiltonian<HilbertField>::build_Q_matrix(
    vector<vector<HilbertField>> &phi
) {
    if(this->dim == 0 or phi.size()==0){
        return Q_matrix<HilbertField>(0,0);
    }
  
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
