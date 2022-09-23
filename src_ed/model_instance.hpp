#ifndef model_instance_h
#define model_instance_h

#include <fstream>
#include <memory>
#include "model_instance_base.hpp"
#include "Q_matrix_set.hpp"
#include "continued_fraction_set.hpp"
#include "ED_basis.hpp"
#include "binary_state.hpp"

#include "Hamiltonian/Hamiltonian_base.hpp"
#include "Hamiltonian/Hamiltonian_Dense.hpp"
#include "Hamiltonian/Hamiltonian_CSR.hpp"
#include "Hamiltonian/Hamiltonian_OnTheFly.hpp"
#include "Hamiltonian/Hamiltonian_Factorized.hpp"

extern double max_gap;

void polynomial_fit(
  vector<double> &xa, //!< array of abcissas
  vector<double> &ya, //!< array of values
  const double x,	//!< projected abcissa
  double &y,	//!< requested extrapolated value
  double &dy	//!< requested error on extrapolated value
);



//! template class for an instance of a model
/**
 a model_instance is of type <Complex> if one of the Hermitian operators in the model is complex,
 or if the symmetry group has complex representations. Either of these results in a complex
 Hilbert space.
 Otherwise it is of type <double>
 Note that some of the Hermitian operators could be of type <double> even if the Hilbert space is
 complex.
 */
template<typename HilbertField>
struct model_instance : model_instance_base
{
  // members
  matrix<HilbertField> tc, tcb, tb; //! one-body matrices for cluster, cluster-bath and bath
  matrix<HilbertField> tc_down, tcb_down, tb_down; //! one-body matrices for cluster, cluster-bath and bath (spin down)
  matrix<HilbertField> tcb_nd, tb_nd; //! non diagonalized versions of these matrices (for debugging/printing)
  matrix<HilbertField> tcb_nd_down, tb_nd_down; //! same, for the spin-down version when mixing = 4
  set<shared_ptr<state<HilbertField>>> states; //!< set of states forming the density matrix

  model_instance(size_t _label, shared_ptr<model> _the_model, const map<string,double> _value, const string &_sectors);
  Hamiltonian<HilbertField> create_hamiltonian(
    shared_ptr<model> the_model,
    const map<string, double> &value, 
    sector s
  );

  void build_cf(state<HilbertField> &Omega, bool spin_down);
  void build_qmatrix(state<HilbertField> &Omega, bool spin_down);
  void clear_states();
  void compute_weights();
  void insert_state(shared_ptr<state<HilbertField>> S);
  void set_hopping_matrix(bool spin_down);
  string GS_string() const;
  void build_bases_and_HS_operators(const sector& GS_sector, bool spin_down);
  
  // realization of virtual base class methods
  pair<double, string> low_energy_states();
  pair<double, double> cluster_averages(shared_ptr<Hermitian_operator> h);
  void Green_function_solve();
  pair<double, string> one_body_solve();
  matrix<Complex>  Green_function(const Complex &z, bool spin_down, bool blocks);
  matrix<Complex>  Green_function_average(bool spin_down);
  matrix<Complex>  self_energy(const Complex &z, bool spin_down);
  matrix<Complex>  hopping_matrix(bool spin_down);
  matrix<Complex>  hopping_matrix_full(bool spin_down);
  vector<tuple<int,int,double>> interactions();
  matrix<Complex>  hybridization_function(Complex w, bool spin_down);
  vector<Complex>  susceptibility(shared_ptr<Hermitian_operator> h, const vector<Complex> &w);
  vector<pair<double,double>> susceptibility_poles(shared_ptr<Hermitian_operator> h);
  double fidelity(model_instance<HilbertField>& inst);
  void print(ostream& fout);
  void write(ostream& fout);
  void read(istream& fin);
  double tr_sigma_inf();
  pair<vector<double>, vector<complex<double>>> qmatrix(bool spin_down);
  pair<vector<double>, vector<complex<double>>> hybridization(bool spin_down);
  void print_wavefunction(ostream& fout);
  pair<matrix<Complex>, vector<uint64_t>>  density_matrix_mixed(vector<int> sites);
  pair<matrix<Complex>, vector<uint64_t>>  density_matrix_factorized(vector<int> sites);
};



//==============================================================================
// implementation of model_instance



/**
 constructor
 */
template<typename HilbertField>
model_instance<HilbertField>::model_instance(size_t _label, shared_ptr<model> _the_model, const map<string,double> _value, const string &_sectors)
: model_instance_base(_label, _the_model, _value, _sectors)

{
  if(typeid(HilbertField) == typeid(Complex)) complex_Hilbert = true;
  else complex_Hilbert = false;

  // allocating matrices
  tc.set_size(n_mixed*the_model->n_sites);
  tb.set_size(n_mixed*the_model->n_bath);
  tcb.set_size(n_mixed*the_model->n_sites, n_mixed*the_model->n_bath);
  if(mixing&HS_mixing::up_down){
    tc_down.set_size(n_mixed*the_model->n_sites);
    tb_down.set_size(n_mixed*the_model->n_bath);
    tcb_down.set_size(n_mixed*the_model->n_sites,n_mixed*the_model->n_bath);
  }
  try{
    set_hopping_matrix(false);
  } catch(const string& s) {qcm_ED_catch(s);}
}





/**
 returns the GS average of a Hermitian operator op
 */
template<typename HilbertField>
pair<double, double> model_instance<HilbertField>::cluster_averages(shared_ptr<Hermitian_operator> h)
{
  if(value.find(h->name)==value.end()) qcm_ED_throw("operator "+h->name+" is not activated in instance "+to_string(label));
  double ave=0.0, var=0.0;
  if(!is_correlated or gf_read){
    ave = h->uncorrelated_average(M,false);
    if(mixing&HS_mixing::up_down) ave += h->uncorrelated_average(M_down,true);
    if(mixing == HS_mixing::normal) ave *= 2.0;
    if(mixing == HS_mixing::anomalous) ave += h->nambu_correction;
    if(mixing == HS_mixing::full) {
      ave += h->nambu_correction_full;
      ave *= 0.5;
    }
  }
  else{
    for(auto& gs : states) h->expectation_value(*gs, ave, var);
    var -= ave*ave;
  }
  if(h->target == 1){
    ave *= h->norm;
    var *= h->norm*h->norm;
  }
  return {ave, var};
}



/**
 erases all the states vectors in the density matrix
 */
template<typename HilbertField>
void model_instance<HilbertField>::clear_states(){
  for(auto& x : states) erase(x->psi);
  gs_solved = false;
}


/**
 Create Hamiltonian in the right format
 */
template<typename HilbertField>
Hamiltonian<HilbertField> model_instance<HilbertField>::create_hamiltonian(
    shared_ptr<model> the_model,
    const map<string, double> &value, 
    sector s
) {
    Hamiltonian<HilbertField> H;
    //enforced Hamiltonian format
    if(the_model->is_factorized) {
        H = Hamiltonian_Factorized<HilbertField>(the_model, value, s);
    }
    else if (the_model->provide_basis(s)->dim < global_int("max_dim_full")) {
        H = Hamiltonian_Dense<HilbertField>(the_model, value, s);
    }
    else {
        switch(Hamiltonian_format) {
            case H_format_csr:
                H = Hamiltonian_CSR<HilbertField>(the_model, value, s);
                break;
            case H_format_dense:
                H = Hamiltonian_Dense<HilbertField>(the_model, value, s);
                break;
            case H_format_onthefly:
                H = Hamiltonian_OnTheFly<HilbertField>(the_model, value, s);
                break;
        }
    }
    return H;
}


/**
 Computes the low-energy states
 */
template<typename HilbertField>
pair<double, string> model_instance<HilbertField>::low_energy_states()
{
  bool is_complex = (typeid(HilbertField) == typeid(Complex));

  // if(!is_correlated or gf_read){ // doing this causes a bug if U=0, in move_submatrix...
  if(gf_read){
    averages.reserve(value.size());
    for(auto& x : value){
      auto X = cluster_averages(the_model->term.at(x.first));
      averages.push_back(tuple<string,double,double>(x.first, X.first, X.second));
    }
  }
  
  if(gs_solved or gf_read) return {GS_energy, GS_string()};
  if(!is_correlated) return one_body_solve();

  // building a set of trial sectors according to the target sector;
  
  if(global_bool("verb_ED")) cout << "computing low-energy state for cluster instance " << full_name() << endl;
  
  sector_set = target_sectors;
  
  GS_energy = 1.0e12; // some large number
  
  // loop over sectors to find the ground state sector
  
  for(auto& s:sector_set){
    the_model->build_HS_operators(s, is_complex);
    Hamiltonian<HilbertField> H = create_hamiltonian(the_model, value, s);
    if(H.dim == 0) continue;
    
    vector<shared_ptr<state<HilbertField>>> gs = H.states(GS_energy); // finds the low-energy states for this sector and adds them to the list
    for(auto& x : gs) states.insert(x);
  }
  
  compute_weights();
  
  // sector_set is now restricted to the sectors containing low-energy states
  sector_set.clear();
  for(auto &x : states) sector_set.insert(x->sec);
  
  // This is one of the places where one switches on the flag sector::up_down in the non-mixed case.
  bool up_down = false;
  for(auto &x : states) if(x->sec.S) up_down = true;
  if(mixing == HS_mixing::normal and up_down) mixing = HS_mixing::up_down;
  
  // average energy and sectors
  GS_energy = (*states.begin())->energy;
  E0 = 0.0;
  for(auto &x : states){
    E0 += x->weight * x->energy;
    if(x->energy < GS_energy) GS_energy = x->energy;
    ostringstream sout;
    sout << *x;
    if(global_bool("verb_ED")) cout << sout.str() << endl;
  }
  
  // average and variances of operators
  averages.reserve(value.size());
  for(auto& x : value){
    auto X = cluster_averages(the_model->term.at(x.first));
    averages.push_back(tuple<string,double,double>(x.first, X.first, X.second));
  }

  gs_solved = true;
  return {GS_energy, GS_string()};
}




/**
 Computes the Boltzmann weights of the different states once they are all known
 */
template<typename HilbertField>
void model_instance<HilbertField>::compute_weights(){
  
  double temperature = global_double("temperature");
  
  if(temperature < SMALL_VALUE){
    // looping over states to keep only the lowest-energy states
    double E0 = 1e12;
    for(auto &x : states) if(x->energy < E0) E0 = x->energy;
    
    // erasing the excited states
    auto it = states.begin();
    while(it != states.end()){
      if((*it)->energy - E0 > SMALL_VALUE){
        ((shared_ptr<state<HilbertField>>)*it).reset();
        states.erase(it++);
      }
      else ++it;
    }

    // all remnaining states are degenerate
    double Z = 1.0/states.size();
    for(auto &x : states){
      x->weight = Z;
    }
    
    if(states.size()==0) qcm_ED_throw("instance " + full_name() + " has no low-energy state!");
    E0 = (*states.begin())->energy;
    return;
  }
  
  // compute the partition function
  double Z = 0.0;
  for(auto &x : states){
    x->weight = exp(-(x->energy-GS_energy)/temperature);
    Z += x->weight;
  }
  E0 = GS_energy-temperature*log(Z);
  
  Z = 1.0/Z;
  for(auto &x : states){
    x->weight *= Z;
  }
  
}




/**
 Calculates the matrices \a tc, \a tb and \a tcb
 */
template<typename HilbertField>
void model_instance<HilbertField>::set_hopping_matrix(bool spin_down)
{
  if(hopping_solved) return;
  matrix<HilbertField>& my_tc = (spin_down and mixing&HS_mixing::up_down)? tc_down :tc;
  matrix<HilbertField>& my_tcb = (spin_down and mixing&HS_mixing::up_down)? tcb_down :tcb;
  matrix<HilbertField>& my_tb = (spin_down and mixing&HS_mixing::up_down)? tb_down :tb;
  
  my_tc.zero();
  for(auto& x : value){
    Hermitian_operator& op = *the_model->term[x.first];
    switch(op.target){
      case 1:
        op.set_hopping_matrix(value.at(x.first), my_tc, spin_down, mixing);
        break;
      case 2:
        op.set_hopping_matrix(value.at(x.first), my_tb, spin_down, mixing);
        break;
      case 3:
        op.set_hopping_matrix(value.at(x.first), my_tcb, spin_down, mixing);
        break;
    }
  }

	if(spin_down == false) SEF_bath = 0.0;
	
  // diagonalize the bath
  if(the_model->n_bath){
    if(spin_down) {
      tb_nd_down = tb_down;
      tcb_nd_down = tcb_down;
    }
    else{
      tb_nd = tb;
      tcb_nd = tcb;
    }
    if(my_tb.is_diagonal() == 0){
      matrix<HilbertField> tcb_tmp(my_tcb);
      vector<double> d(my_tb.r);
      matrix<HilbertField> Sb(my_tb);
      my_tb.eigensystem(d, Sb);
      my_tb.zero();
      my_tb.diagonal(d);
      my_tcb.product(tcb_tmp,Sb);
    }
		
		// contribution of the  bath to the Potthoff functional
		for(size_t i=0; i<the_model->n_bath; ++i) if(realpart(my_tb(i,i)) < 0) SEF_bath -= realpart(my_tb(i,i));
		if(mixing&HS_mixing::up_down and spin_down){
			for(size_t i=0; i<the_model->n_bath; ++i) if(realpart(my_tb(i,i)) < 0) SEF_bath -= realpart(my_tb(i,i));
		}
		else if(mixing==HS_mixing::normal) SEF_bath *= 2.0;
	}
  if((mixing & HS_mixing::up_down) and spin_down == false) set_hopping_matrix(true);
  hopping_solved = true;
}



/**
 inserts a state in the list of low-energy states
 */
template<typename HilbertField>
void model_instance<HilbertField>::insert_state(shared_ptr<state<HilbertField>> S){
  
  // action depends on energy of states compared to previously minimum energy in the set
  if(S->energy < GS_energy){
    GS_energy = S->energy;
    GS_sector = S->sec;
    // removing states that have too high an energy
    auto it = states.begin();
    while(it != states.end()){
      if((*it)->energy - GS_energy > max_gap){
        ((shared_ptr<state<HilbertField>>)*it).reset();
        states.erase(it++);
      }
      else ++it;
    }
    states.insert(S);
  }
  else if(S->energy - GS_energy > max_gap) {
    return;
  }
  else{
    states.insert(S);
  }
}






/**
 computes the representation of the Green function
 */
template<typename HilbertField>
void model_instance<HilbertField>::Green_function_solve()
{
  if(gf_solved or gf_read) return;
  if(!is_correlated){
    one_body_solve();
  	gf_solved = true;
    return;
  }
	if(!gs_solved) low_energy_states();
  
  if(global_bool("verb_ED")) cout << "GF solver for cluster instance " << full_name() << endl;

  for(auto& x : states){

    bool is_complex = (typeid(HilbertField) == typeid(Complex));

    if(GF_solver == GF_format_CF) build_cf(*x, false);
    else build_qmatrix(*x,false);
    
    if(mixing&HS_mixing::up_down){
      if(GF_solver == GF_format_CF) build_cf(*x,true);
      else build_qmatrix(*x,true);
    }
  }
	gf_solved = true;
}






/**
 Evaluates the Green function matrix g (column-order format) at complex frequency z
 */
template<typename HilbertField>
matrix<complex<double>> model_instance<HilbertField>::Green_function(const Complex &z, bool spin_down, bool blocks)
{
  #pragma omp master
  {
    if(spin_down and !(mixing&HS_mixing::up_down or mixing==0))
      qcm_ED_throw("spin_down=True impossible with Hilbert space mixing "+to_string(mixing));
    if(!gf_solved) Green_function_solve();
  }
  block_matrix<Complex> gf_block_matrix(the_model->group->site_irrep_dim*n_mixed);
  if(spin_down and mixing&HS_mixing::up_down){
    for(auto& x : states) x->gf_down->Green_function(z, gf_block_matrix);
  }
  else{
    for(auto& x : states) x->gf->Green_function(z, gf_block_matrix);
  }

  matrix<Complex> G(gf_block_matrix.r);

  if((mixing & HS_mixing::anomalous) && global_bool("strip_anomalous_self")){
    gf_block_matrix.inverse();
    the_model->group->to_site_basis(gf_block_matrix, G, n_mixed);
    auto H = hybridization_function(z, spin_down).v;
    if(the_model->n_bath) G.v += H;
    int R = G.r/2;
    for(int r = R; r<G.r; r++){
      for(int c = 0; c<R; c++){
        G(r,c) = 0.0;
        G(c,r) = 0.0;
      }
    }
    if(the_model->n_bath) G.v -= H;
    G.inverse();
  }
  else if(blocks){
    G = gf_block_matrix.build_matrix();
  }
  else{
    the_model->group->to_site_basis(gf_block_matrix, G, n_mixed);
  }

  return G;
}



/**
 Evaluates the Green function average G (column-order format)
 */
template<typename HilbertField>
matrix<Complex> model_instance<HilbertField>::Green_function_average(bool spin_down)
{
  if(spin_down and !(mixing&HS_mixing::up_down)) qcm_ED_throw("spin_down=True impossible with Hilbert space mixing "+to_string(mixing));
  if(!gf_solved and !gf_read) Green_function_solve();

  matrix<Complex> Gint(n_mixed*(the_model->n_sites));
  block_matrix<Complex> G(the_model->group->site_irrep_dim*n_mixed);
  if(spin_down and mixing&HS_mixing::up_down){
    for(auto& x : states) x->gf_down->integrated_Green_function(G);
  }
  else{
    for(auto& x : states) x->gf->integrated_Green_function(G);
  }
  the_model->group->to_site_basis(G, Gint, n_mixed);
  return Gint;
}





/**
 Evaluates the self energy matrix g (column-order format) at complex frequency z
 */
template<typename HilbertField>
matrix<Complex> model_instance<HilbertField>::self_energy(const Complex &z, bool spin_down)
{
  block_matrix<Complex> gf_block_matrix(the_model->group->site_irrep_dim*n_mixed);
  if(spin_down and mixing&HS_mixing::up_down){
    for(auto& x : states) x->gf_down->Green_function(z, gf_block_matrix);
  }
  else{
    for(auto& x : states) x->gf->Green_function(z, gf_block_matrix);
  }
  gf_block_matrix.inverse();
  matrix<Complex> S(gf_block_matrix.r);
  the_model->group->to_site_basis(gf_block_matrix, S, n_mixed);
  S += (spin_down and mixing&HS_mixing::up_down) ? tc_down : tc;
  if(the_model->n_bath) S.v += hybridization_function(z, spin_down).v;
  S.v *= -1;
  S += z;
  return S;
}



/**
 Constructs the Q_matrix (Lehmann representation) from the Band Lanczos method,
 or full diagonalization if the dimension is small enough.
 @param Omega		state on which the GF is built
 @param spin_down			true if we are building the spin down component of the GF
 */
template<typename HilbertField>
void model_instance<HilbertField>::build_qmatrix(state<HilbertField> &Omega, bool spin_down)
{
  auto& sym_orb = the_model->sym_orb[mixing];
  Q_matrix_set<HilbertField> Qm(the_model->group, mixing);
  Q_matrix_set<HilbertField> Qp(the_model->group, mixing);

  if(global_bool("verb_ED")){
    cout << "\ncomputing Q-matrix for state of energy " << Omega.energy << " in sector " << Omega.sec.name() << endl;
    if(spin_down) cout << "spin down part" << endl;
  }
 
  // building the Q matrices
  int ns = 2*sym_orb.size();
  // #pragma omp parallel for
  for(int s=0; s< ns; s++){
    int r = s/2;
    int pm = 2*(s%2)-1;

    if(sym_orb[r].size() == 0) continue; // irrep not present
    int spin = (spin_down)? 1:sym_orb[r][0].spin;
    sector target_sec = the_model->group->shift_sector(Omega.sec, pm, spin, r);
    if(!the_model->group->sector_is_valid(target_sec)) continue; // target sector is null
    vector<vector<HilbertField>> phi(sym_orb[r].size());
    bool skip_sector=false;
    for(size_t i=0; i< sym_orb[r].size(); i++){
      symmetric_orbital sorb = sym_orb[r][i];
      if(spin_down) sorb.spin =1;
      if(!the_model->create_or_destroy(pm, sorb, Omega, phi[i], HilbertField(1.0))) skip_sector = true;
    }
    if(skip_sector) continue;
    // Assembling the Hamiltonian and Band Lanczos procedure
    Hamiltonian<HilbertField> H = create_hamiltonian(the_model, value, target_sec);
    if(H.dim==0) continue;
    
    Q_matrix<HilbertField> Qtmp;
    Qtmp = H.build_Q_matrix(phi);
    
    Qtmp.e -= Omega.energy; // adjust the eigenvalues by adding/subtracting the GS energy
    if(pm == -1){
      Qtmp.e *= -1.0;
    }
    Qtmp.streamline();
    if(pm==-1) 
      Qm.q[r] = Qtmp;
    else 
      Qp.q[r] = Qtmp;
      Qp.q[r].v.cconjugate(); // IMPORTANT. Source of bug found 2021-08-14
  }

  auto Q = make_shared<Q_matrix_set<HilbertField>>(the_model->group, mixing);
  if(spin_down) Omega.gf_down = Q;
  else Omega.gf = Q;
  for(size_t r=0; r<sym_orb.size(); r++){
    Q->q[r].append(Qm.q[r]);
    Q->q[r].append(Qp.q[r]);
    Q->q[r].check_norm(global_double("accur_Q_matrix"));
    Q->q[r].v.v *= sqrt(Omega.weight);
  }
}







/**
 Computes the hybridization matrix Gamma
 Gamma is \f$-G'^{-1}-\Sigma-t'+\omega\f$.
 It is a projection onto the cluster sites only (bath excluded).
 */
template<typename HilbertField>
matrix<Complex> model_instance<HilbertField>::hybridization_function(Complex w, bool spin_down){
  set_hopping_matrix(false);
  size_t n_bathg = tb.r;
  matrix<HilbertField>& my_tb = (spin_down and mixing&HS_mixing::up_down)? tb_down : tb;
  matrix<HilbertField>& my_tcb = (spin_down and mixing&HS_mixing::up_down)? tcb_down : tcb;
  
  matrix<Complex> Gamma(dim_GF);
  for(size_t i=0; i<n_bathg; ++i){
    Complex z = 1.0/(w - my_tb(i,i));
    for(size_t a=0; a < Gamma.r; ++a){
      for(size_t b=0; b < Gamma.r; ++b) Gamma(a,b) +=  my_tcb(a,i) * conjugate(my_tcb(b,i)) * z;
    }
  }
  return Gamma;
}



/**
 Computes the susceptibility of the operator h
w : array of complex frequencies
 */
template<typename HilbertField>
vector<Complex> model_instance<HilbertField>::susceptibility(shared_ptr<Hermitian_operator> h, const vector<Complex> &w){
  
  if(value.find(h->name)==value.end()) qcm_ED_throw("operator "+h->name+" is not activated in instance "+to_string(label));
  if(!gs_solved) low_energy_states();
  
  vector<Complex> chi(w.size(),0.0);
  
  for(auto& sec : sector_set){
    Hamiltonian<HilbertField> H = create_hamiltonian(the_model, value, sec);
    if(H.dim==0) continue;
    for(auto& gs : states){
      if(gs->sec != sec) continue;
      vector<vector<HilbertField>> b(1);
      b[0].resize(H.dim);
      h->HS_operator.at(sec)->multiply_add(gs->psi, b[0], 1.0);
      Q_matrix<HilbertField> Q = H.build_Q_matrix(b);
      Q.e -= gs->energy; 	// adjust the eigenvalues by subtracting the GS energy
      matrix<Complex> g(1);
      for(size_t j=0; j<w.size(); j++){
        g.zero();
        Q.Green_function(Complex(-w[j].real(),w[j].imag()),g);
        g.v *= -1.0;
        Q.Green_function(w[j],g);
        chi[j] += gs->weight*g(0,0);
      }
    }
  }
  return chi;
}



/**
 Computes the susceptibility of the operator h in the form of a Lehmann representation
 returns an array of poles and residues
 */
template<typename HilbertField>
vector<pair<double,double>> model_instance<HilbertField>::susceptibility_poles(shared_ptr<Hermitian_operator> h){
  
  if(value.find(h->name)==value.end()) qcm_ED_throw("operator "+h->name+" is not activated in instance "+to_string(label));
  if(!gs_solved) low_energy_states();
  
  vector<pair<double,double>> chi;
  chi.reserve(20);
  
  for(auto& sec : sector_set){
    Hamiltonian<HilbertField> H = create_hamiltonian(the_model, value, sec);
    if(H.dim==0) continue;
    for(auto& gs : states){
      if(gs->sec != sec) continue;
      vector<vector<HilbertField>> b(1);
      b[0].resize(H.dim);
      h->HS_operator.at(sec)->multiply_add(gs->psi, b[0], 1.0);
      Q_matrix<HilbertField> Q = H.build_Q_matrix(b);
      Q.e -= gs->energy; 	// adjust the eigenvalues by subtracting the GS energy
      for(size_t i = 0; i<Q.M; i++){
        if(Q.e[i] < 1e-8) continue;
        double r = abs(Q.v(0,i));
        r *= r*gs->weight;
        if(r > 1e-6) chi.push_back({Q.e[i],r});
      }
    }
  }
  // consolidation
  vector<pair<double,double>> chi2;
  chi2.reserve(20);
  for(size_t i = 0; i< chi.size(); i++){
    size_t j;
    for(j = 0; j< chi2.size(); j++){
      if(fabs(chi2[j].first - chi[i].first) < 1e-4){
        chi2[j].second += chi[i].second; break;
      }
    }
    if(j==chi2.size()) chi2.push_back(chi[i]);
  }
  return chi2;
}




/**
 prints the instance details on a file
 */
template<typename HilbertField>
void model_instance<HilbertField>::print(ostream& fout)
{
  fout << "\n\nmodel instance no " << label << " based on " << the_model->name;
  if(complex_Hilbert) fout << "  (complex)";
  fout << endl;
  fout << "Symmetric orbitals:\n";
  size_t i = 1;
  for(auto& x : the_model->sym_orb[mixing]){
    for(auto& y : x) {fout << i << '\t' << y << endl; i++;}
  }
  fout << endl;
  for(auto& x: value) fout << x.first << " :\t" << x.second << endl;

  set_hopping_matrix(false);
  fout << "\nhopping matrix:\n" << tc << endl;
  if(the_model->n_bath){
    fout << "\nbath hopping matrix:\n" << tb << endl;
    fout << "\nbath hybridization matrix:\n" << tcb << endl;
    fout << "\nbath hopping matrix (non diagonalized):\n" << tb_nd << endl;
    fout << "\nbath hybridization matrix (non diagonalized):\n" << tcb_nd << endl;
  }
  if(mixing&HS_mixing::up_down){
    fout << "\nhopping matrix (spin down):\n" << tc_down << endl;
    if(the_model->n_bath){
      fout << "\nbath hopping matrix (spin down):\n" << tb_down << endl;
      fout << "\nbath hybridization matrix (spin down):\n" << tcb_down << endl;
      fout << "\nbath hopping matrix (spin_down, non diagonalized):\n" << tb_nd_down << endl;
      fout << "\nbath hybridization matrix (spin_down, non diagonalized):\n" << tcb_nd_down << endl;
    }
  }
}




/**
 Constructs a ground state and Green function for the one-body problem only
 Useful when the number of sites exceeds the limit of ED, for studying impurity problems
 */
template<typename HilbertField>
pair<double, string> model_instance<HilbertField>::one_body_solve()
{
  if(the_model->group->g != 1) qcm_ED_throw("The symmetry group must be trivial when using 'one_body_solve()'");

  GS_energy = 0.0;
  set_hopping_matrix(false);
  GF_solver = GF_format_BL;
  auto S = make_shared<state<HilbertField>>();
  auto Qset = make_shared<Q_matrix_set<HilbertField>>(the_model->group, mixing);
  S->gf =  Qset;
  int d = tc.r+tb.r;
  Q_matrix<HilbertField> Q(tc.r, d);
  matrix<HilbertField> H(d);
  tc.move_sub_matrix(tc.r,tc.r, 0, 0, 0, 0, H);
  if(tb.r>0){
    tb.move_sub_matrix(tb.r,tb.r, 0, 0, tc.r, tc.r, H);
    tcb.move_sub_matrix(tcb.r,tcb.c, 0, 0, 0, tc.r, H);
    tcb.move_sub_matrix_HC(tcb.r,tcb.c, 0, 0, tc.r, 0, H);
  }

  matrix<HilbertField> U(d);
  H.eigensystem(Q.e, U);
  U.move_sub_matrix(tc.r, d, 0, 0, 0, 0, Q.v);
  Q.v.cconjugate();
  Qset->q[0] = Q;
  for(auto& x:Q.e) if(x < 0.0) GS_energy += x;


  // computing M, the average M_{ab} = \L c^\dagger_b c_a \R
  M.set_size(d);
  for(int a=0; a<d; a++){
    for(int b=0; b<d; b++){
      for(int r=0; r<Q.e.size(); r++){
        if(Q.e[r] < 0.0) M(b,a) += U(a,r)*conjugate(U(b,r));
      } 
    }
  }
  Q.check_norm(1e-6);

  if(mixing&HS_mixing::up_down){
    auto Qset_down = make_shared<Q_matrix_set<HilbertField>>(the_model->group, mixing);
    S->gf_down =  Qset_down;
    Q_matrix<HilbertField> Q_down(tc.r, d);
    matrix<HilbertField> H_down(d);
    tc_down.move_sub_matrix(tc.r,tc.r, 0, 0, 0, 0, H_down);
    if(tb.r>0){
      tb_down.move_sub_matrix(tb.r,tb.r, 0, 0, tc.r, tc.r, H_down);
      tcb_down.move_sub_matrix(tcb.r,tcb.c, 0, 0, 0, tc.r, H_down);
      tcb_down.move_sub_matrix_HC(tcb.r,tcb.c, 0, 0, tc.r, 0, H_down);
    }
    matrix<HilbertField> U(d);
    H_down.eigensystem(Q_down.e, U);
    U.move_sub_matrix(tc.r, d, 0, 0, 0, 0, Q_down.v);
    Q_down.v.cconjugate();
    Qset_down->q[0] = Q_down;
    for(auto& x:Q.e) if(x < 0.0) GS_energy += x;

    // computing M_down, the average M_{ab} = \L c^\dagger_b c_a \R
    M_down.set_size(d);
    for(int a=0; a<d; a++){
      for(int b=0; b<d; b++){
        for(int r=0; r<Q.e.size(); r++){
          if(Q_down.e[r] < 0.0) M_down(b,a) += U(a,r)*conjugate(U(b,r));
        } 
      }
    }

  }

/*
  // Nambu correction to the ground state energy
  double nambu_corr = 0.0;
  if(mixing == HS_mixing::anomalous){
    for(auto& x : value) nambu_corr += the_model->term.at(x.first)->nambu_correction * x.second;
  }
  else if(mixing == HS_mixing::full){
    for(auto& x : value) nambu_corr += the_model->term.at(x.first)->nambu_correction_full * x.second;
  }


  GS_energy += nambu_corr;
  if(mixing == HS_mixing::normal) GS_energy *= 2.0;
  else if(mixing == HS_mixing::full) GS_energy *= 0.5;
  S->energy = GS_energy;
  S->weight = 1.0;
  */
  states.insert(S);
  gf_solved = true;
  gs_solved = true;

  return {GS_energy, "uncorrelated"};
}




template<typename HilbertField>
double model_instance<HilbertField>::fidelity(model_instance<HilbertField>& inst)
{
  low_energy_states();
  inst.low_energy_states();
  double z(0);
  for(auto& S1 : states){
    for(auto& S2 : inst.states){
      if(S1->sec != S2->sec) continue;
      HilbertField zz = S1->psi * S2->psi;
      z += abs(zz)*abs(zz) * S1->weight * S2->weight;
    }
  }
  return z;
}


#define NW 3
template<typename HilbertField>
double model_instance<HilbertField>::tr_sigma_inf()
{
  int d1 = dim_GF;
  int d2 = d1;
  if(mixing&1) d2 /= 2;
  
  vector<double> iw(NW);
  vector<double> S(NW);
  
  for(int i=0; i<iw.size(); i++) iw[i] = (i+1)*1.0e-5;
  for(int i=0; i<iw.size(); i++){
    complex<double> z(0.0);
    auto sigma = self_energy(complex<double>(0.0, 1.0/iw[i]), false);
    for(int j=0; j<d2; j++) z += sigma(j,j);
    if(mixing&1) for(int j=d2; j<d1; j++) z -= sigma(j,j);
    S[i] = real<double>(z);
  }
  if(mixing&HS_mixing::up_down){
    for(int i=0; i<iw.size(); i++){
      complex<double> z(0.0);
      auto sigma = self_energy(complex<double>(0.0, 1.0/iw[i]), true);
      for(int j=0; j<d2; j++) z += sigma(j,j);
      S[i] += real<double>(z);
    }
  }
  if(mixing == HS_mixing::normal){
    for(int i=0; i<iw.size(); i++) S[i] *= 2;
  }
  
  // extrapolation to infinite frequency
  double Sinf, Sinf_err;
  polynomial_fit(iw, S, 0.0, Sinf, Sinf_err);
  return Sinf;
}




template<typename HilbertField>
string model_instance<HilbertField>::GS_string() const
{
  ostringstream sout;
  map<sector, double> weight;
  for(auto& s : states){
    weight[s->sec] += s->weight;
  }
  for(auto& s : weight){
    if(s.second < 0.001) continue;
    sout << s.first << ':' << setprecision(3) << s.second << '/';
  }
  string out = sout.str();
  out.pop_back();
  return out;
}


/**
 writes the Green function information to a stream.
 This can be read instead of solving the model
 */
template<typename HilbertField>
void model_instance<HilbertField>::write(ostream& fout)
{
  if(!gf_solved) Green_function_solve();
  
  // writing the info line
  // fout << "cluster: " << the_model->name << '\n';
  for (auto &x : value) fout << x.first << '\t' << x.second << '\n';
  fout << "\nGS_energy: " << GS_energy << " GS_sector: " << GS_string() << '\n';
  if(GF_solver == GF_format_CF) fout << "GF_format: cf\n";
  else fout << "GF_format: bl\n";
  fout << "mixing\t" << mixing << endl;
  
  for(auto& x: states){
    x->write(fout);
  }
}




template<typename HilbertField>
void model_instance<HilbertField>::read(istream& fin)
{
  while(true){
    vector<string> input = read_strings(fin);
    if(input.size()==0) break;
    if(input.size()!=2) qcm_ED_throw("failed to read a parameter in input. Need two columns per parameter");
    // if(value.find(input[0])==value.end()) qcm_ED_throw("unkown parameter "+input[0]+" in solutions file");
    // value[input[0]] = from_string<double>(input[1]);
  }
  
  string tmp;
  fin >> "GS_energy:" >> GS_energy >> tmp >> tmp;
  fin >> "GF_format:" >> tmp;
  if(tmp == "bl") GF_solver = GF_format_BL;
  else GF_solver = GF_format_CF;
  fin >> "mixing" >> mixing;
  n_mixed = 1;
  if(mixing & HS_mixing::anomalous) n_mixed *= 2;
  if(mixing & HS_mixing::spin_flip) n_mixed *= 2;


  // reading the states
  states.clear();
  while(true){
    parser::next_line(fin);
    vector<string> input = read_strings(fin);
    if(input.size()==0) break;
    states.insert(shared_ptr<state<HilbertField>>(new state<HilbertField>(fin, the_model->group, mixing, GF_solver)));
  }
  gf_read = true;
  M = Green_function_average(false);
  if(mixing & HS_mixing::up_down) M_down = Green_function_average(true);
}



template<typename HilbertField>
pair<vector<double>, vector<complex<double>>> model_instance<HilbertField>::qmatrix(bool spin_down)
{
  if(GF_solver != GF_format_BL) qcm_ED_throw("Green function format is not Lehmann! Cannot output the Q matrix.");
  if(!gf_solved) Green_function_solve();
  if(states.size() > 1)  qcm_ED_throw("The ground state is not a pure state! Cannot output the Q matrix.");
  shared_ptr<Green_function_set> gf;
  if(spin_down) gf = (*states.begin())->gf_down;
  else gf = (*states.begin())->gf;
  auto Q = dynamic_pointer_cast<Q_matrix_set<HilbertField>>(gf)->consolidated_qmatrix();
  return {Q.e, to_complex(Q.v.v)};
}


template<typename HilbertField>
pair<vector<double>, vector<complex<double>>> model_instance<HilbertField>::hybridization(bool spin_down)
{
  if(the_model->n_bath == 0) qcm_ED_throw("System has no bath. Hybridization not defined.");
  if(!hopping_solved) set_hopping_matrix(spin_down);
  vector<double> E(tb.r);
  if(spin_down){
    for(size_t i=0; i<tb.r; i++) E[i] = real<double>(tb_down(i,i));
    return {E, to_complex(tcb_down.v)};
  }
  for(size_t i=0; i<tb.r; i++) E[i] = real<double>(tb(i,i));
  return {E, to_complex(tcb.v)};
}
#endif /* model_instance_hpp */




template<typename HilbertField>
void model_instance<HilbertField>::print_wavefunction(ostream& fout)
{
  for (auto& s : states){
    if(the_model->is_factorized){
      auto B = the_model->provide_factorized_basis(s->sec);
      s->write_wavefunction(fout, *B);
    }
    else{
      auto B = the_model->provide_basis(s->sec);
      s->write_wavefunction(fout, *B);
    }
  }
}




template<typename HilbertField>
matrix<Complex> model_instance<HilbertField>::hopping_matrix_full(bool spin_down)
{
  matrix<Complex> H(tc.r+tb.r);
  if(spin_down){
    tc_down.move_sub_matrix(tc.r, tc.c, 0, 0, 0, 0, H);
    tb_down.move_sub_matrix(tb.r, tb.c, 0, 0, tc.r, tc.c, H);
    tcb_down.move_sub_matrix(tcb.r, tcb.c, 0, 0, 0, tc.c, H);
    tcb_down.move_sub_matrix_HC(tcb.r, tcb.c, 0, 0, tc.r, 0, H);
  }
  else{
    tc.move_sub_matrix(tc.r, tc.c, 0, 0, 0, 0, H);
    tb.move_sub_matrix(tb.r, tb.c, 0, 0, tc.r, tc.c, H);
    tcb.move_sub_matrix(tcb.r, tcb.c, 0, 0, 0, tc.c, H);
    tcb.move_sub_matrix_HC(tcb.r, tcb.c, 0, 0, tc.r, 0, H);
  }
  return H;
}



template<typename HilbertField>
vector<tuple<int,int,double>> model_instance<HilbertField>::interactions()
{
  vector<tuple<int,int,double>> V;
  V.reserve(the_model->n_sites);
  for(auto& x : value){
    if(the_model->term.at(x.first)->is_interaction){
      auto elem = the_model->term.at(x.first)->matrix_elements();
      for(auto& y : elem){
        V.push_back(tuple<int,int,double>(y.r, y.c, y.v.real()*x.second));
      }
    }
  }
  return V;
}




template<typename HilbertField>
pair<matrix<Complex>, vector<uint64_t>> model_instance<HilbertField>::density_matrix_mixed(vector<int> sites)
{
  if(the_model->group->g > 1) qcm_ED_throw("computing the density matrix requires symmetries to be turned off!");

  if(!gs_solved) low_energy_states();
  state<HilbertField> psi = **states.begin();
  shared_ptr<ED_mixed_basis> base = the_model->basis[psi.sec]; 

  // psi.write_wavefunction(cout, *base); // ***

  int nA = sites.size();
  int nB = the_model->n_orb - sites.size();
  vector<int> sitesB(0);
  for(int i=0; i<the_model->n_orb; i++){
    int j;
    for(j=0; j<sites.size(); j++) if(sites[j] == i) break;
    if(j == sites.size()) sitesB.push_back(i);
  }

  // constructing the masks for subsystems A and B
  uint64_t mask_A = 0L;
  uint64_t mask_B = 0L;

  for(int i=0; i<sites.size(); i++)
    mask_A += (1 << sites[i]);
  mask_A = mask_A + (mask_A << 32);
  vector<int> S(sites); S += 1;
  S = sitesB; S+=1;

  // First step: decomposing the basis into the two subsystems
  map<uint64_t, int> index_A;
  map<uint64_t, int> index_B;
  vector<uint64_t> binA;
  vector<uint64_t> binB;
  vector<vector<int>> listA; 
  int i=0;
  int j=0;
  for(int k=0; k<base->dim; k++){
    uint64_t A = base->binlist[k].b & mask_A;
    uint64_t B = base->binlist[k].b ^ A;
    if (index_A.find(A) == index_A.end()){
      index_A[A] = i++;
      // binA.push_back(collapse(A,sites));
      binA.push_back(A);
    }
    if (index_B.find(B) ==  index_B.end()){
      index_B[B] = j++; 
      listA.push_back(vector<int>(0));
      // binB.push_back(collapse(B,sitesB));
      binB.push_back(B);
    }
    listA[index_B[B]].push_back(index_A[A]);
  }
  int dim_A = i;
  int dim_B = j;

  matrix<Complex> rho(dim_A);
  // loop over states
  for(int j=0; j<dim_B; j++){
    vector<int> &v = listA[j];
    // cout << j << " : " << v << endl;
    for(int i=0; i<v.size(); i++){
      auto bI = binary_state(binA[v[i]] | binB[j]);
      int I = base->index(bI);
      for(int ip=0; ip<v.size(); ip++){
        auto bIp = binary_state(binA[v[ip]] | binB[j]);
        int Ip = base->index(bIp);
        rho(v[i], v[ip]) += conj(psi.psi[I])*psi.psi[Ip];
      }
    }
  }

  vector<uint64_t> binA_collapsed(dim_A);
  for(int i=0; i<dim_A; i++) binA_collapsed[i] = collapse(binA[i], sites);

  return {rho, binA_collapsed};
}




template<typename HilbertField>
pair<matrix<Complex>, vector<uint64_t>> model_instance<HilbertField>::density_matrix_factorized(vector<int> sites)
{
    // loop over states
    return {matrix<Complex>(0), vector<uint64_t>(0)};
}
