#include "model.hpp"
#include "Hermitian_operator.hpp"
#include <fstream>
#ifdef _OPENMP
  #include <omp.h>
#endif

/**
 Constructor
 @param _name   name of system (cluster). Used to distinguish different cluster shapes within a larger calculation
 @param _n_orb	number of orbitals (bath included) in the model
 @param _n_bath	number of bath orbitals in the model
 @param gen sequence of generators (permutations of sites)
  */

model::model(const string &_name, const size_t _n_orb, const size_t _n_bath, const vector<vector<int>> &gen, bool bath_irrep)
: name(_name),
n_orb(_n_orb),
n_bath(_n_bath),
n_sites(_n_orb-_n_bath),
is_closed(false),
is_factorized(false)
{
  group = make_shared<symmetry_group>(n_orb, n_sites, gen, bath_irrep);
  in_bath.assign(2*n_orb, false);
  for(size_t i=n_sites; i<n_orb; i++){
    in_bath[i] = true;
    in_bath[i+n_orb] = true;
  }

  // builds symmetric orbitals for all mixings:
  sym_orb.reserve(6);
  for(int m=0; m<6; m++) sym_orb.push_back(group->build_symmetric_orbitals(m));
}







/**
 Prints the definition of the model into a stream
 */
void model::print(ostream& fout)
{
  console::banner('-', "cluster " + name, fout);
  fout << *group << endl;
  for(auto& x: term) x.second->print(fout);
}




/**
 Builds the bases and operators needed for the computation of the ground state
 @param GS_sector sector of the state 
 @param complex true if the HS operators must be complex-valued
*/
void model::build_HS_objects(const sector& sec, bool is_complex)
{
  if(is_factorized){
    if(factorized_basis.find(sec) == factorized_basis.end()) factorized_basis[sec] = make_shared<ED_factorized_basis>(sec, n_orb);
  }
  else{
    if(basis.find(sec) == basis.end()) basis[sec] = make_shared<ED_mixed_basis>(sec, group);
  }
  if(Hamiltonian_format == H_FORMAT::H_format_onthefly) return;

  // building the operators
  for(auto& x : term){
    if(!x.second->is_active) continue;
    if(x.second->HS_operator.find(sec) == x.second->HS_operator.end())
      x.second->HS_operator[sec] = x.second->build_HS_operator(sec, is_complex);
  }
}



/**
 Builds the bases and operators needed for the computation of the Green function
 @param GS_sector sector of the state on which the Green function is to be built
 @param spin_down true if we are building the spin down component of the GF
 @param complex true if the HS operators must be complex-valued
 */
void model::build_HS_objects_GF(const sector& GS_sector, int mixing, bool spin_down, bool is_complex)
{
  auto needed_secs = needed_sectors_GF(GS_sector, mixing, spin_down);
  int n_sec = needed_secs.size();

  // building the bases
  vector<shared_ptr<ED_mixed_basis>> needed_bases(n_sec);
  vector<shared_ptr<ED_factorized_basis>> needed_factorized_bases(n_sec);

  if(is_factorized){
    #pragma omp parallel for
    for(int i=0; i<n_sec; i++){
      if(basis.find(needed_secs[i].first) == basis.end())
        needed_factorized_bases[i] = make_shared<ED_factorized_basis>(needed_secs[i].first, n_orb);
    }

    for(int i=0; i<n_sec; i++){
      if(factorized_basis.find(needed_secs[i].first) == factorized_basis.end())
        factorized_basis[needed_secs[i].first] = needed_factorized_bases[i];
      // else needed_factorized_bases[i] = factorized_basis[needed_secs[i].first];
    }
  }
  else{
    #pragma omp parallel for
    for(int i=0; i<n_sec; i++){
      if(basis.find(needed_secs[i].first) == basis.end())
        needed_bases[i] = make_shared<ED_mixed_basis>(needed_secs[i].first, group);
    }

    for(int i=0; i<n_sec; i++){
      if(basis.find(needed_secs[i].first) == basis.end())
        basis[needed_secs[i].first] = needed_bases[i];
      // else needed_bases[i] = basis[needed_secs[i].first];
    }
  }

  if(Hamiltonian_format == H_FORMAT::H_format_onthefly) return;

  if(!is_factorized){
    // building the creation and annihilation operators
    // list of needed destruction identifiers
    vector<destruction_identifier> D;
    vector<shared_ptr<ED_mixed_basis>> needed_B_bases;
    vector<shared_ptr<ED_mixed_basis>> needed_T_bases;
    D.reserve(16);
    needed_B_bases.reserve(16);
    needed_T_bases.reserve(16);

    int ns = 2*sym_orb[mixing].size();
    for(int s=0; s< ns; s++){
      int r = s/2;
      int pm = 2*(s%2)-1;
      if(sym_orb[mixing][r].size() == 0) continue; // irrep not present
      int spin = (spin_down)? 1:sym_orb[mixing][r][0].spin;
      sector target_sec = group->shift_sector(GS_sector, pm, spin, r);
      if(!group->sector_is_valid(target_sec)) continue; // target sector is null
      for(size_t i=0; i< sym_orb[mixing][r].size(); i++){
        symmetric_orbital sorb = sym_orb[mixing][r][i];
        if(spin_down) sorb.spin =1;
        if((pm == 1 and sorb.nambu==0) or (pm == -1 and sorb.nambu==1)){
          D.push_back({target_sec, sorb});
          needed_B_bases.push_back(basis.at(target_sec));
          needed_T_bases.push_back(basis.at(GS_sector));
        }  
        else{
          D.push_back({GS_sector, sorb});
          needed_T_bases.push_back(basis.at(target_sec));
          needed_B_bases.push_back(basis.at(GS_sector));
        }
      }
    }

    // cout << "needed destructions operators:\n";
    // for(auto& x : D) cout << x.secB << ":" << x.sorb << ", "; cout << endl;

    if(is_complex){
      vector<shared_ptr<destruction_operator<Complex>>> needed_destruction(D.size());
      #pragma omp parallel for
      for(int i=0; i<D.size(); i++){
        if(destruction_complex.find(D[i])==destruction_complex.end())
          needed_destruction[i] = make_shared<destruction_operator<Complex>>(needed_B_bases[i], needed_T_bases[i], D[i].sorb);
      }
      for(int i=0; i<D.size(); i++){
        if(destruction_complex.find(D[i])==destruction_complex.end()) destruction_complex[D[i]] = needed_destruction[i];
      }
    }
    else{
      vector<shared_ptr<destruction_operator<double>>> needed_destruction(D.size());
      #pragma omp parallel for
      for(int i=0; i<D.size(); i++){
        if(destruction.find(D[i])==destruction.end())
          needed_destruction[i] = make_shared<destruction_operator<double>>(needed_B_bases[i], needed_T_bases[i], D[i].sorb);
      }
      for(int i=0; i<D.size(); i++){
        if(destruction.find(D[i])==destruction.end()) destruction[D[i]] = needed_destruction[i];
      }
    }
  }
  // building the operators
  for(auto& x : term){
    Hermitian_operator& op = *x.second;
    if(!op.is_active) continue;
    vector<shared_ptr<HS_Hermitian_operator>> needed_ops(n_sec);
    #pragma omp parallel for
    for(int i=0; i<n_sec; i++){
      if(op.HS_operator.find(needed_secs[i].first) == op.HS_operator.end())
        needed_ops[i] = op.build_HS_operator(needed_secs[i].first, is_complex or group->complex_irrep[needed_secs[i].first.irrep]);
    }
    for(int i=0; i<n_sec; i++){
      if(op.HS_operator.find(needed_secs[i].first) == op.HS_operator.end())
        op.HS_operator[needed_secs[i].first] = needed_ops[i];
    }
  }
}




/**
 Builds the list of needed sectors for computing the Green function
 @param GS_sector sector of the state on which the Green function is to be built
 @param spin_down true if we are building the spin down component of the GF
 */
vector<pair<sector,int>> model::needed_sectors_GF(const sector& GS_sector, int mixing, bool spin_down)
{
  vector<pair<sector,int>> needed_secs;

  // building the list of needed sectors
  needed_secs.reserve(16);

  for(size_t r=0; r<sym_orb[mixing].size(); r++){
    for(int pm = -1; pm<2; pm += 2){ // loop over destruction (pm = -1) and creation (pm = +1)
      
      // constructing the target sector
      if(sym_orb[mixing][r].size() == 0) continue; // irrep not present
      int spin = (spin_down)? 1:sym_orb[mixing][r][0].spin;
      sector target_sec = group->shift_sector(GS_sector, pm, spin, r);
      if(!group->sector_is_valid(target_sec)) continue;
      needed_secs.push_back({target_sec, pm});
    }
  }
  return needed_secs;
}


/**
Prints a graph representation of the model, using the dot language.
The fixed positions of the cluster sites per se (not the bath) are provided in 'pos'
*/
void model::print_graph(const vector<vector<double>> &pos){
  ofstream fout(name+".dot");

  fout << "graph {\nK=1.3;\n";
  for(int i=0; i<pos.size(); i++){
    fout << i+1 << " [shape=square color=\"blue\" pos=\"" << pos[i][0] << "," << pos[i][1] << "!\"]\n";
  }
  for(int i=pos.size(); i<n_orb; i++){
    fout << i+1 << " [shape=circle color=\"red\"]\n";
  }
  for(auto& x : term){
    if(x.second->is_interaction) continue;
    auto elem = x.second->matrix_elements();
    for(auto& e : elem){
      if(e.r >= n_orb or e.c >= n_orb) continue;
      if(e.r >= e.c) continue;
      string lab = x.second->name;
      if(abs(e.v - 1.0)<1e-4) lab = lab;
      else if(abs(e.v + 1.0)<1e-4) lab = '-'+lab;
      else lab = to_string(e.v)+lab;
      if(e.c >= n_sites) fout << e.r+1 << " -- " << e.c+1 << " [color = \"green\" headlabel=\"" << lab << "\" labeldistance=3 labelangle=0 fontsize=10 fontcolor=\"#990000\"];\n";
      else fout << e.r+1 << " -- " << e.c+1 << " [label=\"" << lab << "\" fontsize=10];\n";
    }
  }
  fout << "}\n";
  fout.close();
}



/**
 Applies a creation or annihilation operator
 |y> += z * c_a |x>  ou  |y> += z * c^+_a |x>
 @param pm		if = 1: creation; if = -1: destruction.
 @param a		symmetric orbital
 @param x		in state
 @param y		out state (allocated before)
 @param z		coefficient of the out state
 */
bool model::create_or_destroy(int pm, const symmetric_orbital &a, state<double> &x, vector<double> &y, double z)
{
  if(a.nambu) pm = -pm;
  sector target_sec = group->shift_sector(x.sec, pm, a.spin, a.irrep);

  if(is_factorized){
    auto B = factorized_basis.at(x.sec);
    auto T = factorized_basis.at(target_sec);
    if(y.size() != T->dim) y.resize(T->dim);
    if(pm==1){ // creation
      for(uint32_t I=0; I<T->dim; ++I){
        auto P = Destroy(a.orb+n_orb*a.spin, I, *T, *B);
        if(!get<3>(P)) continue;
        get<1>(P) = 1-2*(get<1>(P)%2);
        y[I] += z*get<1>(P)*x.psi[get<0>(P)];
      }
    }
    else{ // destruction
      for(uint32_t I=0; I<B->dim; ++I){
        auto P = Destroy(a.orb+n_orb*a.spin, I, *B, *T);
        if(!get<3>(P)) continue;
        get<1>(P) = 1-2*(get<1>(P)%2);
        y[get<0>(P)] += z*get<1>(P)*x.psi[I];
      }
    }
  }
  else{
    auto B = basis.at(x.sec);
    sector target_sec = group->shift_sector(x.sec, pm, a.spin, a.irrep);
    auto T = basis.at(target_sec);
    if(y.size() != T->dim) y.resize(T->dim);
    int n = group->N;

    if(pm==1){ // creation
      if(Hamiltonian_format == H_FORMAT::H_format_onthefly){
        for(uint32_t I=0; I<T->dim; ++I){
          auto P = Destroy(a.orb+n_orb*a.spin, I, *T, *B);
          if(!get<3>(P)) continue;
          get<1>(P) = get<1>(P)%(2*group->g);
          y[I] += z*group->phaseX<double>(get<1>(P))*x.psi[get<0>(P)];
        }
      }
      else{
        destruction_identifier D(target_sec,a);
        if(destruction.find(D)==destruction.end()){
          cout << "destruction operator [" << target_sec << ", " << a << "] not found!" << endl;
        }
        destruction.at(D)->multiply_add_conjugate(x.psi,y,z);
      }
      // cout << "creation at " << a.label << '\n' << x.psi << '\n' << y << "\n\n"; // tempo
    }
    else{ // destruction
      if(Hamiltonian_format == H_FORMAT::H_format_onthefly){
        for(uint32_t I=0; I<B->dim; ++I){
          auto P = Destroy(a.orb+n_orb*a.spin, I, *B, *T);
          if(!get<3>(P)) continue;
          get<1>(P) = get<1>(P)%(2*group->g);
          y[get<0>(P)] += z*group->phaseX<double>(get<1>(P))*x.psi[I];
        }
      }
      else{
        destruction_identifier D(x.sec,a);
        if(destruction.find(D)==destruction.end()){
          cout << "destruction operator [" << x.sec << ", " << a << "] not found!" << endl;
        }
        destruction.at(D)->multiply_add(x.psi,y,z);
      }
      // cout << "destruction at " << a.label << '\n' << x.psi << '\n' << y << "\n\n"; // tempo
    }
  }
  // cout << x.psi << '\n' << y << "\n\n"; // tempo
  return true;
}






/**
 Applies a creation or annihilation operator
 |y> += z * c_a |x>  ou  |y> += z * c^+_a |x>
 @param pm		if = 1: creation; if = -1: destruction.
 @param a		symmetric orbital
 @param x		in state
 @param y		out state (allocated before)
 @param z		coefficient of the out state
 */
bool model::create_or_destroy(int pm, const symmetric_orbital &a, state<Complex> &x, vector<Complex> &y, Complex z)
{
  if(a.nambu) pm = -pm;
  sector target_sec = group->shift_sector(x.sec, pm, a.spin, a.irrep);

  if(is_factorized){
    auto B = factorized_basis.at(x.sec);
    auto T = factorized_basis.at(target_sec);
    if(y.size() != T->dim) y.resize(T->dim);
    if(pm==1){ // creation
      for(uint32_t I=0; I<T->dim; ++I){
        auto P = Destroy(a.orb+n_orb*a.spin, I, *T, *B);
        if(!get<3>(P)) continue;
        get<1>(P) = get<1>(P)%(2*group->g);
        y[I] += z*group->phaseX<Complex>(get<1>(P))*x.psi[get<0>(P)];
      }
    }
    else{ // destruction
      for(uint32_t I=0; I<B->dim; ++I){
        auto P = Destroy(a.orb+n_orb*a.spin, I, *B, *T);
        if(!get<3>(P)) continue;
        get<1>(P) = get<1>(P)%(2*group->g);
        y[get<0>(P)] += z*group->phaseX<Complex>(get<1>(P))*x.psi[I];
      }
    }
  }
  else{
    auto B = basis.at(x.sec);
    sector target_sec = group->shift_sector(x.sec, pm, a.spin, a.irrep);
    auto T = basis.at(target_sec);
    if(y.size() != T->dim) y.resize(T->dim);
    int n = group->N;

    if(pm==1){ // creation
      if(Hamiltonian_format == H_FORMAT::H_format_onthefly){
        for(uint32_t I=0; I<T->dim; ++I){
          auto P = Destroy(a.orb+n_orb*a.spin, I, *T, *B);
          if(!get<3>(P)) continue;
          get<1>(P) = get<1>(P)%(2*group->g);
          y[I] += z*group->phaseX<Complex>(get<1>(P))*x.psi[get<0>(P)];
        }
      }
      else{
        destruction_identifier D(target_sec,a);
        if(destruction_complex.find(D)==destruction_complex.end()){
          cout << "destruction operator [" << target_sec << ", " << a << "] not found!" << endl;
        }
        destruction_complex.at(D)->multiply_add_conjugate(x.psi,y,z);
      }
      // cout << "creation at " << a.label << '\n' << x.psi << '\n' << y << "\n\n"; // tempo
    }
    else{ // destruction
      if(Hamiltonian_format == H_FORMAT::H_format_onthefly){
        for(uint32_t I=0; I<B->dim; ++I){
          auto P = Destroy(a.orb+n_orb*a.spin, I, *B, *T);
          if(!get<3>(P)) continue;
          get<1>(P) = get<1>(P)%(2*group->g);
          y[get<0>(P)] += z*group->phaseX<Complex>(get<1>(P))*x.psi[I];
        }
      }
      else{
        destruction_identifier D(x.sec,a);
        if(destruction_complex.find(D)==destruction_complex.end()){
          cout << "destruction operator [" << x.sec << ", " << a << "] not found!" << endl;
        }
        destruction_complex.at(D)->multiply_add(x.psi,y,z);
      }
      // cout << "destruction at " << a.label << '\n' << x.psi << '\n' << y << "\n\n"; // tempo
    }
  }
  return true;
}

