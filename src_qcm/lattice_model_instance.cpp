#include <iomanip>
#include <fstream>

#include "lattice_model_instance.hpp"
#include "integrate.hpp"
#include "parser.hpp"
#include "console.hpp"
#include "QCM.hpp"
#ifdef _OPENMP
  #include <omp.h>
#endif

extern shared_ptr<parameter_set> param_set;

string strip_at(const string& s)
{
  string name = s;
  size_t pos1 = name.rfind('@');
  size_t pos2 = name.rfind('_');
  if(pos1 != string::npos){
    int label =  from_string<int>(name.substr(pos2+1));
    name.erase(pos1);
    name += '_' + to_string(label);
  }
  return name;
}



//==============================================================================
/** 
 default constructor
 @param _model [in] : lattice model
 @param _params [in] : values of the model parameters, as a map
 @param _sectors [in] : sectors to look for the ground state 
 @param _label [in] : label given to the instance 
 */

lattice_model_instance::lattice_model_instance(shared_ptr<lattice_model> _model, const map<string, double>& _params, const vector<string>& _sectors, int _label)
: label(_label), model(_model), params(_params), sectors(_sectors)
{
  n_clus = model->clusters.size();
  vector<map<string,double>> cluster_values(n_clus);
  gs_solved = false;
  gf_solved = false;
  average_solved = false;
  SEF_solved = false;
  E_pot = 0.0;
  E_kin = 0.0;
  console::level = (int)global_int("verbose");
  if(sectors.size() == 0) qcm_throw("target sectors were not specified!");
  
  model->close_model();
  
  for(auto& x : params){
    string name = x.first;
    auto P = model->name_and_label(name, true);
    // cout << "---->  " << P.first << '\t' << P.second << endl; // TEMPO
    if(P.second){
      if(x.second != 0.0){
        cluster_values[P.second-1][P.first] = x.second;
        // cout << P.second-1 << '\t' << P.first << '\t' << x.second << endl; // TEMPO
      }
    }
    else model->term.at(name)->is_active = true;
  }
  for(size_t i=0; i<n_clus; i++){
    // cout << "cluster no " << i+1 << " name = " << model->clusters[i].name << endl; // TEMPO
    if(cluster_values[i].size() == 0) qcm_throw("cluster "+to_string(i)+" has no nonzero operators");
    ED::new_model_instance(model->clusters[i].name, cluster_values[i], sectors[i], label*n_clus+i);
    model->clusters[i].mixing = ED::mixing(label*n_clus+i);
  }
	if(model->GF_offset.size() == 0) model->post_parameter_consolidate(label);
}



//==============================================================================
/** 
 finds the ground states of all clusters
 @returns a vector of (double, string) giving the ground state energy and sector for each cluster
 */
vector<pair<double,string>> lattice_model_instance::ground_state()
{
  if(gs_solved) return gs;
  static bool first_time = true;

  clus_ave.resize(n_clus);
  gs.resize(n_clus);
	GS_energy.resize(n_clus+1);
	GS_energy[0] = 0.0;


  if(first_time){
    first_time = false;
    for(size_t i = 0; i<model->inequiv.size(); i++){
      auto I = model->inequiv[i];
      gs[I] = ED::ground_state_solve(n_clus*label+I);
      clus_ave[I] = ED::cluster_averages(n_clus*label+I);
    }
  }
  else{
    #ifdef _OPENMP
      extern int max_num_threads;
      int nt = max_num_threads;
      if(nt > model->inequiv.size()) nt = model->inequiv.size();
      #pragma omp parallel for num_threads(nt)
    #endif  
    for(size_t i = 0; i<model->inequiv.size(); i++){
      auto I = model->inequiv[i];
      gs[I] = ED::ground_state_solve(n_clus*label+I);
      clus_ave[I] = ED::cluster_averages(n_clus*label+I);
    }
  }
  for(size_t i = 0; i<n_clus; i++){
    if(model->clusters[i].ref != i){
      gs[i] = gs[model->clusters[i].ref];
      clus_ave[i] = clus_ave[model->clusters[i].ref];
    }
  }
  for(size_t i = 0; i<n_clus; i++){
    GS_energy[i+1] = gs[i].first;
    GS_energy[0] += gs[i].first;
    for(auto& x : clus_ave[i]) get<0>(x) += '_' + to_string(i+1);
  }
  gs_solved = true;
  
  return gs;
}



//==============================================================================
/** 
 building the k-independent part of the one-body matrices
 */
void lattice_model_instance::build_H()
{
  H.set_size(model->dim_GF);
  if(model->mixing == HS_mixing::up_down) H_down.set_size(model->dim_GF);
  for(auto& x : model->term){
    lattice_operator& op = *x.second;
    for(auto& e : op.GF_elem) H(e.r, e.c) += e.v*params[op.name];
    if(model->mixing == HS_mixing::up_down) for(auto& e : op.GF_elem_down) H_down(e.r, e.c) += e.v*params[op.name];
  }
  
  // building the cluster one-body matrices
  build_cluster_H();
  
  gf_solved = true;
}

//==============================================================================
/** 
 finds the Green function representation of all clusters
 */
void lattice_model_instance::Green_function_solve()
{
  static bool first_time = true;
  if(gf_solved) return;
  if(!gs_solved) ground_state();
  if(first_time){
    first_time = false;
    for(size_t i = 0; i<model->inequiv.size(); i++){
      auto I = model->inequiv[i];
      ED::Green_function_solve(n_clus*label+I);
    }
  }
  else{
    #ifdef _OPENMP
      extern int max_num_threads;
      int nt = max_num_threads;
      if(nt > model->inequiv.size()) nt = model->inequiv.size();
      #pragma omp parallel for num_threads(nt)
    #endif
    for(size_t i = 0; i<model->inequiv.size(); i++){
      auto I = model->inequiv[i];
      ED::Green_function_solve(n_clus*label+I);
    }
  } 
  build_H();
}



//==============================================================================
/** 
 building the cluster one-body matrices
 */
void lattice_model_instance::build_cluster_H()
{
  Hc.block.assign(n_clus, matrix<Complex>());
  for(size_t i = 0; i<n_clus; i++){
    Hc.block[i].assign(cluster_hopping_matrix(i));
  }
  Hc.set_size();
  if(model->mixing == HS_mixing::up_down){
    Hc_down.block.assign(n_clus, matrix<Complex>());
    for(size_t i = 0; i<n_clus; i++){
      Hc_down.block[i].assign(cluster_hopping_matrix(i, true));
    }
    Hc_down.set_size();
  }
}


//==============================================================================
/** 
 returns the cluster Green function for cluster # i at frequency w
 @param i [in] index of the cluster (1 to the number of clusters)
 @param w [in] complex frequency
 @param spin_down [in] true if we are asking for the spin down part (mixing = 4)
 @returns a complex-valued matrix containing the cluster Green function
 */
matrix<complex<double>> lattice_model_instance::cluster_Green_function(size_t i, complex<double> w, bool spin_down)
{
  if(i >= model->clusters.size()) qcm_throw("cluster label out of range");
  if(!gf_solved) Green_function_solve();
  int I = n_clus*label+model->clusters[i].ref;
  matrix<Complex> g = ED::Green_function(w, spin_down, I);
  matrix<Complex> G;
  int mix = model->clusters[i].mixing;
  if(model->mixing == mix){
    return g;
  }

  // combinaisons 0:2, 0:4, 1:3, 1:5
  else if(((model->mixing&1) == 0 and (mix&1) == 0) or ((model->mixing&1)==1 and (mix&1)==1)) {
    G = upgrade_cluster_matrix(model->mixing, mix, g);
  }
  // combinaisons 0:1, 0:3, 0:5, 2:3
  else if((model->mixing&1) == 1 and (mix&1) == 0){
    auto gm = ED::Green_function(-w, false, I);
    G = upgrade_cluster_matrix_anomalous(model->mixing, mix, g, gm);
  }
  else{
    qcm_throw("undefined mixing combinations in cluster_Green_function()");
  }
  return G;
}

//==============================================================================
/** 
 returns the cluster Green function for cluster # i at frequency w
 @param i [in] index of the cluster (1 to the number of clusters)
 @param w [in] complex frequency
 @param spin_down [in] true if we are asking for the spin down part (mixing = 4)
 @returns a complex-valued matrix containing the cluster self-energy
 */
matrix<complex<double>> lattice_model_instance::cluster_self_energy(size_t i, complex<double> w, bool spin_down)
{
  if(i >= model->clusters.size()) qcm_throw("cluster label out of range");
  if(!gf_solved) Green_function_solve();
  int I = n_clus*label+model->clusters[i].ref;
  matrix<Complex> g = ED::self_energy(w, spin_down, I);
  int mix = model->clusters[i].mixing;
  if(model->mixing == mix) return g;

  // combinaisons 0:2, 0:4, 1:3, 1:5
  else if(((model->mixing&1) == 0 and (mix&1) == 0) or ((model->mixing&1)==1 and (mix&1)==1)) {
    return upgrade_cluster_matrix(model->mixing, mix, g);
  }
  // combinaisons 0:1, 0:3, 0:5, 2:3
  else if((model->mixing&1) == 1 and (mix&1) == 0){
    auto gm = ED::self_energy(-w, false, I);
    return upgrade_cluster_matrix_anomalous(model->mixing, mix, g, gm);
  }
  else{
    qcm_throw("undefined mixing combinations in cluster_self_energy()");
    matrix<complex<double>> empty;
    return empty;
  }
}

//==============================================================================
/** 
 returns the hybridization function for cluster # i at frequency w
 @param i [in] index of the cluster (1 to the number of clusters)
 @param w [in] complex frequency
 @param spin_down [in] true if we are asking for the spin down part (mixing = 4)
 @returns a complex-valued matrix containing the cluster hybridization function
 */
matrix<complex<double>> lattice_model_instance::hybridization_function(size_t i, complex<double> w, bool spin_down)
{
  if(i >= model->clusters.size()) qcm_throw("cluster label out of range");
  int I = n_clus*label+model->clusters[i].ref;
  matrix<Complex> g = ED::hybridization_function(w, spin_down, I);
  int mix = model->clusters[i].mixing;
  if(model->mixing == mix) return g;

  // combinaisons 0:2, 0:4, 1:3, 1:5
  else if(((model->mixing&1) == 0 and (mix&1) == 0) or ((model->mixing&1)==1 and (mix&1)==1)) {
    return upgrade_cluster_matrix(model->mixing, mix, g);
  }
  // combinaisons 0:1, 0:3, 0:5, 2:3
  else if((model->mixing&1) == 1 and (mix&1) == 0){
    auto gm = ED::hybridization_function(-w, false, I);
    return upgrade_cluster_matrix_anomalous(model->mixing, mix, g, gm);
  }
  else{
    qcm_throw("undefined mixing combinations in hybridization_function()");
    matrix<complex<double>> empty;
    return empty;
  }
}


//==============================================================================
/** 
 returns the cluster hopping matrix for cluster # i
 @param i [in] index of the cluster (1 to the number of clusters)
 @param spin_down [in] true if we are asking for the spin down part (mixing = 4)
 @returns a complex-valued matrix containing the cluster non interacting matrix (hopping + pairing)
 */
matrix<complex<double>> lattice_model_instance::cluster_hopping_matrix(size_t i, bool spin_down)
{
  if(i >= model->clusters.size()) qcm_throw("cluster label out of range");
  int I = n_clus*label+model->clusters[i].ref;
  matrix<Complex> g = ED::hopping_matrix(spin_down, I);
  int mix = model->clusters[i].mixing;

  if(model->mixing == mix) return g;

  // combinations 0:2, 0:4, 1:3, 1:5
  else if(((model->mixing&1) == 0 and (mix&1) == 0) or ((model->mixing&1)==1 and (mix&1)==1)) {
    return upgrade_cluster_matrix(model->mixing, mix, g);
  }
  // combinations 0:1, 0:3, 0:5, 2:3
  else if((model->mixing&1) == 1 and (mix&1) == 0){
    return upgrade_cluster_matrix_anomalous(model->mixing, mix, g, g);
  }
  else{
    qcm_throw("undefined mixing combinations in cluster_hopping_matrix()");
    matrix<complex<double>> empty;
    return empty;
  }
}


//==============================================================================
/** 
 computes the self-energy for the Green function object G
 @param G [in] Green function object
 */
void lattice_model_instance::cluster_self_energy(Green_function& G)
{
  if(!gf_solved) Green_function_solve();
	G.sigma.block.assign(n_clus, matrix<Complex>());
  for(size_t i = 0; i<n_clus; i++){
    auto S =  cluster_self_energy(i, G.w, G.spin_down);
    G.sigma.block[i].assign(S);
  }
  G.sigma.set_size();
}


//==============================================================================
/** 
 computes cluster Green function and puts it in the structure G
 @param w [in] complex frequency
 @param sig [in] if true, computes also the self-energy
 @param spin_down [in] if true, computes the spin-down part
 @returns a Green_function object
 */
Green_function lattice_model_instance::cluster_Green_function(Complex w, bool sig, bool spin_down)
{
  if(!gf_solved) Green_function_solve();
	Green_function G;
	G.w = w;
	G.spin_down = spin_down;
	G.G.block.resize(n_clus);

  for(size_t i = 0; i<n_clus; i++){
    G.G.block[i].assign(cluster_Green_function(i, w, spin_down));
  }
  G.G.set_size();
  if(sig){
    G.sigma.block.resize(n_clus);
    for(size_t i = 0; i<n_clus; i++){
      G.sigma.block[i].assign(cluster_self_energy(i, w, spin_down));
    }
    G.sigma.set_size();
  }
  
  if(model->bath_exists){
    G.gamma.block.resize(n_clus);
    for(size_t i = 0; i<n_clus; i++){
      G.gamma.block[i].assign(hybridization_function(i, w, spin_down));
    }
    G.gamma.set_size();
  }
	return G;
}


//==============================================================================
/** 
 Calculates the self-energy functional (SEF) = Potthoff functional
 @returns the Potthoff functional
 */
double lattice_model_instance::Potthoff_functional()
{
  if(SEF_solved) return omega;
  if(!gf_solved) Green_function_solve();
	
	double omega_clus=0.0;
	for(size_t i=0; i<n_clus; i++){
    int I = model->clusters[i].ref;
    omega_clus += ED::Potthoff_functional(label*n_clus+I);
  }
	
	vector<double> Iv(1);
	Iv[0] = 0.0;
  
  // lambda function
  auto F = [this] (Complex w, vector3D<double> &k, const int *nv, double *I) {SEF_integrand(w, k, nv, I);};
  QCM::wk_integral(model->spatial_dimension, F, Iv, global_double("accur_SEF"));
	double omega_trace = -Iv[0];
	
	
  //-------------------------------------------------------------------------------
	// contribution of the diagonal terms
	double omega_diag=0.0;
	for(size_t i=0; i<model->dim_GF; i++){
		if(model->is_nambu_index[i]) omega_diag += realpart(H(i,i)-Hc(i,i));
		else omega_diag -= realpart(H(i,i)-Hc(i,i));
	}
	
	if(model->mixing == HS_mixing::up_down){
		for(size_t i=0; i<model->dim_GF; i++) omega_diag -= realpart(H_down(i,i)-Hc_down(i,i));
	}
	if(model->mixing == HS_mixing::normal) omega_diag *= 2.0;
	
	omega_diag *= 0.5;
	
  //-------------------------------------------------------------------------------
  // Putting it all together
  if(model->n_mixed == 4){ // If Nambu doubling, then your really have to divide by 2.
    omega_diag *= 0.5;
    omega_trace *= 0.5;
  }
  omega = omega_trace - omega_diag + omega_clus;


  //-------------------------------------------------------------------------------
	// We want the energy density (or per atom), so we divide by the number of cluster orbitals
	omega /= model->sites.size();
  SEF_solved = true;

  // writing to file
  static bool first_print = true;
  ofstream fout("sef.tsv",ios::app);
  print_info();
  if(first_print){
    fout << line_info_names << endl;
    first_print = false;
  }
  fout << line_info_values << endl;
  
	return(omega);
}




//==============================================================================
/** 
 Integrand of the SEF integral
 @param w [in] complex frequency
 @param k [in] wavevector
 @param nv number of values in the integrand
 @param I values of these components
 */
void lattice_model_instance::SEF_integrand(Complex w, vector3D<double> &k, const int *nv, double *I){
	
	matrix<Complex> VG(model->dim_GF);
	
  check_signals();
	Green_function G = cluster_Green_function(w, false);
	Green_function_k K(G,k);
  set_V(K);
	G.G.mult_left(K.V, VG);
	VG.add(Complex(-1.0,0));
	I[0] = log(abs(VG.determinant()));
	
	if(model->mixing == HS_mixing::up_down){
		VG.zero();
    Green_function G_down = cluster_Green_function(w, false, true);
		Green_function_k K_down(G_down,k);
		set_V(K_down);
		G_down.G.mult_left(K_down.V, VG);
		VG.add(Complex(-1.0,0));
		double I_down = log(abs(VG.determinant()));
    I[0] += I_down;
	}
	else if(model->mixing == HS_mixing::normal) I[0] *= 2.0;
}


//==============================================================================
/**
 Prints model_instance parameters to a stream
  @param out [in] output stream
  @param format [in] print format, as described in the enumeration print_format
 */
void lattice_model_instance::print_parameters(ostream& out, print_format format)
{
  bool print_all = global_bool("print_all");

  out << setprecision((int)global_int("print_precision"));
  if(param_set == nullptr){
    switch(format){
      case print_format::names_noave :
        for(auto& x : params) out << x.first << '\t';
        break;
      case print_format::names :
        for(auto& x : params) out << x.first << '\t' << "ave_" << x.first << '\t';
        for(size_t i = 0; i<n_clus; i++){
          // if(model->clusters[i].ref != i) continue;
          for(auto& x : clus_ave[i]) out  << "ave_" << strip_at(get<0>(x)) << '\t';
        }
        for(size_t i = 0; i<n_clus; i++){
          // if(model->clusters[i].ref != i) continue;
          for(auto& x : clus_ave[i]) out  << "var_" << strip_at(get<0>(x)) << '\t';
        }
        break;
      case print_format::values_noave :
        for(auto& x : params) out << chop(x.second, 1e-10) << '\t';
        break;
      case print_format::values :
        for(auto& x : params) out << chop(x.second, 1e-10) << '\t';
        for(auto& x : ave) out << chop(x.second, 1e-10) << '\t';
        for(size_t i = 0; i<n_clus; i++){
          // if(model->clusters[i].ref != i) continue;
          for(auto& x : clus_ave[i]) out << chop(get<1>(x), 1e-10) << '\t';
        }
        for(size_t i = 0; i<n_clus; i++){
          // if(model->clusters[i].ref != i) continue;
          for(auto& x : clus_ave[i]) out << chop(get<2>(x), 1e-10) << '\t';
        }
        break;
      case print_format::namesvalues :
        for(auto& x : params) out << strip_at(x.first) << " = " << chop(x.second, 1e-10) << ", ";
        for(size_t i = 0; i<n_clus; i++){
          // if(model->clusters[i].ref != i) continue;
          for(auto& x : clus_ave[i]) out << "ave_" << strip_at(get<0>(x)) << " = " << chop(get<1>(x), 1e-10) << ", ";
        }
        out << endl;
        break;
      default:
        break;
    }
  }
  else{
    switch(format){
      case print_format::names_noave :
        for(auto& x : params) if(param_set->is_dependent(x.first) ==  false or print_all==true) out << x.first << '\t';
        break;
      case print_format::names :
        for(auto& x : params) if(param_set->is_dependent(x.first) ==  false or print_all==true) out << x.first << '\t';
        for(auto& x : ave) out << "ave_" << x.first << '\t';
        for(size_t i = 0; i<n_clus; i++){
          if(model->clusters[i].ref != i) continue;
          for(auto& x : clus_ave[i]) out  << "ave_" << strip_at(get<0>(x)) << '\t';
        }
        for(size_t i = 0; i<n_clus; i++){
          if(model->clusters[i].ref != i) continue;
          for(auto& x : clus_ave[i]) out  << "var_" << strip_at(get<0>(x)) << '\t';
        }
        break;
      case print_format::values_noave :
        for(auto& x : params) if(param_set->is_dependent(x.first) ==  false or print_all==true) out << chop(x.second, 1e-10) << '\t';
        break;
      case print_format::values :
        for(auto& x : params) if(param_set->is_dependent(x.first) ==  false or print_all==true) out << chop(x.second, 1e-10) << '\t';
        for(auto& x : ave) out << chop(x.second, 1e-10) << '\t';
        for(size_t i = 0; i<n_clus; i++){
          if(model->clusters[i].ref != i) continue;
          for(auto& x : clus_ave[i]) out << chop(get<1>(x), 1e-10) << '\t';
        }
        for(size_t i = 0; i<n_clus; i++){
          if(model->clusters[i].ref != i) continue;
          for(auto& x : clus_ave[i]) out << chop(get<2>(x), 1e-10) << '\t';
        }
        break;
      case print_format::namesvalues :
        for(auto& x : params) if(param_set->is_dependent(x.first) ==  false or print_all==true)  out << strip_at(x.first) << " = " << chop(x.second, 1e-10) << ", ";
        for(size_t i = 0; i<n_clus; i++){
          if(model->clusters[i].ref != i) continue;
          for(auto& x : clus_ave[i]) out << "ave_" << strip_at(get<0>(x)) << " = " << chop(get<1>(x), 1e-10) << ", ";
        }
        out << endl;
        break;
      default:
        break;
    }
  }
}



//==============================================================================
/**
 Prints instance info on an internal line
 */
void lattice_model_instance::print_info()
{
  std::ostringstream ostr1;
  std::ostringstream ostr2;
  int print_precision = (int)global_int("print_precision");
  ostr1 << "model\t";
  ostr2 << model->name << '\t';
  if(SEF_solved){
    ostr1 << "omega\t";
    ostr2 << setprecision(print_precision) << omega << '\t';
  }
  if(average_solved){
    ostr1 << "E_kin\t";
    ostr2 << setprecision(print_precision) << E_kin << '\t';
    if(global_bool("potential_energy")){
      ostr1 << "E_pot\t";
      ostr2 << setprecision(print_precision) << E_pot << '\t';
    }
  }
  ground_state();
  print_parameters(ostr1, print_format::names);
  print_parameters(ostr2, print_format::values);
  for(size_t i = 0; i<n_clus; i++){
    if(model->clusters[i].ref != i) continue;
    ostr1 << "E0_" << i+1 << "\tsector_"  << i+1 << '\t';
    ostr2 << setprecision(print_precision) << gs[i].first << '\t' << gs[i].second << '\t';
  }
  ostr1 << "githash\t";
  ostr2 << QCM::git_hash() << '\t';
  ostr1 << "githash_ED\t";
  ostr2 << ED::git_hash() << '\t';
  line_info_names = ostr1.str();
  line_info_values = ostr2.str();
}



//==============================================================================
/** 
 Produces the Lehmann form of an array of Green functions at fixed frequency 
 @param k [in] array of wavevectors
 @param band [in] band (=orbital) label
 @param spin_down [in] true if the spin-down part is to be produced (mixing = 4)
 @returns for each wavevector, a pair of arrays giving the location of the pole and the residue
 */
vector<pair<vector<double>, vector<double>>> lattice_model_instance::Lehmann_Green_function(vector<vector3D<double>> &k, int band, bool spin_down)
{
  if(band >= model->n_mixed*model->n_band) qcm_throw("the band number is out of range in Lehmann_Green_function");
  vector<pair<vector<double>, vector<double>>> res(k.size());
  auto G = cluster_Green_function(Complex(0., 1.0), false, spin_down);
  vector<double> Lambda;
  block_matrix<Complex> QB;
  size_t nclus = model->clusters.size();
  for(size_t c=0; c<nclus; c++){
    auto q = ED::qmatrix(spin_down, nclus*label+c);
    Lambda.insert(Lambda.end(), q.first.begin(), q.first.end());
    matrix<Complex> Q(model->GF_dims[c], q.first.size());
    Q.v = q.second;
    QB.block.push_back(Q);
  }
  QB.set_size();
  size_t m = Lambda.size();
  const double threshold = 1e-6;
	
  for(int i=0; i<k.size(); i++){
    console::message(3, "k = "+to_string<vector3D<double>>(k[i]));
    res[i].first.resize(m);
    res[i].second.resize(m);
    Green_function_k K(G, k[i]);
	  matrix<Complex> Qk(QB.r,QB.c);
	  matrix<Complex> Lk(m);
	  matrix<Complex> Uk(m);
    set_V(K, true);
	  QB.simil(Lk,K.V);
    Lk.add_to_diagonal(Lambda);
	  Lk.eigensystem(res[i].first,Uk);
	  QB.mult_left(Uk, Qk);
	
    vector<Complex> qk(QB.r);
    vector<Complex> psi(model->n_band);
    for(size_t a=0; a<m; ++a){
      Qk.extract_column(a,qk);
      auto psi = model->periodize(K.k, qk);
      res[i].second[a] = abs(psi[band]*psi[band])*model->Lc;
    }
    // cleaning up
    vector<double>& e = res[i].first;
    for(size_t a=1; a<m; ++a){
		  if(abs(e[a]-e[a-1]) < 1e-6){
			  e[a] = (e[a]+e[a-1])/2;
			  res[i].second[a] += res[i].second[a-1];
			  res[i].second[a-1] = 0.0;
		  }
    }
    vector<double> ep, rp;
    ep.reserve(m);
    rp.reserve(m);
    for(size_t a=1; a<m; ++a){
      if(res[i].second[a] > threshold){
        ep.push_back(res[i].first[a]);
        rp.push_back(res[i].second[a]);
      }
    }
    res[i].first = ep;
    res[i].second = rp;
	}
  return res;
}



//==============================================================================
/** 
 upgrades a cluster Green function g from its native mixing clus_mix to the 
 lattice model's mixing latt_mix. No additional anomalous component required. 
 @param latt_mix [in] mixing state of the lattice model
 @param clus_mix [in] mixing state of the cluster
 @param g [in] the cluster Green function 
 @returns the upgraded Green function
 */
matrix<Complex> lattice_model_instance::upgrade_cluster_matrix(int latt_mix, int clus_mix, matrix<Complex> &g)
{
  int d = g.r;
  if(clus_mix == HS_mixing::normal and latt_mix == HS_mixing::spin_flip){
    matrix<Complex> G(d*2);
    g.move_sub_matrix(d, d, 0, 0, 0, 0, G);
    g.move_sub_matrix(d, d, 0, 0, d, d, G, 1.0);
    return G;
  }
  else if(clus_mix == HS_mixing::normal and latt_mix == HS_mixing::up_down){
    return g;
  }
  else if(clus_mix == HS_mixing::anomalous and latt_mix == HS_mixing::full){
    matrix<Complex> G(d*2);
    int h = d/2;
    g.move_sub_matrix(h, h, 0, 0, 0, 0, G);
    g.move_sub_matrix(h, h, 0, 0, h, h, G);
    g.move_sub_matrix(h, h, h, h, d, d, G);
    g.move_sub_matrix(h, h, h, h, d+h, d+h, G);
    g.move_sub_matrix(h, h, h, 0, d+h, 0, G);
    g.move_sub_matrix(h, h, h, 0, d, h, G, -1.0);
    g.move_sub_matrix(h, h, h, 0, h, d, G, -1.0);
    g.move_sub_matrix(h, h, h, 0, 0, d+h, G);
    return G;
  }
  else{
    qcm_throw("impossible combination of cluster and lattice mixings");
    matrix<Complex> empty;
    return empty;
  }
}



//==============================================================================
/** 
 upgrades a cluster Green function g and its Nambu version gm from its native mixing clus_mix to the 
 lattice model's mixing latt_mix.
 @param latt_mix [in] mixing state of the lattice model
 @param clus_mix [in] mixing state of the cluster
 @param g [in] the cluster Green function 
 @param gm [in] Nambu version of the cluster Green function 
 @returns the upgraded Green function
 */
matrix<Complex> lattice_model_instance::upgrade_cluster_matrix_anomalous(int latt_mix, int clus_mix, matrix<Complex> &g, matrix<Complex> &gm)
{
  int d = g.r;
  if(clus_mix == HS_mixing::normal and latt_mix == HS_mixing::anomalous){
    matrix<Complex> G(d*2);
    g.move_sub_matrix(d, d, 0, 0, 0, 0, G);
    gm.move_sub_matrix_transpose(d, d, 0, 0, d, d, G, -1.0);
    return G;
  }
  else if(clus_mix == HS_mixing::normal and latt_mix == HS_mixing::full){
    matrix<Complex> G(d*4);
    g.move_sub_matrix(d, d, 0, 0, 0, 0, G);
    g.move_sub_matrix(d, d, 0, 0, d, d, G, 1.0);
    gm.move_sub_matrix_transpose(d, d, 0, 0, 2*d, 2*d, G, -1.0);
    gm.move_sub_matrix_transpose(d, d, 0, 0, 3*d, 3*d, G, -1.0);
    return G;
  }
  else if(clus_mix == HS_mixing::spin_flip and latt_mix == HS_mixing::full){
    matrix<Complex> G(d*2);
    g.move_sub_matrix(d, d, 0, 0, 0, 0, G);
    gm.move_sub_matrix_transpose(d, d, 0, 0, d, d, G, -1.0);
    return G;
  }
  else{
    qcm_throw("impossible combination of cluster and lattice mixings");
    matrix<complex<double>> empty;
    return empty;
  }
}