//
//  QCM.hpp
//  qcm3
//
//  Created by David Sénéchal on 16-11-22.
//  Copyright © 2016 David Sénéchal. All rights reserved.
//

#ifndef QCM_hpp
#define QCM_hpp

/**
 Interface for using the QCM library
 */

#include <string>
#include <vector>
#include <complex>
#include <map>
#include <functional>
#include "vector3D.hpp"

using namespace std;

struct lattice_matrix_element;

/**
 Interface per se
 */
namespace QCM{
  double Berry_flux(vector<vector3D<double>>& k, int band=0, int label=0);
  double monopole(vector3D<double>& k, double a, int nk, int band, bool rec, int label=0);
  double Potthoff_functional(int label=0);
  double spectral_average(const string& name, const complex<double> w, int label=0);
  int mixing();
  int spatial_dimension();
  map<string,double> parameters(int label=0);
  matrix<complex<double>> cluster_Green_function(size_t i, complex<double> w, bool spin_down = false, int label=0);
  matrix<complex<double>> cluster_self_energy(size_t i, complex<double> w, bool spin_down = false, int label=0);
  matrix<complex<double>> cluster_hopping_matrix(size_t i, bool spin_down = false, int label=0);
  matrix<complex<double>> CPT_Green_function(const complex<double> w, const vector3D<double> &k, bool spin_down = false, int label=0);
  matrix<complex<double>> hybridization_function(complex<double> w, bool spin_down = false, size_t i=0, int label=0);
  matrix<complex<double>> periodized_Green_function(const complex<double> w, const vector3D<double> &k, bool spin_down = false, int label=0);
  matrix<complex<double>> band_Green_function(const complex<double> w, const vector3D<double> &k, bool spin_down = false, int label=0);
  matrix<complex<double>> projected_Green_function(const complex<double> w, bool spin_down = false, int label=0);
  matrix<complex<double>> tk(const vector3D<double> &k, bool spin_down = false, int label=0);
  matrix<complex<double>> V_matrix(const complex<double> w, const vector3D<double> &k, bool spin_down = false, int label=0);
  pair<string,string> properties(int label=0);
  pair<vector<array<double,9>>, vector<array<complex<double>, 11>>> site_and_bond_profile(int label=0);
  size_t Green_function_dimension();
  size_t reduced_Green_function_dimension();
  string git_hash();
  vector<double> Berry_curvature(vector3D<double>& k1, vector3D<double>& k2, int nk, int band=0, bool rec=false, int dir=3, int label=0);
  vector<double> dos(const complex<double> w, int label=0);
  vector<double> momentum_profile(const string& op, const vector<vector3D<double>> &k_set, int label=0);
  vector<matrix<complex<double>>> CPT_Green_function_inverse(const complex<double> w, const vector<vector3D<double>> &k, bool spin_down = false, int label=0);
  vector<matrix<complex<double>>> CPT_Green_function(const complex<double> w, const vector<vector3D<double>> &k, bool spin_down = false, int label=0);
  vector<matrix<complex<double>>> periodized_Green_function(const complex<double> w, const vector<vector3D<double>> &k, bool spin_down = false, int label=0);
  vector<matrix<complex<double>>> band_Green_function(const complex<double> w, const vector<vector3D<double>> &k, bool spin_down = false, int label=0);
  vector<complex<double>> periodized_Green_function_element(int r, int c, const complex<double> w, const vector<vector3D<double>> &k, bool spin_down = false, int label=0);
  vector<matrix<complex<double>>> self_energy(const complex<double> w, const vector<vector3D<double>> &k, bool spin_down = false, int label=0);
  vector<matrix<complex<double>>> tk(const vector<vector3D<double>> &k, bool spin_down = false, int label=0);
  vector<pair<double,string>> ground_state(int label=0);
  vector<pair<string,double>> averages(int label=0, bool print=true);
  vector<pair<vector<double>, vector<double>>> Lehmann_Green_function(vector<vector3D<double>> &k, int band=0, bool spin_down = false, int label=0);
  vector<tuple<string, int, int>> cluster_info();
  vector<vector<double>> dispersion(const vector<vector3D<double>> &k, bool spin_down = false, int label=0);
  void add_cluster(const string &name, const vector3D<int64_t> &cpos, const vector<vector3D<int64_t>> &pos, int ref=0);
  void anomalous_operator(const string &name, vector3D<int64_t> &link, complex<double> amplitude, int band1, int band2, const string& type);
  void density_wave(const string &name, vector3D<int64_t> &link, complex<double> amplitude, int band, vector3D<double> Q, double phase, const string& type);
  void explicit_operator(const string &name, const string &type, const vector<tuple<vector3D<int64_t>, vector3D<int64_t>, complex<double>>> &elem, int tau=1, int sigma=0);
  void global_parameter_init();
  void hopping_operator(const string &name, vector3D<int64_t> &link, double amplitude, int band1, int band2, int tau, int sigma);
  void interaction_operator(const string &name, vector3D<int64_t> &link, double amplitude, int band1, int band2, const string &type);
  void k_integral(int dim, function<void (vector3D<double> &k, const int *nv, double I[])> f, vector<double> &Iv, const double accuracy);
  void new_lattice_model(const string &name, vector<int64_t> &superlattice, vector<int64_t> &lattice);
  void new_model_instance(int label=0);
  void print_model(const string& filename, bool asy_operators=false, bool asy_labels=false, bool asy_band=false, bool asy_neighbors=false, bool asy_working_basis=false);
  void qcm_init();
  void set_basis(vector<double> &basis);
  void set_parameter(const string& name, double value);
  void set_parameters(vector<pair<string,double>>&, vector<tuple<string, double, string>>&);
  void wk_integral(int dim, function<void (Complex w, vector3D<double> &k, const int *nv, double I[])> f, vector<double> &Iv, const double accuracy);
  void Green_function_solve(int label=0);
  void CDMFT_variational_set(vector<string>& varia);
  void CDMFT_host(const vector<double>& freqs, const vector<double>& weights, int label=0);
  double CDMFT_distance(const vector<double>& p, int label=0);
};

#endif /* QCM_hpp */







