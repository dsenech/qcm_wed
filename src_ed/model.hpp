#ifndef model_h
#define model_h

#include <map>
#include <set>
#include <string>

#include "ED_basis.hpp"
#include "state.hpp"
#include "destruction_operator.hpp"


#define MIN_GAP 1.0e-6

struct Hermitian_operator;

//! class for a cluster model, independent of the precise value of its parameters
struct model
{
  bool is_closed; //!< true if operators can no longer be added, as the first instance of the model was created
  bool is_factorized; //!< true if the Hamiltonian can be expressed as H = K_up\otimes 1 + 1\otimes K_down + V
  bool has_complex_HS; //!< true if the model requires a complex Hilbert space
  map<destruction_identifier, shared_ptr<destruction_operator<Complex>>> destruction_complex; //!< set of destruction operators
  map<destruction_identifier, shared_ptr<destruction_operator<double>>> destruction; //!< set of destruction operators
  map<sector, shared_ptr<ED_mixed_basis>> basis; //!< list of bases
  map<sector, shared_ptr<ED_factorized_basis>> factorized_basis; //!< list of bases
  map<string, shared_ptr<Hermitian_operator>> term; //!< list of operators in the Hamiltonian, by name
  shared_ptr<symmetry_group> group; //!< contains data on symmetry operations
  size_t n_bath; //!< number of sites considered as 'bath' (no Green function computation for these)
  size_t n_orb; //!< total number of sites (=L+nb)
  size_t n_sites; //!< number of sites in the cluster per se (for Green function computations)
  string name; //!< name of the system (i.e. cluster name)
  vector<bool> in_bath; //! indicates whether an orbital (from 0 to 2*n_orb) is in bath or not
  vector<vector<vector<symmetric_orbital>>> sym_orb;

  model(const string &_name, const size_t _n_orb, const size_t _n_bath, const vector<vector<int>> &gen, bool bath_irrep=false);
  void add_chemical_potential();
  bool create_or_destroy(int pm, const symmetric_orbital &a, state<double> &x, vector<double> &y, double z);
  bool create_or_destroy(int pm, const symmetric_orbital &a, state<Complex> &x, vector<Complex> &y, Complex z);
  void print(ostream& fout);
  void build_HS_objects(const sector& GS_sector, bool is_complex);
  void build_HS_objects_GF(const sector& GS_sector, int mixing, bool spin_down, bool is_complex);
  vector<pair<sector,int>> needed_sectors_GF(const sector& GS_sector, int mixing, bool spin_down);
  void print_graph(const vector<vector<double>> &pos);
};





#endif
  
