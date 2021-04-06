#ifndef state_h
#define state_h

#include <cstdio>

#include "sector.hpp"
#include "Q_matrix_set.hpp"
#include "continued_fraction_set.hpp"
#include "parser.hpp"
#include "ED_basis.hpp"

//! state (e.g. the ground state) from which a contribution to the Green function is computed
template<typename HilbertField>
struct state
{
	sector sec; //!< sector to which the state belongs
	vector<HilbertField> psi; 	//!< the Hilbert-space vector representing the state
	double energy; //!< energy of the state
	double weight; //!< its weight in the density matrix
	
  shared_ptr<Green_function_set> gf;
  shared_ptr<Green_function_set> gf_down;
  
  /**
   default constructor
   */
  state() : sec("R0") {}

  
  
  /**
	 constructor from sector and dimension
	 */
	state(sector _sec, size_t dim) : sec(_sec)
	{
		psi.resize(dim);
	}
  
  
  /**
   constructor from ASCII file
   */
  state(istream& fin, shared_ptr<symmetry_group> group , int mixing, GF_FORMAT GF_solver)
  {
    vector<string> input = read_strings(fin);
    if(input.size()!=3) qcm_ED_throw("failed to read a state from input file. Need sector, energy and weight in header line");
    sec = sector(input[0]);
    energy = from_string<double>(input[1]);
    weight = from_string<double>(input[2]);
    if(GF_solver == GF_format_BL) gf = shared_ptr<Q_matrix_set<HilbertField>>(new Q_matrix_set<HilbertField>(fin, group, mixing));
    else gf = shared_ptr<continued_fraction_set>(new continued_fraction_set(fin, sec, group, mixing));
    if(mixing&HS_mixing::up_down){
      if(GF_solver == GF_format_BL) gf_down = shared_ptr<Q_matrix_set<HilbertField>>(new Q_matrix_set<HilbertField>(fin, group, mixing));
      else gf_down = shared_ptr<continued_fraction_set>(new continued_fraction_set(fin, sec, group, mixing));
    }
  }
  

  /**
  constructor from a continued fraction solution from a previous computation
   */
  state(shared_ptr<symmetry_group> group , int mixing, const double ener, const double _weight, const string& _sec, const vector<vector<double>> &a, const vector<vector<double>> &b)
  : energy(ener), weight(_weight), sec(sector(_sec))
  {
    gf = shared_ptr<continued_fraction_set>(new continued_fraction_set(sec, group, mixing, a, b));
    if(mixing&HS_mixing::up_down){
      // TEMPO : completer pour down. Arguments supplementaires ou integres dans A et B ?
      gf_down = shared_ptr<continued_fraction_set>(new continued_fraction_set(sec, group, mixing, a, b));
    }
  }


  /**
  constructor from a q-matrix solution from a previous computation
   */
  state(shared_ptr<symmetry_group> group , int mixing, const double ener, const double _weight, const string& _sec, const vector<vector<double>> &e, const vector<matrix<HilbertField>> &q)  : energy(ener), weight(_weight), sec(sector(_sec))
  {
    gf = shared_ptr<Q_matrix_set<HilbertField>>(new Q_matrix_set<HilbertField>(sec, group, mixing, e, q));
    if(mixing&HS_mixing::up_down){
      // TEMPO : completer pour down. Arguments supplementaires ou integres dans A et B ?
      gf_down = shared_ptr<Q_matrix_set<HilbertField>>(new Q_matrix_set<HilbertField>(sec, group, mixing, e, q));
    }
  }


  /**
   writing the Green function representation to an ASCII file
   */
  void write(ostream& fout)
  {
    fout << "state\n" << sec << '\t' << energy << '\t' << weight << endl;
    if(gf != nullptr) gf->write(fout);
    if(gf_down != nullptr) gf_down->write(fout);
  }


  /**
   writing the wavefunction to an ASCII file
   */
  void write_wavefunction(ostream& fout, const ED_basis &B)
  {
    fout << "state\n" << sec << '\t' << energy << '\t' << weight << endl;
    if(B.dim <= global_int("dim_max_print")){
      for(int i=0; i<B.dim; i++){
        fout << abs(psi[i])*abs(psi[i]) << '\t' << psi[i] << '\t';
        B.print_state(fout,i);
        fout << '\n';
      }
    }
  }


};

template<typename HilbertField>
std::ostream& operator<<(std::ostream &flux, const state<HilbertField> &x)
{
  flux << "E = " << x.energy << " (" << x.sec << ") weight = " << x.weight << endl;
  return flux;
}



#endif

