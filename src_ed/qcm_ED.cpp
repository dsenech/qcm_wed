/** @file
 This file is the main header for the library, included in other projects
 */

#include <iostream>
#include <fstream>
#include <memory>
#include <tuple>
#include "qcm_ED.hpp"
#include "model_instance.hpp"
#include "one_body_operator.hpp"
#include "anomalous_operator.hpp"
#include "interaction_operator.hpp"
#include "Hund_operator.hpp"
#include "Heisenberg_operator.hpp"
#ifdef _OPENMP
  #include <omp.h>
#endif
using namespace ED;

// global variables
map<string, shared_ptr<model>> models;
map<size_t, shared_ptr<model_instance_base>> model_instances;

double max_gap; //!< maximum excitation energy compatible with minimum_weight and temperature


//------------------------------------------------------------------------------


// template specializations exposed
// template class matrix_element<double>;
// template class matrix_element<Complex>;

// template class model_instance<double>;
// template class model_instance<Complex>;

namespace ED{
  
  void new_model(const string &name, const size_t L, const size_t nb, const vector<vector<int>> &gen, bool bath_irrep)
  {
    console::level = (int)global_int("verbose");
    if(models.find(name) != models.end()){
      qcm_ED_throw("The name "+name+" has already been used for a model");
    }
    if(nb>32) qcm_ED_throw("number of bath orbitals exceeds limits or negative!");
    if(L>1024) qcm_ED_throw("number of cluster sites  exceeds limits or negative!");
    
    models[name] = make_shared<model>(name, L+nb, nb, gen, bath_irrep);
    
    max_gap = -log(global_double("minimum_weight"))*global_double("temperature");
    if(max_gap < MIN_GAP) max_gap = MIN_GAP;
  }
  
  
  
  
  
  void new_operator(const string &model_name, const string &_name, const string &_type, const vector<matrix_element<double>> &elements)
  {
    if(!elements.size()) return;
    if(models.find(model_name) == models.end())
      qcm_ED_throw("The model "+model_name+" is not defined. Check spelling.");
    shared_ptr<model> M = models.at(model_name);
    if(M->is_closed) qcm_ED_throw("model " + model_name + " has already been instantiated and is closed for modifications");
    if(_type == "one-body") M->term[_name] = make_shared<one_body_operator<double>>(_name, M, elements);
    else if(_type == "anomalous") M->term[_name] = make_shared<anomalous_operator<double>>(_name, M, elements);
    else if(_type == "interaction") M->term[_name] = make_shared<interaction_operator>(_name, M, elements);
    else if(_type == "Hund") M->term[_name] = make_shared<Hund_operator>(_name, M, elements);
    else if(_type == "Heisenberg") M->term[_name] = make_shared<Heisenberg_operator>(_name, M, elements);
    else if(_type == "X") M->term[_name] = make_shared<Heisenberg_operator>(_name, M, elements, 'X');
    else if(_type == "Y") M->term[_name] = make_shared<Heisenberg_operator>(_name, M, elements, 'Y');
    else if(_type == "Z") M->term[_name] = make_shared<Heisenberg_operator>(_name, M, elements, 'Z');
    else console::message(0, "ED_WARNING : type of operator "+_name+" is not yet implemented");
  }
  
  
  void new_operator(const string &model_name, const string &_name, const string &_type, const vector<matrix_element<Complex>> &elements)
  {
    if(!elements.size()) return;
    shared_ptr<model> M = models.at(model_name);
    if(M->is_closed) qcm_ED_throw("model " + model_name + " has already been instantiated and is closed for modifications");
    if(_type == "one-body") M->term[_name] = make_shared<one_body_operator<Complex>>(_name, M, elements);
    else if(_type == "anomalous") M->term[_name] = make_shared<anomalous_operator<Complex>>(_name, M, elements);
    else console::message(0, "ED_WARNING : type of operator "+_name+" is not yet implemented");
  }
  
  
  
  bool exists(const string &model_name, const string &name)
  {
    if(models.find(model_name) == models.end()) qcm_ED_throw("model "+model_name+" does not exist!");
    model& M = *models.at(model_name);
    if(M.term.find(name) != M.term.end()) return true;
    else return false;
  }
  
  
  
  
  pair<size_t,size_t>  model_size(const string &name)
  {
    if(models.find(name) == models.end()) qcm_ED_throw("model "+name+" does not exist!");
    return {models.at(name)->n_sites, models.at(name)->n_bath};
  }
  
  
  
  
  
  void new_model_instance(const string &model_name, map<string, double> &param, const string &sec, size_t label)
  {
    bool need_complex = false;
    if(models.find(model_name) == models.end()) qcm_ED_throw("The model "+model_name+" is not defined and so no model instance based on it is allowed.");
    model& mod = *models.at(model_name);
    if(mod.is_closed == false) mod.is_closed = true;
    
    // first, remove values associated with non existent
    auto it = param.begin();
    while(it != param.end()){
      if(mod.term.find(it->first) == mod.term.end()){
        // console::message(5, "ED WARNING : operator "+it->first+" in instance "+to_string(label)+" does not exist in cluster model "+mod.name);
        param.erase(it++);
      }
      else ++it;
    }
    
    // need to know whether the instance is complex or real
    if(model_instances.find(label) != model_instances.end()) need_complex = model_instances[label]->complex_Hilbert;
    else{
      for(auto& v : param){
        if(v.second != 0 and mod.term.at(v.first)->is_complex){
          need_complex = true;
          break;
        }
      }
    }
    // decides whether the sector set requires a complex Hilbert space
    if(mod.group->has_complex_irrep) need_complex = true;
    
    if(model_instances.find(label) != model_instances.end()) model_instances[label].reset();
    if(need_complex) model_instances[label] = make_shared<model_instance<Complex>>(label, models.at(model_name), param, sec);
    else model_instances[label] = make_shared<model_instance<double>>(label, models.at(model_name), param, sec);
  }
  
  
  
  int mixing(size_t label)
  {
    return model_instances.at(label)->mixing;
  }
  
  
    
  bool complex_HS(size_t label)
  {
    return model_instances.at(label)->complex_Hilbert;
  }
  


  pair<double, string> ground_state_solve(size_t label)
  {
    if(model_instances.find(label) == model_instances.end()) qcm_ED_throw("The cluster instance label "+to_string(label)+" is out of range.");
    
    pair<double, string> gs;
    gs = model_instances.at(label)->low_energy_states();
    return gs;
  }
  
  
  
  vector<tuple<string,double,double>>  cluster_averages(size_t label)
  {
    if(model_instances.find(label) == model_instances.end()) qcm_ED_throw("The cluster instance label "+to_string(label)+" is out of range.");
    
    auto& M = model_instances.at(label);
    if(!M->gs_solved) M->low_energy_states();
    return M->averages;
  }
  
  
  
  
  void Green_function_solve(size_t label)
  {
    if(model_instances.find(label) == model_instances.end()) qcm_ED_throw("The label "+to_string(label)+" is out of range.");
    model_instances.at(label)->Green_function_solve();
  }
  
  
  
  matrix<complex<double>> Green_function(const Complex &z, bool spin_down, const size_t label, bool blocks)
  {
    return model_instances.at(label)->Green_function(z, spin_down, blocks);
  }

  

  matrix<complex<double>> Green_function_average(bool spin_down, const size_t label)
  {
    return model_instances.at(label)->Green_function_average(spin_down);
  }

  
  

  matrix<complex<double>> self_energy(const Complex &z, bool spin_down, const size_t label)
  {
    return model_instances.at(label)->self_energy(z, spin_down);
  }
  
  
  matrix<complex<double>> hopping_matrix(bool spin_down, const size_t label)
  {
    return model_instances.at(label)->hopping_matrix(spin_down);
  }
  matrix<complex<double>> hopping_matrix_full(bool spin_down, const size_t label)
  {
    return model_instances.at(label)->hopping_matrix_full(spin_down);
  }
  vector<tuple<int,int,double>> interactions(const size_t label)
  {
    return model_instances.at(label)->interactions();
  }

  
  
  matrix<complex<double>> hybridization_function(const Complex w, bool spin_down, const size_t label)
  {
    return model_instances.at(label)->hybridization_function(w, spin_down);
  }
  
  
  
  vector<complex<double>> susceptibility(const string &op, const vector<Complex> &w, const size_t label)
  {
    if(model_instances.find(label) == model_instances.end()) qcm_ED_throw("The label "+to_string(label)+" is out of range.");
    auto& M = model_instances.at(label);
    if(M->the_model->term.find(op) == M->the_model->term.end()) qcm_ED_throw("Operator "+op+" is not defined in the model.");
    return model_instances.at(label)->susceptibility(M->the_model->term.at(op), w);
  }
  
  
  
  vector<pair<double,double>> susceptibility_poles(const string &op, const size_t label){
    if(model_instances.find(label) == model_instances.end()) qcm_ED_throw("The label "+to_string(label)+" is out of range.");
    auto& M = model_instances.at(label);
    if(M->the_model->term.find(op) == M->the_model->term.end()) qcm_ED_throw("Operator "+op+" is not defined in the model.");
    return model_instances.at(label)->susceptibility_poles(M->the_model->term.at(op));
  }
  
  
  
  double Potthoff_functional(const size_t label){
    model_instance_base& M = *model_instances.at(label);

    if(!M.is_correlated) qcm_throw("The Potthoff functional cannot be computed in the noninteracting case!");

    if(!M.gf_solved) M.Green_function_solve();
    return M.SEF_bath + M.E0;
  }
  
  
  
  size_t Green_function_dimension(size_t label)
  {
    return model_instances.at(label)->dimension();
  }
  
  
  
  double tr_sigma_inf(const size_t label)
  {
    return model_instances.at(label)->tr_sigma_inf();
  }
  
  
  void print_models(ostream& fout)
  {
    for(auto& x : models) x.second->print(fout);
    console::banner('=', "model instances",fout);
    for(auto& x : model_instances) x.second->print(fout);
  }
  
  
  
  
  double fidelity(const string& model_name, map<string, double> &param1, map<string, double> &param2, const string &sec)
  {
    bool need_complex = false;
    if(models.find(model_name) == models.end()) qcm_ED_throw("The model "+model_name+" is not defined and so no model instance based on it is allowed.");
    model& mod = *models.at(model_name);
    if(mod.is_closed == false) mod.is_closed = true;
    
    // first, remove values associated with non existent operators
    auto it = param1.begin();
    while(it != param1.end()){
      if(mod.term.find(it->first) == mod.term.end()){
        console::message(0, "ED WARNING : operator "+it->first+" does not exist in cluster model "+mod.name);
        param1.erase(it++);
      }
      else ++it;
    }
    it = param2.begin();
    while(it != param2.end()){
      if(mod.term.find(it->first) == mod.term.end()){
        console::message(0, "ED WARNING : operator "+it->first+" does not exist in cluster model "+mod.name);
        param2.erase(it++);
      }
      else ++it;
    }
    
    // need to know whether the instance is complex or real
    for(auto& v : param1){
      if(v.second != 0 and mod.term.at(v.first)->is_complex){
        need_complex = true;
        break;
      }
    }
    if(need_complex){
      auto I1 = model_instance<Complex>(0, models.at(model_name), param1, sec);
      auto I2 = model_instance<Complex>(1, models.at(model_name), param2, sec);
      return I1.fidelity(I2);
    }
    else{
      auto I1 = model_instance<double>(0, models.at(model_name), param1, sec);
      auto I2 = model_instance<double>(1, models.at(model_name), param2, sec);
      return I1.fidelity(I2);
    }
  }
  
  
  
  void write_instance(ostream& fout, int label)
  {
    if(model_instances.find(label) == model_instances.end()) qcm_ED_throw("The label "+to_string(label)+" is out of range.");
    auto& M = model_instances.at(label);
    M->write(fout);
  }
  
  
  
  void read_instance(istream& fin, int label)
  {
    if(model_instances.find(label) == model_instances.end()) qcm_ED_throw("The label "+to_string(label)+" is out of range.");
    auto& M = model_instances.at(label);
    M->read(fin);
  }

  string git_hash()
  {
#ifdef GITHASH
    return string(GITHASH);
#endif
    return "";
  }
  
  
  pair<vector<double>, vector<complex<double>>> qmatrix(bool spin_down, const size_t label)
  {
    if(model_instances.find(label) == model_instances.end()) qcm_ED_throw("The label "+to_string(label)+" is out of range.");
    auto& M = model_instances.at(label);
    if(!M->complex_Hilbert){
      return dynamic_pointer_cast<model_instance<double>>(M)->qmatrix(spin_down);
    }
    return dynamic_pointer_cast<model_instance<complex<double>>>(M)->qmatrix(spin_down);
  }

  
  pair<vector<double>, vector<complex<double>>> hybridization(bool spin_down, const size_t label)
  {
    if(model_instances.find(label) == model_instances.end()) qcm_ED_throw("The label "+to_string(label)+" is out of range.");
    auto& M = model_instances.at(label);
    if(!M->complex_Hilbert){
      dynamic_pointer_cast<model_instance<double>>(M)->hybridization(spin_down);
    }
    return dynamic_pointer_cast<model_instance<complex<double>>>(M)->hybridization(spin_down);
  }

  string print_wavefunction(const size_t label)
  {
    if(model_instances.find(label) == model_instances.end()) qcm_ED_throw("The label "+to_string(label)+" is out of range.");
    auto& M = model_instances.at(label);
    ostringstream sout;
    M->print_wavefunction(sout);
    return sout.str();
  }

  pair<string, vector<matrix_element<Complex>>> matrix_elements(const string& model_name, const string& op_name)
  {
    if(models.find(model_name) == models.end())
      qcm_ED_throw("The model "+model_name+" is not defined. Check spelling.");
    shared_ptr<model> M = models.at(model_name);
    if(M->term.find(op_name) == M->term.end())
      qcm_ED_throw("operator "+op_name+" is not defined in model . Check spelling.");
    return {M->term[model_name]->type(), M->term[model_name]->matrix_elements()};
  }
}

