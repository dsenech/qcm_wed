#ifndef parser_h
#define parser_h

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <cstdlib>
#include <iomanip>
#include <limits>

#include "console.hpp"

void qcm_ED_throw(const std::string& s);

using namespace std;

istream & operator==(istream &input, const char* search);
istream & operator>>(istream &input, const char* search);
istream & operator==(istream &input, const string &search);
istream & operator>>(istream &input, const string &search);

vector<string> read_strings(istream &s);

namespace parser{
  extern	bool no_rewind;
  extern int display_accur;
  inline void next_line(std::istream &flux){flux.ignore(numeric_limits<streamsize>::max(), '\n');};
  istream & find_next(istream &flux, const char* search);
  bool has_only_digits(const string &s);
  bool find_string(const string& s, vector<string>& input);
}

//! converts a string to a generic type
template<typename T>
T from_string(const string &s){
  istringstream sin(s);
  T x;
  sin >> x;
  if(sin.fail()){
    cerr << "Fatal error in 'from_string()', string '" << s << "' cannot be interpreted as type requested\n";
    exit(1);
  }
  return x;
}

//! converts a generic type to a string
template<typename T>
string to_string(const T &x){
  ostringstream sout;
  sout << x;
  return sout.str();
}

int cluster_index_from_string(string& S);
void check_name(const string& S);

#endif
