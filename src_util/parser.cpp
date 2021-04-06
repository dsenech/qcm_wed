#include <iomanip>
#include <cstring>
#include <cstdio>

#include "parser.hpp"

using namespace std;


bool parser::no_rewind = false;
int parser::display_accur = 8;

/**
 Looks for the string \a search in the input stream and positions the stream right after the string if found.
 If not, returns the stream at EOF.
 Used to input parameters (located after a keyword) or options that are not mandatory
 */
istream & operator==(istream &flux, const char* search)
{
  string buff;
  
  // resets the stream to the beginning
  if(!parser::no_rewind){
    flux.clear();
    flux.seekg(0);
  }
  
  //recherche la chaine
  while (flux.peek() != EOF)
  {
    flux >> buff;
    if(buff.compare(search)==0){
      return flux;
    }
    flux.ignore(numeric_limits<streamsize>::max(), '\n');
  }
  flux.setstate(ios::failbit);
  return flux;
}






/**
 Looks for the string \a search in the input stream and positions the stream right after the string if found.
 If not, returns the stream at EOF.
 Used to input parameters (located after a keyword) or options that are not mandatory
 */
istream & operator==(istream &flux, const string &search)
{
  string buff;
  
  // resets the stream to the beginning
  if(!parser::no_rewind){
    flux.clear();
    flux.seekg(0);
  }
  
  //recherche la chaine
  while (flux.peek() != EOF)
  {
    flux >> buff;
    if(buff.compare(search)==0){
      return flux;
    }
    flux.ignore(numeric_limits<streamsize>::max(), '\n');
  }
  flux.setstate(ios::failbit);
  return flux;
}






/**
 Looks for the string \a search in the input stream and positions the stream right after the string if found.
 If not, terminate the program.
 Used to input mandatory parameters (located after a keyword)
 */
istream & operator>>(istream &flux, const char* search)
{
  string buff;
  
  // resets the stream to the beginning
  if(!parser::no_rewind){
    flux.clear();
    flux.seekg(0);
  }
  
  
  //recherche la chaine
  while (flux.peek() != EOF)
  {
    flux >> buff;
    if(buff.compare(search)==0){
      return flux;
    }
    flux.ignore(numeric_limits<streamsize>::max(), '\n');
  }
  cerr << "ATTENTION: " << search << " not found in input file.";
  exit(1);
  return flux;
}








/**
 Looks for the string \a search in the input stream and positions the stream right after the string if found.
 If not, terminate the program.
 Used to input mandatory parameters (located after a keyword)
 */
istream & operator>>(istream &flux, const string &search)
{
  string buff;
  
  // resets the stream to the beginning
  if(!parser::no_rewind){
    flux.clear();
    flux.seekg(0);
  }
  
  
  //recherche la chaine
  while (flux.peek() != EOF)
  {
    flux >> buff;
    if(buff.compare(search)==0){
      return flux;
    }
    flux.ignore(numeric_limits<streamsize>::max(), '\n');
  }
  cerr << "ATTENTION: " << search << " not found in input file.";
  exit(1);
  return flux;
}








/**
 Looks for the string \a search in the input stream and positions the stream right after the string if found.
 If not, returns the stream at EOF.
 Does not start from the beginning of file, but from current position, contrary to the operators defined above.
 */
istream & parser::find_next(istream &flux, const char* search)
{
  string buff;
  
  //recherche la chaine
  while (flux.peek() != EOF)
  {
    flux >> buff;
    if(buff.compare(search)==0){
      return flux;
    }
    flux.ignore(numeric_limits<streamsize>::max(), '\n');
  }
  cerr << "ATTENTION: " << search << " not found in input file.";
  exit(1);
  return flux;
}

bool parser::has_only_digits(const string &s){
  return (s.find_first_not_of("0123456789") == string::npos);
}

/**
 reads a vector of string from a line. Stops at the first comment character (#)
 clears the vector before reading
 */

vector<string> read_strings(istream &s)
{
  vector<string> X;
  string line;

  do{
    getline(s,line);
    istringstream sin(line);
    while(sin.good()){
      string tmp;
      sin >> skipws >> tmp;
      if(tmp[0]=='#' or tmp.size() == 0) break;
      X.push_back(tmp);
    }
  } while(line.length() > 0 and line[0]=='#');
  return X;
}




/**
 searches for a match in a vector of strings and returns true is there is one
 */
bool parser::find_string(const string& s, vector<string>& input)
{
  for(auto& x : input) if(s == x) return true;
  return false;
}



ostream & operator<<(ostream &s, vector<string> &X)
{
  for(auto& x: X) s << x << "  ";
  return s;
}


/**
 returns the suffix number of a string. E.g. "abc_12" yields 12.
 ED WARNING : string argument is modified: the suffix is removed.
 */
int cluster_index_from_string(string& S)
{
  int label = 0;
  size_t pos = S.rfind('_');
  if(pos != string::npos){
    label =  from_string<int>(S.substr(pos+1));
    S.erase(pos);
  }
  return label;
}




void check_name(const string& S)
{
  if(S.rfind('_') != string::npos) qcm_ED_throw("the separator '_' is forbidden within operator names!" );
}





