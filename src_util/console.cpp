#include <cstring>
#ifdef _OPENMP
  #include <omp.h>
#endif

#include "console.hpp"
#include "global_parameter.hpp"



/**
 splits a string by delimiter and returns a vector of substrings
 */
vector<std::string> split_string(const string &s, char delim) {
  vector<string> elems;
  std::stringstream ss;
  ss.str(s);
  string item;
  while (getline(ss, item, delim)) elems.push_back(item);
  return elems;
}





/**
 string Console routines
 */

double console::precision=1e-6;

/**
 Prints a banner-like message with message \a s, padded with character \a c
 @param c character for padding
 @param s string to print
 @param fout output stream
 */
void console::banner(const char c, const char s[128], ostream &fout)
{
	size_t i,l,l2,l3;
	
  l = (int)strlen(s);
  if(l==0){
    for(i=0; i<80; ++i) fout << c;
    fout << endl;
    return;
  }
  if(l < 76){
    l2 = 80 - l - 4;
    if(l2%2) l3 = l2/2 + 1;
    else l3 = l2/2;
    l2 = l2/2;
    fout << "\n";
    for(i=0; i<l2; ++i) fout << c;
    fout << "  ";
    fout << s;
    fout << "  ";
    for(i=0; i<l3; ++i) fout << c;
    fout << endl;
  }
  else{
    for(i=0; i<80; ++i) fout << c;
    fout << endl;
    fout << s << endl;
    for(i=0; i<80; ++i) fout << c;
    fout << endl;
  }
}





/**
 Prints a banner-like message with message \a s, padded with character \a c
 @param c character for padding
 @param s string to print
 @param fout output stream
 */
void console::banner(const char c, const string &s, ostream &fout)
{
	size_t i,l,l2,l3;
	
  l = (int)s.size();
  if(l==0){
    for(i=0; i<80; ++i) fout << c;
    fout << endl;
    return;
  }
  if(l < 76){
    l2 = 80 - l - 4;
    if(l2%2) l3 = l2/2 + 1;
    else l3 = l2/2;
    l2 = l2/2;
    fout << "\n";
    for(i=0; i<l2; ++i) fout << c;
    fout << "  ";
    fout << s;
    fout << "  ";
    for(i=0; i<l3; ++i) fout << c;
    fout << endl;
  }
  else{
    for(i=0; i<80; ++i) fout << c;
    fout << endl;
    fout << s << endl;
    for(i=0; i<80; ++i) fout << c;
    fout << endl;
  }
}





void console::normal_stop(){
	banner('/',"Program exited normally");
	exit(0);
}

void qcm_throw(const std::string& s)
{
    console::banner('*', s, std::cerr);
    throw(s);
}

void qcm_ED_throw(const std::string& s)
{
    console::banner('*', s, std::cerr);
    throw(s);
}

