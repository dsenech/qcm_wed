/**
Controls output to the screen of various messages, depending on a preset verbosity level
Performs other formatting-related tasks
 */


#ifndef console_h
#define console_h

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <complex>

#include "types.hpp"

using namespace std;

#define SHORT_DISPLAY 7
#define NORMAL_DISPLAY 9
#define LONG_DISPLAY 14

vector<std::string> split_string(const string &s, char delim);

/**
 Chops a real number \a x if its absolute value is smaller than the limit \a c
 */
inline double chop(double x, double c=1e-6){return (fabs(x)<c) ? 0 : x;}

void check_signals();
void qcm_ED_throw(const std::string& s);
void qcm_ED_catch(const std::string& s);
void qcm_throw(const std::string& s);
void qcm_catch(const std::string& s);


/**
 Chops a complex number \a x separately for its real and imaginary parts
 */
inline Complex chop(Complex x, double c=1e-6){return Complex(fabs(x.real()) < c ? 0 : x.real(),fabs(x.imag()) < c ? 0:x.imag());}

namespace console{
	extern double precision;
	extern int level;
	
	void message(int level, const string &str);
	void banner(const char c, const char s[128], std::ostream &fout = std::cout);
	void banner(const char c, const string &s, ostream &fout = std::cout);
	void normal_stop();
}

#endif
