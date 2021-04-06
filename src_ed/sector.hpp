#ifndef sector_h
#define sector_h

#include <cstring>
#include "console.hpp"

namespace HS_mixing {
  const int normal=0;
  const int anomalous=1;
  const int spin_flip=2;
  const int full=3;
  const int up_down=4;
}

void qcm_ED_throw(const std::string& s);

//! represents a sector of the Hilbert space
struct sector{
	const static int not_conserved = 9999; //!< value of N or S if not conserved
	const static int unspecified = 10000; //!< value of N or S or irrep if unspecified
	int N; //!< the total number of particles in the sector (N_1 + N_2).  = not_conserved is particle number is not conserved
	int S; //!< the total spin in the sector (N_1 - N_2). = not_conserved if spin is not conserved.
	size_t irrep; //!< the label of the point group representation used (from 0 to g-1)
	
	sector(): N(0), S(0), irrep(0) {}
	
	sector(int _N, int _S, size_t _irrep): N(_N), S(_S), irrep(_irrep) {}
	
  // returns the number of up spins
  // N = Nup+Ndw, S = Nup-Ndw ---> Nup = (N+S)/2
  int Nup() const{
    if(N==not_conserved or S==not_conserved) return not_conserved;
    else return (N+S)/2;
  }

  // returns the number of down spins
  // N = Nup+Ndw, S = Nup-Ndw ---> Ndw = (N-S)/2
  int Ndw() const{
    if(N==not_conserved or S==not_conserved) return not_conserved;
    else return (N-S)/2;
  }


	/**
   constructor from a string
   */
	sector(const string &str) {
		
		if(str.find("R") == string::npos) irrep = 0;
		if(str.find("S") == string::npos) S = not_conserved;
		if(str.find("N") == string::npos) N = not_conserved;
		
    bool valid = true;
    int nel=0;
    if(str.find("R") == string::npos){
      if(str.find("S") == string::npos){
        if(str.find("N") == string::npos){
          valid = false;
        }
        else{
          nel = sscanf(str.c_str(),"N%d",&N);
          if(nel!=1) valid = false;
        }
      }
      else{
        if(str.find('N') == string::npos){
          nel = sscanf(str.c_str(),"S%d",&S);
          if(nel!=1) valid = false;
        }
        else{
          nel = sscanf(str.c_str(),"N%d:S%d",&N,&S);
          if(nel!=2) valid = false;
        }
      }
    }
    else{
      if(str.find("S") == string::npos){
        if(str.find("N") == string::npos){
          nel = sscanf(str.c_str(),"R%ld",&irrep);
          if(nel!=1) valid = false;
        }
        else{
          nel = sscanf(str.c_str(),"R%ld:N%d",&irrep,&N);
          if(nel!=2) valid = false;
        }
      }
      else{
        if(str.find("N") == string::npos){
          nel = sscanf(str.c_str(),"R%ld:S%d",&irrep,&S);
          if(nel!=2) valid = false;
        }
        else{
          nel = sscanf(str.c_str(),"R%ld:N%d:S%d",&irrep,&N,&S);
          if(nel!=3) valid = false;
          if((S+N)%2) qcm_ED_throw("sector string " + str + " makes no sense: N + S should be even!");
        }
      }
    }
    
    if(!valid){
      qcm_ED_throw("sector string " + str + " does not conform to standard!");
    }
  }
		
	
	
	
	friend std::ostream & operator<<(std::ostream &flux, const sector &s){
		flux << "R" << s.irrep;
		if(s.N != s.not_conserved) flux << ":N" << s.N;
		if(s.S != s.not_conserved) flux << ":S" << s.S;
		return flux;
	}
	
	
	
	
	
	friend std::istream & operator>>(std::istream &flux, sector &s){
		char tmp_str[32];
		flux >> tmp_str;
		
		if(strchr(tmp_str,'R')!=nullptr){
			if(strchr(tmp_str,'S')==nullptr){
				s.S = s.not_conserved;
				if(strchr(tmp_str,'N')==nullptr) s.N = s.not_conserved;
				else sscanf(tmp_str,"R:%ld,N:%d",&s.irrep,&s.N);
			}
			else{
				if(strchr(tmp_str,'N')==nullptr){
					s.N = s.not_conserved;
					sscanf(tmp_str,"R:%ld,S:%d",&s.irrep,&s.S);
				}
				else sscanf(tmp_str,"R:%ld,N:%d,S:%d",&s.irrep,&s.N,&s.S);
			}
		}
		else{
			qcm_ED_throw("sector string " + s.name() + " does not conform to standard!");
		}
		return flux;
	}
	
  
  
  
	string name() const{
		ostringstream s;
		s << *this;
		return s.str();
	}
	
};


namespace std
{
	template<>
	struct less<sector>{
		bool operator()(const sector &x, const sector &y) const{
			if(x.S < y.S) return true;
			else if(x.S > y.S) return false;
			else if(x.N < y.N) return true;
			else if(x.N > y.N) return false;
			else if(x.irrep < y.irrep) return true;
			else return false;
		}
	};
}



bool operator!=(const sector &S1, const sector &S2);
bool operator==(const sector &S1, const sector &S2);
bool operator>(const sector &S1, const sector &S2);
bool operator<(const sector &S1, const sector &S2);




#endif
