/**
 *  Package HmapII
 *  File: submatrix.h
 *  Desc: Substitution matrix classes
 *
 *  Created on 3/13/2003
 *  Author: cltang @ honig lab
 *
 */

#ifndef _SUBMATRIX
#define _SUBMATRIX

#include <map>
#include <string>

using namespace std;

class SubstitutionMatrix  {
  
  friend ostream& operator<< (ostream& os, SubstitutionMatrix& m);
  
protected:
  
  string alphabet;
  map <char,map<char,float> > sub_matrix;
  
public:
  
  // check if letter exists in alphabet of submatrix
  inline bool hasLetter (char x) const {
    return alphabet.find (x) != string::npos;
  };
  
  // score the alignment of letters a to b
  inline float score (char a, char b) const { 
    return sub_matrix.find(a)->second.find(b)->second; 
  };
  
};

class BlosumMatrix : public SubstitutionMatrix {
  
public:
  
  BlosumMatrix (const char* filename);
  
};

#endif // _SUBMATRIX
