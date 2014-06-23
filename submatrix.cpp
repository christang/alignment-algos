/**
 *  Package HmapII
 *  File: submatrix.cpp
 *  Desc: Substitution matrix classes
 *
 *  Created on 3/13/2003
 *  Author: cltang @ honig lab
 *
 */

#include "submatrix.h"
#include <fstream>

// Constructor for blosum formatted matrix

BlosumMatrix::BlosumMatrix (const char* filename) 
{
  // modified from VCTroll (dspetrey @ honig lab)
  
  const int buffer_size = 256;
  char buffer[buffer_size];
  
  ifstream sfile(filename);
  if(!sfile.good()) {
    throw string ("File not found (substitution matrix) ")+filename;
  }		         
  
  while(1) {
    sfile.getline(buffer,buffer_size);
    if(buffer[0]!='#') break;
  }

  // Define residue types.  First line after comments.

  int nletters=0;
  for(int i=0; i<buffer_size; ++i) { 
    if(buffer[i]=='\0')
      break;
    if(buffer[i]!=' ') 
      if(buffer[i]!='\n') {
	alphabet.append(1,buffer[i]);
	++nletters;
      }
  }

  // Read in substitution values.
  
  for(int i=0; i<nletters; i++) {
    sfile >> buffer;
    for(int j=0;j<nletters;j++) {
      sfile >> sub_matrix[alphabet[i]][alphabet[j]];
    }
  }
}

// Print readout of a substitution matrix

ostream& operator<<(ostream& os, 
		    SubstitutionMatrix& m) 
{
  int nletters = m.alphabet.length();
  for(int i=0;i<nletters;i++) {
    for(int j=0;j<nletters;j++) {
      os << m.alphabet[i] << m.alphabet[j] << ":";
      os << m.sub_matrix[m.alphabet[i]][m.alphabet[j]] << endl;
    }
  }
  return os;
}
