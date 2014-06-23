/**
 *  Package HMAP2.1
 *  File: sflags.h
 *  Desc: Handle suboptimal region flags.
 *
 *  Created on 11/26/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_SFLAGS
#define _HMAP2_SFLAGS

#include <iostream>
#include <valarray>

#include "sequence.h"

using namespace std;

class SuboptFlags : public Sequence<SequenceElem*>
{
  typedef Sequence<SequenceElem*> Base;
public:
  SuboptFlags (bool f, size_t len);
  inline bool operator[] (unsigned int i) const { return subopt[i]; }
  void append (const string& s);
  void append (const char* cs);
  void Set(unsigned int,bool);
  size_t size () const;
  void setString();
private:
  valarray<bool> subopt;
  unsigned int last_elem;
};

//istream& operator>> (istream& is, SuboptFlags& s);
//ostream& operator<< (ostream& os, SuboptFlags& s);

#endif  //_HMAP2_SFLAGS
