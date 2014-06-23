/**
 *  Package HMAP2.1
 *  File: sequence.h
 *  Desc: Base class describing sequences.
 *
 *  Created on 9/11/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_SEQUENCE
#define _HMAP2_SEQUENCE

#include <sstream>
#include <string>
#include <vector>

using namespace std;

class SequenceElem {
 
 public:
  int   index;
  char  olc;
  
  SequenceElem ();
  SequenceElem (const SequenceElem& e);
  SequenceElem (int i, char o);
  SequenceElem& operator= (const SequenceElem& e);

  inline bool isHead() { return olc==Head; }
  inline bool isTail() { return olc==Tail; }

  static const char Head;
  static const char Tail;
};

template <class elem_t>
class Sequence : public vector <elem_t> {

 public:
  // seq_length = length of sequence without head/tail characters
  unsigned int seq_length;
  string seq_name;

  // NOTE: below is original version that did not seem to work (AK)
  //  inline char olc(int i) { return vector<elem_t>::at(i).olc; }

  // this is the new version, which assumes we have a vector of pointers
  inline char olc(int i) const { return vector<elem_t>::at(i)->olc; }

  const string* getString() const;
  //  bool same_sequence( Sequence& );
  //  void remove_double_start_ends();

 protected:
  mutable string seq_string;

 private:
  void buildSequenceString() const;

};

template <class elem_t>
const string* Sequence<elem_t>::getString() const
{
  if (seq_string=="") 
     buildSequenceString();
  return &seq_string;
}


/*
template <class elem_t>
bool Sequence<elem_t>::same_sequence( Sequence& seq2 )
{
  unsigned int i1(0), i2(0);
  string s1, s2;

  while( ( i1 < size() ) || ( i2 < seq2.size() ) ) {

    // skip gap characters
    if(      olc( i1 ) == '-' ) { i1++; continue; }
    if( seq2.olc( i2 ) == '-' ) { i2++; continue; }

    if( olc(i1) != seq2.olc(i2) ) { 

      cerr << olc(i1) << "  " << seq2.olc(i2) << endl;  // print the conflicting residues

      // include several characters past the problem
      int n = 0;
      while( ( i1 < size() ) && ( i2 < seq2.size() ) && ( n < 10 ) ) {
	s1.append( 1, olc( i1++ ) );
	s2.append( 1, seq2.olc( i2++ ) );
	n++;
      }

      cerr << "seq 1: " << s1 << endl;  // print the sequences
      cerr << "seq 2: " << s2 << endl;

      return false;
    }
    else {
      s1.append( 1, olc( i1 ) );
      s2.append( 1, seq2.olc( i2 ) );
    }

    i1++;
    i2++;
  }

  return true;
}


template <class elem_t>
void Sequence<elem_t>::remove_double_start_ends()
{
  if( ( olc(0) == '^' ) && ( olc(1) == '^' ) ) {
    erase( begin() );
  }
  if( ( olc( size()-1 ) == '$' ) && ( olc( size()-2 ) == '$' ) ) {
    vector<elem_t>::iterator it = end();
    it--;
    erase( it );
  }
}
*/

template <class elem_t>
void Sequence<elem_t>::buildSequenceString() const
{
  stringbuf buffer ("");
  for (typename vector<elem_t>::const_iterator it=vector<elem_t>::begin(); 
       it!=vector<elem_t>::end(); ++it) {
    buffer.sputc ((*it)->olc);
  }
  seq_string = buffer.str();
}

#endif // _HMAP2_SEQUENCE
