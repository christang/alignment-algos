#ifndef _ALISTRANDEVAL
#define _ALISTRANDEVAL

#include <iostream>
#include <string>
#include <vector>
#include <list>

#include "ssss_shared_defs.h"

using namespace std;

class Alignment_Strand_Evaluator {

 private:

  // data
  int num_sses;
  bool** SSE_contacts;

  // strand vectors
  vector<int> All_Strands;
  vector<int> Edge_Strands;
  vector<int> Core_Strands;

  // rule vectors
  vector<list<int> > All_Strands_Paired;
  vector<list<int> > No_Missing_Cores;

 public:

  // data

  // FUNCTIONS
  Alignment_Strand_Evaluator();
  ~Alignment_Strand_Evaluator();

  // general
  void load_SSE_contacts( int size, bool** contacts );
  void load_All_Strands( vector<SSE_Data> );
  void determine_rules();
  bool ali_passes_rules( const list<int>& );
  bool list_contains( const list<int>&, int );


  // access

  // debugging

};

#endif  //_ALISTRANDEVAL
