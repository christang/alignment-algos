#ifndef _SSEFRAGSET
#define _SSEFRAGSET

#include <iostream>
#include <fstream>
#include <vector>
#include "math.h"

#include "ali_frag.h"

using namespace std;

class SSE_Frag_Set {

 private:

  // DATA

  // flags for frag status
  int active_flag; // = 1
  int available_flag; // = 0
  int redundant_flag; // = -1

 public:

  // DATA
  vector<Ali_Frag> Frags;

  int sse_id;
  int ss_type; // helix = 329, strand = 330
  int t0;
  int t1;
  int qt_shift_lo;
  int qt_shift_hi;
  int sse_len;
  int query_len;
  int templ_len;
  int min_aligned_residues;

  float average_frag_score;
  float standev_frag_score;

  // FUNCTIONS
  SSE_Frag_Set( int, int, int, int, int,
		int, int, int,
		vector<Ali_Frag>, int );

  ~SSE_Frag_Set();

  // OPERATORS

  // GENERAL
  void find_biggest_gap( int&, int&, int& );
  void fill_gap( int, int );

  vector<Ali_Frag*> get_ordered_frags();
  vector<Ali_Frag*> get_all_frags_qt_sorted();
  vector<Ali_Frag*> get_active_frags();
  vector<Ali_Frag*> get_available_frags();
  float get_highest_available_frag_zscore();
  void activate_top_available_frag();
  void set_frag_zscores();
  void activate_frag( int );
  vector<int> find_available_neighbors( int, int );
  bool an_available_frag_exists() const;

  vector<Ali_Frag*> find_shift_neighbors( float, int );
  int get_frag_status( Ali_Frag* );
  int get_num_active_frags();

  // ACCESS
  inline Ali_Frag* get_frag( int f ) { return &Frags[f]; }

  //PRINT
  void print_frag_scores( ofstream& ) const;
  void print_frag_scores( ostream& ) const;
  void print_sse_info( ostream& = cerr ) const;
  void print_sse_info( string, ostream& = cerr ) const;

  // DEBUGGING
  void get_statistics( int&, int&, float&, int& );
};

#endif  //_SSEFRAGSET
