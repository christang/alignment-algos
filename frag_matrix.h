#ifndef _FRAGMATRIX
#define _FRAGMATRIX

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <list>

#include "ssss_shared_defs.h"
#include "frag_set.h"
#include "ali_str_info.h"
#include "sse_frag_set.h"
#include "ali_dist.h"
#include "skel_set.h"

using namespace std;

class Frag_Matrix {

 private:

  // DATA
  string query_seq, templ_seq;

  Frag_Set* Main_FS;
  Ali_Str_Info* Str_Data;
  Ali_Dist* Compare_to_Native;

  int num_sses;
  int query_len, templ_len;

  int max_in_betw_shift;

  int ali_mode; // 0 = full sse ali, 1 = partial

  int min_aligned_residues;

  // debugging
  bool debug;
  stringstream msg;

  // FUNCTIONS

  // connection
  bool connection_is_valid( Frag_Set*, Ali_Frag*, Ali_Frag* );
  bool loop_spans_gap( int, int, int, int );
  void get_connection_info( Frag_Set*, Frag_ID, Frag_ID, int&, int&, float& );
  void get_inter_frag_endpoints( Ali_Frag*, Ali_Frag*, int&, int&, int&, int& );

  // filling the matrix
  void find_biggest_gap_in_Matrix( Frag_Set*, int&, int&, int& );

  // general
  bool it_is_valid_starting_frag( Frag_Set*, Frag_ID, int );
  bool sse_is_native( int, int ) const;
  int find_min_ali_len( int ) const;

 public:

  // DATA

  // FUNCTIONS

  // constructor and destructor
  Frag_Matrix( int,
	       Frag_Set*, Ali_Str_Info*,
	       int, int, Ali_Dist* );
  ~Frag_Matrix();

  // setting up fragment set
  void create_all_fragments( Frag_Set* );
  void find_fragment_connections( Frag_Set* );
  unsigned long long int get_number_of_alis_to_search( Frag_Set* );
  void find_N_terminal_connections( Frag_Set*, int );

  // filling the matrix
  float fill_Frag_Set_by_zscore( Frag_Set* );
  void fill_gaps_in_Frag_Matrix( Frag_Set* );
  bool activate_next_fragment( float&, unsigned long long int&, Frag_Set* );
  //  void add_in_between_frags( Frag_Set*, float );
  void add_in_between_frags2( Frag_Set* );
  void get_in_between_frag_pairs( Frag_Set*, vector<Frag_ID>&, vector<Frag_ID>& );

  // access
  Ali_Frag* get_frag( Frag_Set*, Frag_ID );
  Ali_Frag* get_frag( Frag_Set*, int, int );
  Frag_ID make_Frag_ID( int, int );

  // debugging
  inline void set_debug( bool b ) { debug = b; }
  inline void debug_on() { debug = true; }
  inline void debug_off() { debug = false; }
  inline void print( ostream& os = cerr ) { 
    if( debug ) { fprint( os ); }
    else { msg.clear(); }
  }
  inline void fprint( ostream& os = cerr ) { os << msg.str(); msg.clear(); }
  void report_frag_quality( Frag_Set* );
  void report_full_sse_frag_set_info( Frag_Set* );
  void pause();

};

#endif  //_FRAGMATRIX
