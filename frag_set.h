#ifndef _FRAGSET
#define _FRAGSET

#include <iostream>
#include <string>
#include <vector>
#include <list>

#include "ssss_shared_defs.h"
#include "sse_frag_set.h"

using namespace std;

//class Ali_Dist;

class Frag_Set {

 private:

  // DATA
  vector<SSE_Frag_Set> Frag_Columns; // holds all fragments!

  string query_seq, templ_seq;

  int num_sses;
  int query_len, templ_len;

  bool verbose;

 public:

  // DATA

  // FUNCTIONS
  Frag_Set();
  ~Frag_Set();

  // operators
  vector<Ali_Frag*> operator-( Frag_Set& );

  // general
  void clear_all();
  void add_column( SSE_Frag_Set );
  void set_num_sses();
  void activate_terminal_caps();
  void initialize_all_zscores();
  void seed_all_columns( int );
  void count_frag_children();
  float activate_next_best_available_frag();
  bool an_available_frag_exists();
  vector<Ali_Frag*> get_active_frags( int );
  vector<Ali_Frag*> get_available_frags( int );
  vector<Ali_Frag*> export_all_frags();

  bool frags_in_order( int, int, int, int );
  bool frags_in_order( Ali_Frag*, Ali_Frag* );

  // region-focused
  vector<Ali_Frag*> get_active_frags_in_region( int, int, int, int );
  vector<Ali_Frag*> get_available_frags_in_region( int, int, int, int );


  // access
  inline SSE_Frag_Set* get_col( int i ) { return &Frag_Columns[i]; }
  int num_frags_in_sse( int );
  Ali_Frag* get_frag( Frag_ID );
  Ali_Frag* get_frag( int, int );
  inline void set_verbose( bool b ) { verbose = b; }

  // debugging
  //  void print_Active_Frags_matrix() const;
  //  void print_Available_Frags_matrix() const;
  void print_all_frags();
  void print_all_sses();
  void print_all_frag_connections();
  void report_statistics();
  //  void pause();

};

#endif  //_FRAGSET
