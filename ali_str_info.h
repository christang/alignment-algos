#ifndef _ALI_STR_INFO
#define _ALI_STR_INFO

#include <iostream>
#include <vector>
#include <string>

#include "ssss_shared_defs.h"

using namespace std;

class Ali_Str_Info {

 private:

  // DATA
  int templ_len, query_len;
  int num_templ_SSEs;
  string templ_seq, query_seq;

  float** sims;
  float** cb_dists;
  bool** templ_contacts;
  bool* query_predicted_loops;
  vector<SSE_Data> templ_SSEs;
  int *TSR_to_N, *TSR_to_C;

 public:

  // DATA

  // FUNCTIONS
  Ali_Str_Info();
  ~Ali_Str_Info();

  // OPERATORS

  // SETUP
  void load_seq_lengths( int, int );
  void load_seq_strings( string, string );
  void load_sims( float** );
  void load_cb_dists( float** );
  void load_contacts( bool** );
  void load_query_predicted_loops( bool* );
  void load_SSE_Data( vector<SSE_Data> );
  void load_TSRs( int*, int* );

  // ACCESS
  inline float get_sim( int i, int j ) { return sims[i][j]; }
  inline float get_cb_dist( int i, int j ) { return cb_dists[i][j]; }
  inline bool get_contact( int i, int j ) { return templ_contacts[i][j]; }
  inline bool get_query_predicted_loop( int i ) { return query_predicted_loops[i]; }
  inline SSE_Data get_SSE_Data( int i ) { return templ_SSEs[i]; }
  inline int get_TSR_to_N( int i ) { return TSR_to_N[i]; }
  inline int get_TSR_to_C( int i ) { return TSR_to_C[i]; }
  inline int get_num_templ_SSEs() { return num_templ_SSEs; }
  inline string get_templ_seq() { return templ_seq; }
  inline string get_query_seq() { return query_seq; }
  inline int get_templ_len() { return templ_len; }
  inline int get_query_len() { return query_len; }

  // DEBUGGING

};

#endif  //_ALI_STR_INFO
