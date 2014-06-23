#ifndef _SKELALI
#define _SKELALI

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <list>

#include "ssss_shared_defs.h"

#include "ali_frag.h"
#include "ali_str_info.h"
#include "frag_set.h"

using namespace std;

class Skel_Ali {

 private:

  // data
  vector<Frag_Connection> connections; // list of 'connections' between fragments

  float score;
  float shift;
  float param;

  int num_aligned_residues;  // number of aligned residues in the skeleton alignment
  float SSE_CO;
  vector<int> contacting_residues;
  int num_contacting_residues;

  int templ_len;

  Ali_Str_Info* Str_Data;
  Frag_Set* Frags;

 public:

  // data

  // FUNCTIONS
  Skel_Ali( Frag_Connection, Ali_Str_Info*, Frag_Set*, int );
  Skel_Ali( Ali_Str_Info*, Frag_Set* );
  Skel_Ali( Skel_Ali* );
  ~Skel_Ali();

  // general
  void add_connection( Frag_Connection );
  void calc_skel_SSE_CO();
  void update_contacted_residues();
  int get_last_templ_res_idx();

  vector<Res_Pair> export_vrp() const;

  // operators
  inline bool operator<( const Skel_Ali& sa ) { return( param < sa.get_param() ); }
  bool operator==( const Skel_Ali& );
  bool operator!=( const Skel_Ali& );

  // access
  inline bool last_frag_is_C_terminal() { return get_frag( connections.back().next_frag )->frag_is_C_terminal; }
  inline int num_connections() const { return connections.size(); }
  inline Frag_Connection get_connection( int i ) const { return connections[i]; }
  inline Frag_Connection get_last_connection() const { return connections.back(); }

  inline int get_contact_status( int i ) const { return contacting_residues[i]; }
  inline int get_num_contacts() const { return num_contacting_residues; }

  inline int get_templ_len() const { return templ_len; };
  //  inline int get_Cterm_templ_res() const { return Frags->get_frag( connections.back().next_frag )->core_t1(); }

  inline Ali_Frag* get_frag( Frag_ID f ) const { return Frags->get_frag( f ); }

  inline Frag_Set* get_frag_set() const { return Frags; }
  inline Ali_Str_Info* get_str_info() const { return Str_Data; }

  inline void  set_param( float s ){ param = s; }
  inline float get_param() const { return param; }

  inline void  set_score( float s ){ score = s; }
  inline float get_score() const { return score; }

  inline void  set_shift( float s ){ shift = s; }
  inline float get_shift() const { return shift; }

  inline int get_num_aligned() const { return num_aligned_residues; }
  inline float get_contact_order() const { return SSE_CO; }


  list<int> get_sse_id_list();

  // debugging
  void print( string, string, int, ostream& = cerr ) const;
  void print() const;
  string print2() const;
  void print_plain() const;

};

#endif  //_SKELALI
