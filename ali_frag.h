#ifndef _ALIFRAG
#define _ALIFRAG

#include <iostream>
#include <sstream>
#include <vector>

#include "ssss_shared_defs.h"

using namespace std;

class Ali_Frag {

 private:

  // DATA
  int t_sse_beg; // first residue of the template sse
  int t_sse_end; // last residue of the template sse
  int t_core_beg; // first template residue of the fragment
  int t_core_end; // last template residue of the fragment
  int qt_shift; // q - t, uniquely defines the diagonal of this fragment

  vector<Frag_Connection> next_frags; // possible connections to downstream frags

  int status;  // active(1), available(0), redundant(-1)

 public:

  // DATA
  int sse_id;  // the template SSE in the fragment
  int frag_id; // unique identifier for each alignment of the query to that template SSE

  bool frag_is_N_terminal, frag_is_C_terminal;

  float score;
  float z_score;

  unsigned long long int num_children;

  // FUNCTIONS
  Ali_Frag( int, int, int, float, bool, bool );
  Ali_Frag( int, int, int, int, int, float, bool, bool );
  ~Ali_Frag();

  // OPERATORS
  inline bool operator<( const Ali_Frag& af ) { return( score < af.ss() ); }

  // GENERAL
  void make_connection( Frag_ID, int, int, float );

  // ACCESS
  inline int core_t0() const { return t_core_beg; };
  inline int core_t1() const { return t_core_end; };
  inline int core_q0() const { return t_core_beg + qt_shift; };
  inline int core_q1() const { return t_core_end + qt_shift; };
  inline int sse_t0() const { return t_sse_beg; };
  inline int sse_t1() const { return t_sse_end; };
  inline int sse_q0() const { return t_sse_beg + qt_shift; };
  inline int sse_q1() const { return t_sse_end + qt_shift; };
  inline int q( int t ) const { return t + qt_shift; }
  inline int qt() const { return qt_shift; }
  inline int core_len() const { return ( t_core_end - t_core_beg ) + 1; };
  inline int sse_len() const { return ( t_sse_end - t_sse_beg ) + 1; };
  inline float ss() const { return score; }; // ss = "sim score"
  inline float zs() const { return z_score; }; // zs = "z score"
  inline void set_score( float s ) { score = s; };

  inline int get_status() const { return status; }
  inline bool is_active() const { return status == 1; }
  inline bool is_available() const { return status == 0; }
  inline bool is_redundant() const { return status == -1; }
  inline void make_active() { status = 1; }
  inline void make_available() { status = 0; }
  inline void make_redundant() { status = -1; }

  Frag_ID get_id();

  inline int num_next() { return next_frags.size(); };
  //  inline void add_next( Ali_Frag* af ) { next.push_back( af ); };
  inline Frag_Connection get_next( int i ) { return next_frags[i]; };
  inline Frag_Connection get_last_next() { return next_frags.back(); };
  inline void clear_next() { next_frags.clear(); };

  inline unsigned long long int get_num_children() { return num_children; };
  inline void set_num_children( unsigned long long int n ) { num_children = n; };

  // DEBUGGING
  void print( ostream& = cerr ) const;
  void print( string, string, ostream& = cerr ) const;
  string print2() const;
  string print2( string, string ) const;
  void print( string, string, int, int, ostream& = cerr ) const;
  void print_all_info( string, string ) const;
  void print_one_line( ostream& = cerr );
  void print_one_line( string, string, ostream& = cerr );

};

#endif  //_ALIFRAG
