#include "ali_frag.h"

using namespace std;

/*********************************************
* implementation of Frag_Set class *
*********************************************/

// constructor 1
Ali_Frag::Ali_Frag( int t1, int t2, int qt, float s, bool N_term, bool C_term )
{
  // for fragments that are fully aligned to the template sse

  cerr << "ali_frag constructor1: top" << endl;

  t_sse_beg = t1;
  t_sse_end = t2;
  t_core_beg = t_sse_beg;
  t_core_end = t_sse_end;
  qt_shift = qt;
  score = s;
  num_children = 0;
  frag_is_N_terminal = N_term;
  frag_is_C_terminal = C_term;

  cerr << "core_t0: " << core_t0() << endl;
  cerr << "core_t1: " << core_t1() << endl;
  cerr << "core_q0: " << core_q0() << endl;
  cerr << "core_q1: " << core_q1() << endl;

  cerr << "sse_t0: " << sse_t0() << endl;
  cerr << "sse_t1: " << sse_t1() << endl;
  cerr << "sse_q0: " << sse_q0() << endl;
  cerr << "sse_q1: " << sse_q1() << endl;

  cerr << "ali_frag constructor1: end" << endl;
}

// constructor 2
Ali_Frag::Ali_Frag( int t1_sse, int t2_sse, int t1_core, int t2_core, int qt,
		    float s, bool N_term, bool C_term )
{
  // for fragments that are partially aligned to the template sse

  t_sse_beg = t1_sse;
  t_sse_end = t2_sse;
  t_core_beg = t1_core;
  t_core_end = t2_core;
  qt_shift = qt;
  score = s;
  num_children = 0;
  frag_is_N_terminal = N_term;
  frag_is_C_terminal = C_term;
}

Ali_Frag::~Ali_Frag()
{}

// OPERATORS

// GENERAL
void Ali_Frag::make_connection( Frag_ID f_next, int prev_end, int next_beg, float score )
{
  Frag_Connection res;

  Frag_ID tmp_frag;
  tmp_frag.sse_idx = sse_id;
  tmp_frag.frag_idx = frag_id;

  res.prev_frag = tmp_frag;
  res.next_frag = f_next;
  res.prev_end_res_idx = prev_end;
  res.next_beg_res_idx = next_beg;
  res.connection_score = score;

  next_frags.push_back( res );
}

// ACCESS

Frag_ID Ali_Frag::get_id()
{
  Frag_ID res;
  res.sse_idx = sse_id;
  res.frag_idx = frag_id;

  //  cerr << "get_id: sse_id and frag_id: " << sse_id << " and " << frag_id << endl;

  return res;
}


// DEBUGGING
void Ali_Frag::print( ostream& os ) const
{
  os << "Frag: sse id: " << sse_id << ", frag_id: " << frag_id << endl;
  os << "      core: [" << core_t0() << "," << core_q0() << "] - [" << core_t1() << "," << core_q1() << "]" << endl;
  os << "       sse: [" << sse_t0() << "," << sse_q0() << "] - [" << sse_t1() << "," << sse_q1() << "]" << endl;
  os << "        qt: " << qt_shift << endl;
  os << " -- score:   " << score << endl;
  os << " -- z-score: " << z_score << endl;
}


void Ali_Frag::print( string query_seq, string templ_seq, ostream& os ) const
{
  print( os );

  os << "T: ";
  for( int t=t_core_beg; t<=t_core_end; t++ ) {
    os << templ_seq[t];
  }
  os << endl;

  os << "Q: ";
  for( int t=t_core_beg; t<=t_core_end; t++ ) {
    os << query_seq[t+qt_shift];
  }

  os << endl;
}


string Ali_Frag::print2() const
{
  stringstream res;

  res << "Frag: sse id: " << sse_id << ", frag_id: " << frag_id << endl;
  res << "      core: [" << core_t0() << "," << core_q0() 
      << "] - [" << core_t1() << "," << core_q1() << "]" << endl;
  res << "       sse: [" << sse_t0() << "," << sse_q0() 
      << "] - [" << sse_t1() << "," << sse_q1() << "]" << endl;
  res << "        qt: " << qt_shift << endl;
  res << " -- score:   " << score << endl;
  res << " -- z-score: " << z_score << endl;

  return res.str();
}


string Ali_Frag::print2( string query_seq, string templ_seq ) const
{
  stringstream res;

  res << print2();

  res << "T: ";
  for( int t=t_core_beg; t<=t_core_end; t++ ) {
    res << templ_seq[t];
  }
  res << endl;

  res << "Q: ";
  for( int t=t_core_beg; t<=t_core_end; t++ ) {
    res << query_seq[t+qt_shift];
  }

  res << endl;

  return res.str();
}


void Ali_Frag::print( string query_seq, string templ_seq,
		      int t_beg, int t_end, ostream& os ) const
{
  print( os );

  for( int t=t_sse_beg; t<=t_sse_end; t++ ) {
    os << templ_seq[t];
  }
  os << endl;

  for( int t=t_sse_beg; t<=t_sse_end; t++ ) {
    if( ( t >= t_beg ) && ( t <= t_end ) ) {
      os << "|";
    }
    else {
      os << " ";
    }
  }
  os << endl;

  for( int t=t_sse_beg; t<=t_sse_end; t++ ) {
    os << query_seq[t+qt_shift];
  }
  os << endl;
}


void Ali_Frag::print_all_info( string query_seq, string templ_seq ) const
{

  print( query_seq, templ_seq );

  cerr << "# children: " << num_children << endl;
  cerr << "----------------------" << endl;
  cerr << "Next frags:" << endl;
  for( unsigned int i=0; i<next_frags.size(); i++ ) {
    cerr << "Frag (" << next_frags[i].next_frag.sse_idx 
	 << "," << next_frags[i].next_frag.frag_idx << ")" << endl;
  }
  cerr << "----------------------" << endl;
}


void Ali_Frag::print_one_line( ostream& os )
{
  os << qt_shift << ", " << score << ", " << z_score;
}


void Ali_Frag::print_one_line( string templ_seq, string query_seq, ostream& os )
{
  print_one_line( os );

  os << ", "  
     << templ_seq.substr( t_core_beg, 3 ).c_str() << "/" 
     << query_seq.substr( q( t_core_beg ), 3 ).c_str();
}
