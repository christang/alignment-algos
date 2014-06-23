#include "skel_ali.h"

using namespace std;

/*********************************************
* implementation of Skel_Ali class *
*********************************************/

// constructor 1
Skel_Ali::Skel_Ali( Frag_Connection fc, Ali_Str_Info* asi, Frag_Set* fs, int num_ali_init )
{
  connections.clear();
  connections.push_back( fc );

  Str_Data = asi;
  Frags = fs;

  templ_len = Str_Data->get_templ_len();

  score = get_frag( fc.prev_frag )->ss();
  score += fc.connection_score;
  score += get_frag( fc.next_frag )->ss();

  num_aligned_residues = num_ali_init;
  num_aligned_residues += ( ( get_frag( fc.next_frag )->core_t1() - fc.next_beg_res_idx ) + 1 );

  num_contacting_residues = 0;
  SSE_CO = 0;

  contacting_residues.resize( templ_len, -1 );

  // zero out residues in first frag after N-term
  for( int i=fc.next_beg_res_idx; i<=get_frag( fc.next_frag )->core_t1(); i++ ) {
    contacting_residues[i] = 0;
  }

}


// constructor 2
Skel_Ali::Skel_Ali( Ali_Str_Info* asi, Frag_Set* fs )
{
  connections.clear();

  Str_Data = asi;
  Frags = fs;

  templ_len = Str_Data->get_templ_len();

  score = 0;

  num_aligned_residues = 0;

  num_contacting_residues = 0;
  SSE_CO = 0;

  contacting_residues.resize( templ_len, -1 );

}


// constructor 3
Skel_Ali::Skel_Ali( Skel_Ali* sa )
{
  connections.clear();

  for( unsigned int i=0; i<sa->connections.size(); i++ ) {
    connections.push_back( sa->connections[i] );
  }

  score = sa->score;
  num_aligned_residues = sa->get_num_aligned();
  SSE_CO = sa->get_contact_order();
  //  contacts = sa->get_contacts();

  Frags = sa->get_frag_set();
  Str_Data = sa->get_str_info();
  templ_len = sa->get_templ_len();

  contacting_residues.resize( templ_len );
  for( int i=0; i<templ_len; i++ ) {
    contacting_residues[i] = sa->get_contact_status( i );
  }

  num_contacting_residues = sa->get_num_contacts();
}

Skel_Ali::~Skel_Ali()
{}

// GENERAL
void Skel_Ali::add_connection( Frag_Connection fc )
{

  connections.push_back( fc ); // add the new Frag_Connection to the list
  score += get_frag( fc.next_frag )->ss();   // update the score
  //  score += src->get_frag( fc.next_frag )->ss();   // update the score
  score += fc.connection_score;

  if( !get_frag( fc.next_frag )->frag_is_C_terminal ) { // update the coverage for a non-terminal frag
  //  if( !src->get_frag( fc.next_frag )->frag_is_C_terminal ) { // update the coverage for a non-terminal frag
    num_aligned_residues += ( ( fc.prev_end_res_idx - get_frag( fc.prev_frag )->core_t1() ) +
			      ( get_frag( fc.next_frag )->core_t1() - fc.next_beg_res_idx + 1 ) );
  }
  else { // frag is C-terminal
    num_aligned_residues += ( fc.prev_end_res_idx - get_frag( fc.prev_frag )->core_t1() );
  }

  // zero out the contacts of the c-terminal extension of prev_frag
  for( int i=fc.prev_end_res_idx; i>get_frag( fc.prev_frag )->core_t1(); i-- ) {
    contacting_residues[i] = 0;
  }

  // zero out the contacts of the n-terminal extension and core of next_frag
  for( int i=fc.next_beg_res_idx; i<=get_frag( fc.next_frag )->core_t1(); i++ ) {
    contacting_residues[i] = 0;
  }

  update_contacted_residues(); // find new contacts and update 'num_contacting_residues', etc.
}

void Skel_Ali::calc_skel_SSE_CO()
{
  SSE_CO = (float)num_contacting_residues / (float)num_aligned_residues;
}

void Skel_Ali::update_contacted_residues()
{
  // set the values for the most recently added fragment and change the values of previous frags accordingly

  Frag_Connection last = connections.back();

  int t_prev_end = last.prev_end_res_idx;
  int t_prev_core_end = get_frag( last.prev_frag )->core_t1();

  for( int t_new=t_prev_end; t_new>t_prev_core_end; t_new-- ) {

    for( int fc_idx=1; fc_idx<(int)(connections.size()-1); fc_idx++ ) {

      int t_prev_sse_beg = connections[fc_idx-1].next_beg_res_idx;
      int t_prev_sse_end = connections[fc_idx].prev_end_res_idx;

      for( int t_prev=t_prev_sse_beg; t_prev<=t_prev_sse_end; t_prev++ ) {

	//	if( contacts[t_new][t_prev] ) {
	if( Str_Data->get_contact( t_new, t_prev ) ) {

	  if( contacting_residues[t_new] == 0 ) { // this contact not yet recorded in the ith residue in the new frag
	    num_contacting_residues++;
	    contacting_residues[t_new] = 1;
	  }

	  if( contacting_residues[t_prev] == 0 ) { // this contact not yet recorded in the kth residue in the jth frag
	    num_contacting_residues++;
	    contacting_residues[t_prev] = 1;
	  }

	}
      }

    }

  }

  int t_curr_beg = last.next_beg_res_idx;
  int t_curr_core_end = get_frag( last.next_frag )->core_t1();

  for( int t_new=t_curr_beg; t_new<=t_curr_core_end; t_new++ ) {

    for( int fc_idx=1; fc_idx<(int)connections.size(); fc_idx++ ) {

      int t_prev_sse_beg = connections[fc_idx-1].next_beg_res_idx;
      int t_prev_sse_end = connections[fc_idx].prev_end_res_idx;

      for( int t_prev=t_prev_sse_beg; t_prev<=t_prev_sse_end; t_prev++ ) {

	//	if( contacts[t_new][t_prev] ) {
	if( Str_Data->get_contact( t_new, t_prev ) ) {

	  if( contacting_residues[t_new] == 0 ) { // this contact not yet recorded in the ith residue in the new frag
	    num_contacting_residues++;
	    contacting_residues[t_new] = 1;
	  }

	  if( contacting_residues[t_prev] == 0 ) { // this contact not yet recorded in the kth residue in the jth frag
	    num_contacting_residues++;
	    contacting_residues[t_prev] = 1;
	  }

	}
      }

    }

  }

}


int Skel_Ali::get_last_templ_res_idx()
{
  if( connections.size() != 0 ) {
    return Frags->get_frag( connections.back().next_frag )->core_t1();
  }
  else{
    return 0;
  }
}


vector<Res_Pair> Skel_Ali::export_vrp() const
{
  vector<Res_Pair> res;

  // add the endpoints of each connection to 'vrp'
  for( int i=0; i<(int)connections.size(); i++ ) {
    Res_Pair tmp;

    tmp.t = (float)connections[i].prev_end_res_idx;
    tmp.q = (float)( ( get_frag( connections[i].prev_frag ) )->q( connections[i].prev_end_res_idx ) );
    tmp.rel_pos = -2; // default initial value
    res.push_back( tmp );

    tmp.t = (float)connections[i].next_beg_res_idx;
    tmp.q = (float)( ( get_frag( connections[i].next_frag ) )->q( connections[i].next_beg_res_idx ) );
    tmp.rel_pos = -2;
    res.push_back( tmp );
  }

  return res;
}


bool Skel_Ali::operator==( const Skel_Ali& sa )
{

  if( num_connections() != sa.num_connections() ) {
    return false;
  }

  Ali_Frag *af1, *af2; 

  for( int i=0; i<num_connections(); i++ ) {
    af1 = get_frag( get_connection(i).prev_frag );
    af2 = get_frag( sa.get_connection(i).prev_frag );

    if( af1 != af2 ) { return false; } // should be pointing to the same frag
  }

  af1 = get_frag( get_last_connection().next_frag );
  af2 = get_frag( sa.get_last_connection().next_frag );

  if( af1 != af2 ) { return false; } // should be pointing to the same C-terminal cap

  return true; // if we reach this point, the Skel_Alis are equal
}


bool Skel_Ali::operator!=( const Skel_Ali& sa )
{
  return( !(*this == sa) );
}


list<int> Skel_Ali::get_sse_id_list()
{

  list<int> res;

  for( int i=0; i<(int)connections.size()-1; i++ ) {
    res.push_back( connections[i].next_frag.sse_idx );
  }

  return res;
}




// DEBUGGING
void Skel_Ali::print( string qseq, string tseq, int min_ali_res, ostream& os ) const
{
  os << "-----------" <<  endl;
  os << "Skel info:    " <<  endl;
  os << "#frags:       " << connections.size() << endl;
  os << "score:        " << score << endl;
  os << "native shift: " << shift << endl;
  os << "SSE_CO:       " << SSE_CO << endl;
  os << "cov_res:      " << num_aligned_residues << endl;
  os << "Frags:        " << endl;

  // print the N-terminus and its connection score
  os << endl;
  get_frag( connections[0].prev_frag )->print(qseq, tseq, os );
  os << "cnxn score: " << connections[0].connection_score << endl;
  os << endl;

  int beg;
  int end;

  // print each fragment and its connection score
  for( int i=1; i<(int)connections.size(); i++ ) {
    beg = connections[i-1].next_beg_res_idx;
    end = connections[i].prev_end_res_idx;
    os << endl;
    get_frag( connections[i].prev_frag )->print(qseq, tseq, beg, end, os );
    os << endl;
    os << "cnxn score: " << connections[i].connection_score << endl;
  }

  // print the C-terminus (or whatever the last fragment is)
  os << endl;
  beg = connections.back().next_beg_res_idx;
  end = get_frag( connections.back().next_frag)->core_t1();
  get_frag( connections.back().next_frag )->print(qseq, tseq, beg, end, os );
  os << endl;

  os << "-----------" <<  endl;
}

void Skel_Ali::print() const
{
  cerr << "Skel info:    " <<  endl;
  cerr << "#frags:       " << connections.size() << endl;
  cerr << "score:        " << score << endl;
  cerr << "native shift: " << shift << endl;
  cerr << "SSE_CO:       " << SSE_CO << endl;
  cerr << "cov_res:      " << num_aligned_residues << endl;
  cerr << "Frags:        " << endl;

  for( int i=0; i<(int)connections.size(); i++ ) {
    if( i > 0 ) {
      cerr << "Includes template residues: "
	   << connections[i-1].next_beg_res_idx << " to "
	   << connections[i].prev_end_res_idx << endl;
    }

    get_frag( connections[i].prev_frag )->print();
    cerr << endl;
  }

  cerr << "last frag above connected to (starting at " << connections.back().next_beg_res_idx << "): " << endl;
  get_frag( connections.back().next_frag )->print();
  cerr << endl;

}


string Skel_Ali::print2() const
{
  stringstream res;

  res << "Skel info:    " <<  endl;
  res << "#frags:       " << connections.size() << endl;
  res << "score:        " << score << endl;
  res << "native shift: " << shift << endl;
  res << "SSE_CO:       " << SSE_CO << endl;
  res << "cov_res:      " << num_aligned_residues << endl;
  res << "Frags:        " << endl;

  for( int i=0; i<(int)connections.size(); i++ ) {
    res << get_frag( connections[i].prev_frag )->print2();
    res << endl;
  }

  res << get_frag( connections.back().next_frag )->print2();

  return res.str();
}


void Skel_Ali::print_plain() const
{
  for( unsigned int i=0; i<connections.size(); i++ ) {
    cerr << "(" << connections[i].prev_frag.sse_idx << ","
	 << connections[i].prev_frag.frag_idx << ") - ";
  }

  cerr << "(" << connections.back().next_frag.sse_idx << "," 
       << connections.back().next_frag.frag_idx << ") - END" << endl;
}
