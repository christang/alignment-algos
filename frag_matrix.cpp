#include "frag_matrix.h"

using namespace std;

/*********************************************
* implementation of Frag_Matrix class *
*********************************************/

/********************************************/
/************* PRIVATE FUNCTIONS ************/
/********************************************/

// connection
bool Frag_Matrix::connection_is_valid( Frag_Set* fs, Ali_Frag* af1, Ali_Frag* af2 )
{

  bool res = false;
  int t1_prev, q1_prev, t0_next, q0_next;

  get_inter_frag_endpoints( af1, af2, t1_prev, q1_prev, t0_next, q0_next );
  
  if( fs->frags_in_order( t1_prev, q1_prev, t0_next, q0_next ) ) {
    
    if( Str_Data->get_TSR_to_N( t1_prev ) + Str_Data->get_TSR_to_C( t0_next ) > min_aligned_residues ) {
      
      if( loop_spans_gap( t1_prev, q1_prev, t0_next, q0_next ) ) {
	
	// this connection passed all the tests, return true
	res = true;
      }
      else {}
    }
    else {}
  }
  else {}

  return res;
}


bool Frag_Matrix::loop_spans_gap( int t1_prev, int q1_prev, int t0_next, int q0_next )
{
  // number of residues between ends of SSEs = loop_res_len = ( q0_next - q1_prev ) - 1;
  // we are using loop_res_len + 1 here

  return( Str_Data->get_cb_dist( t1_prev, t0_next ) < (float)( q0_next - q1_prev ) * 3.3f );
}


void Frag_Matrix::get_connection_info( Frag_Set* fs,
				       Frag_ID prev_frag_id, Frag_ID next_frag_id, 
				       int& prev_end, int& next_beg,
				       float& connect_score )
{
  // can extend connected fragments toward each other, ...
  // starting at their cores and extending to the end ...
  // of the SSE.  stores the highest-scoring extension ... 
  // in 'prev_end' and 'next_beg' and the corresponding ...
  // score in 'connect_score.

  // case 0: simply use the core values, don't look for extensions.
  // case 1: as described above

  Ali_Frag* prev_frag = fs->get_frag( prev_frag_id );
  Ali_Frag* next_frag = fs->get_frag( next_frag_id );

  switch( ali_mode ) {
    
  case 0:
    {
      prev_end = prev_frag->core_t1();
      next_beg = next_frag->core_t0();
      connect_score = 0.f;
      break;
    }
  case 1:
    {
      int max_prev_end(-1);
      int max_next_beg(-1);
      float max_connect_score(-1000.f);

      for( int t_prev=prev_frag->core_t1();
	   t_prev<=prev_frag->sse_t1(); t_prev++ ) {
	// search forward from last core res to last sse res

	for( int t_next=next_frag->core_t0();
	     t_next>=next_frag->sse_t0(); t_next-- ) {
	  // search backward from first core res to first sse res

	  if( prev_frag->frag_is_N_terminal ||
	      next_frag->frag_is_C_terminal ||
	      ( fs->frags_in_order( t_prev, prev_frag->q( t_prev ),
				    t_next, next_frag->q( t_next ) ) &&
		loop_spans_gap( t_prev, prev_frag->q( t_prev ),
				t_next, next_frag->q( t_next ) ) ) ) {
	  
	    // must match conditions: 
	    // 1) always allow connections to terminii
	    // 2) frags must be in order (next can follow prev)
	    // 3) ther are enough loop residues between next and prev 

	    float curr_score(0.f);

	    // compute the prev side of the connection score
	    for( int t1=(prev_frag->core_t1()+1); t1<=t_prev; t1++ ) {
	      curr_score += Str_Data->get_sim( prev_frag->q(t1), t1 );
	    }

	    // compute the next side of the connection score
	    for( int t2=(next_frag->core_t0()-1); t2>=t_next; t2-- ) {
	      curr_score += Str_Data->get_sim( next_frag->q(t2), t2 );
	    }

	    if( curr_score > max_connect_score ) {
	      max_connect_score = curr_score;
	      max_prev_end = t_prev;
	      max_next_beg = t_next;
	    }
	 
	  }
 
	}

      }

      prev_end = max_prev_end;
      next_beg = max_next_beg;
      connect_score = max_connect_score;

      break;
    }
  default:
    cerr << "Improper choice of mode.  Exiting." << endl;
    exit(-1);
  }

}


void Frag_Matrix::get_inter_frag_endpoints( Ali_Frag* prev_frag, Ali_Frag* next_frag,
					    int& t1_prev, int &q1_prev, int& t0_next, int& q0_next )
{
  t1_prev = prev_frag->core_t1();
  q1_prev = prev_frag->core_q1();
  t0_next = next_frag->core_t0();
  q0_next = next_frag->core_q0();
}


// filling the matrix
void Frag_Matrix::find_biggest_gap_in_Matrix( Frag_Set* fs, int& sse_id, int& q_gap_beg, int& q_gap_end )
{
  int max_gap = -1;
  int cur_gap, cur_q_gap_beg, cur_q_gap_end;

  for( int i=1; i<=num_sses; i++ ) {
    fs->get_col(i)->find_biggest_gap( cur_gap, cur_q_gap_beg, cur_q_gap_end );

    if( cur_gap >= max_gap ) {
      max_gap = cur_gap;
      sse_id = i;
      q_gap_beg = cur_q_gap_beg;
      q_gap_end = cur_q_gap_end;
    }
  }

}


// general
bool Frag_Matrix::it_is_valid_starting_frag( Frag_Set* fs, Frag_ID f, int t_next_beg )
{
  Ali_Frag* af = get_frag( fs, f );

  if( af->frag_is_C_terminal ) { return false; }
  else {
    return( ( af->core_q0() < ( ( query_len - 2 ) - min_aligned_residues ) ) &&
	    ( af->core_t0() < ( ( templ_len - 2 ) - min_aligned_residues ) ) &&
	    ( Str_Data->get_TSR_to_C( t_next_beg ) > min_aligned_residues ) );

  }
}


bool Frag_Matrix::sse_is_native( int t_beg, int t_end ) const
{
  if( Compare_to_Native == NULL ) { return false; } // cannot do this when not in tracking mode

  vector<Res_Pair> sse_ali = Compare_to_Native->get_local_native_ali( t_beg, t_end );

  return( (int)sse_ali.size() >= find_min_ali_len( t_end - t_beg + 1 ) );
}


int Frag_Matrix::find_min_ali_len( int sse_len ) const
{
  if( sse_len <= 4 ) { return sse_len; }
  else if( sse_len <= 6 ) { return 5; }
  else if( sse_len <= 8 ) { return 6; }
  else if( sse_len <= 10 ) { return 7; }
  else if( sse_len <= 14 ) { return 9; }
  else if( sse_len <= 20 ) { return 11; }
  else if( sse_len <= 30 ) { return 15; }
  else { return 20; }
}


/********************************************/
/************* PUBLIC FUNCTIONS *************/
/********************************************/

// constructor 1
Frag_Matrix::Frag_Matrix( int min_cov_res,
			  Frag_Set* fs, Ali_Str_Info* asi,
			  int max_frag_shift, int m, Ali_Dist* ad )
{

  Main_FS = fs;
  Str_Data = asi;

  templ_seq = Str_Data->get_templ_seq();
  query_seq = Str_Data->get_query_seq();

  templ_len = Str_Data->get_templ_len();
  query_len = Str_Data->get_query_len();

  num_sses = Str_Data->get_num_templ_SSEs();

  max_in_betw_shift = max_frag_shift;
  ali_mode = m;
  Compare_to_Native = ad;

  debug_off();

  min_aligned_residues = min_cov_res;
}


// destructor
Frag_Matrix::~Frag_Matrix()
{}


// setting up fragment set
void Frag_Matrix::create_all_fragments( Frag_Set* fs )
{
  cerr << "create_all_fragments: top" << endl;

  fs->clear_all();

  int sse_id=0;
  vector<Ali_Frag> top_sse_frags;
  vector<Ali_Frag> bottom_sse_frags;

  top_sse_frags.push_back( Ali_Frag( 0, 0, 0, 0.f, true, false ) ); // setup the N-terminal frag

  fs->add_column(  SSE_Frag_Set( sse_id, 0, 0, -1, -1, 
				 query_len, templ_len, -1,
				 top_sse_frags, min_aligned_residues ) );

  for( sse_id=1; sse_id<=num_sses; sse_id++ ) {

    list<Ali_Frag> current_sse_frags;

    int t0 = Str_Data->get_SSE_Data(sse_id-1).beg_id; // start of template SSE
    int t1 = Str_Data->get_SSE_Data(sse_id-1).end_id; // end of template SSE

    int SSE_len = t1 - t0 + 1;

    int min_ali_len = find_min_ali_len( SSE_len );

    int q0_lo = max( min_ali_len - SSE_len + 1, t0 + min_aligned_residues - ( templ_len - 2 ) );
    int q0_hi = min( ( query_len - 2 ) - min_ali_len + 1, t0 - min_aligned_residues + ( query_len - 2 ) );

    int qt_shift_lo = q0_lo - t0;
    int qt_shift_hi = q0_hi - t0;

    for( int q0=q0_lo; q0<=q0_hi; q0++ ) {

      int qt_shift = q0 - t0; // query_id - templ_id

      // NOTE: shouldn't need this anymore!
      // skip fragments in the extreme lower-left and upper-right corners
      // they can't participate in alignments that surpass the min coverage
      if( ( qt_shift > ( query_len - 2 ) - min_aligned_residues ) ||
	  ( qt_shift < min_aligned_residues - ( templ_len - 2 ) ) ) {
	continue;
      }
    
      float score, max_score(-1000.f);
      int max_i(-1);

      for( int i=0; i<(SSE_len-min_ali_len+1); i++ ) {

	// check if frag is starting (or ending) before the first query residue ...
	// (or after the last query residue)
	if( ( ( q0+i ) < 1 ) || ( ( q0+i+(min_ali_len-1) ) > ( query_len - 2 ) ) ) {
	  continue;
	}

	bool frag_does_not_contain_loop( true );
	score = 0;

	for( int j=0; j<min_ali_len; j++ ) {

	  score += Str_Data->get_sim( q0+i+j, t0+i+j );

	  frag_does_not_contain_loop = ( frag_does_not_contain_loop &&
					 !Str_Data->get_query_predicted_loop( q0+i+j ) );

	}

	if( ! frag_does_not_contain_loop ) { 
	  //	  cerr << "FRAG ELIMINATED DUE TO PREDICTED LOOP" << endl;
	  //	  pause();
	  //	  continue;
	} // this frag overlaps a well-predicted loop

	if( score > max_score ) {
	  max_score = score;
	  max_i = i;
	}

      }

      if( max_score == -1000.f ) { continue; } // never found a fragment, don't make an Ali_Frag

      Ali_Frag tmp( max( 1, t0+qt_shift ) - qt_shift,
		    min( (query_len-2), t1+qt_shift ) - qt_shift,
		    (t0+max_i), (t0+max_i)+(min_ali_len-1), qt_shift,
		    max_score, false, false );

      current_sse_frags.push_back( tmp );
    }

    current_sse_frags.sort();
    current_sse_frags.reverse();

    top_sse_frags.clear();

    // copy ALL frags to top_sse_frags
    while( !current_sse_frags.empty() ) {
      top_sse_frags.push_back( current_sse_frags.front() );
      current_sse_frags.pop_front();
    }

    bottom_sse_frags.clear();

    // store ALL frags in Available set of this column, not just top N
    fs->add_column( SSE_Frag_Set( sse_id, t0, t1, qt_shift_lo, qt_shift_hi,
				  query_len, templ_len, Str_Data->get_SSE_Data(sse_id-1).ss_type,
				  top_sse_frags, min_aligned_residues ) );

  }

  // the (n+1)th column holds only the (imaginary) C-terminal frag, which is just a marker for (templ_len,query_len)

  // setup the C-terminal frag
  top_sse_frags.clear();
  bottom_sse_frags.clear();
  top_sse_frags.push_back( Ali_Frag( templ_len-1, templ_len-1,
				     (query_len-1) - (templ_len-1),
				     0.f,
				     false, true ) );

  fs->add_column( SSE_Frag_Set( sse_id, templ_len-1, templ_len-1, -1, -1,
				query_len, templ_len, -1,
				top_sse_frags, min_aligned_residues ) );

  fs->activate_terminal_caps();

  cerr << "create_all_fragments: end" << endl;
}


void Frag_Matrix::find_fragment_connections( Frag_Set* fs )
{
  // This function does a lot of the work between the 'finding fragments phase' and the
  // search for high-scoring combinations of fragments phase.  For each fragment, it finds
  // which other fragments it can link to.

  // for each Ali_Frag x, find every Ali_Frag y that could follow x in a valid alignment.
  for( int i=1; i<=num_sses; i++ ) {

    for( int j=0; j<fs->num_frags_in_sse( i ); j++ ) { 

      get_frag( fs, i, j )->clear_next();  // erase all previous connections

      int prev_end, next_beg;
      float connect_score;

      // find all frags in subsequent sses that follow Columns[i].Active[j]
      for( int m=i+1; m<=num_sses; m++ ) {

	for( int n=0; n<fs->num_frags_in_sse( m ); n++ ) {

	  // check frags' relative order, min coverage and loop constraint
	  if( connection_is_valid( fs, get_frag(fs,i,j), get_frag(fs,m,n) ) ) {

	    // passed test, so make the connection
	    get_connection_info( fs, get_frag(fs,i,j )->get_id(), get_frag(fs,m,n)->get_id(),
				 prev_end, next_beg, connect_score );
	    
	    get_frag( fs, i, j)->make_connection( get_frag( fs, m, n)->get_id(),
						  prev_end, next_beg, connect_score );
	  }

	}
      }

      // every frag must connect to the C-terminus
      get_connection_info( fs, get_frag( fs, i, j )->get_id(), get_frag( fs, num_sses+1, 0 )->get_id(),
			   prev_end, next_beg, connect_score );

      get_frag( fs, i, j)->make_connection( get_frag( fs, num_sses+1, 0)->get_id(),
					    prev_end, next_beg, connect_score );

    }
  }

}


unsigned long long int Frag_Matrix::get_number_of_alis_to_search( Frag_Set* fs )
{
  find_N_terminal_connections( fs, 0 );

  return( get_frag( fs, 0, 0 )->get_num_children() );
}


void Frag_Matrix::find_N_terminal_connections( Frag_Set* fs, int mode )
{

  unsigned long long int num_children = 0;

  get_frag( fs, 0, 0 )->clear_next();

  for( int m=1; m<=num_sses; m++ ) {
    for( int n=0; n<fs->num_frags_in_sse( m ); n++ ) {

      Frag_ID curr_frag = get_frag( fs, m, n )->get_id();

      int prev_end, next_beg;
      float connect_score;
      
      get_connection_info( fs, get_frag( fs, 0, 0 )->get_id(), curr_frag,
			   prev_end, next_beg, connect_score );

      if( it_is_valid_starting_frag( fs, curr_frag, next_beg ) ) {

	get_frag( fs, 0, 0 )->make_connection( curr_frag, prev_end, next_beg, connect_score );

	num_children += get_frag( fs, m, n )->get_num_children();
      }

    }
  }

  get_frag( fs, 0, 0 )->set_num_children( num_children );
}


// filling the matrix
float Frag_Matrix::fill_Frag_Set_by_zscore( Frag_Set* fs )
{
  float zsc = fs->activate_next_best_available_frag();

  find_fragment_connections( fs ); // re-do all fragment connections now that we have a new frag

  fs->count_frag_children();  // re-calculate the number of children of each fragment

  return zsc;
}


void Frag_Matrix::fill_gaps_in_Frag_Matrix( Frag_Set* fs )
{
  int sse_id(-1), q_gap_beg(-1), q_gap_end(-1);

  find_biggest_gap_in_Matrix( fs, sse_id, q_gap_beg, q_gap_end );

  fs->get_col(sse_id)->fill_gap( q_gap_beg, q_gap_end );

  find_fragment_connections( fs ); // re-do all fragment connections now that we have a new frag

  fs->count_frag_children();  // re-calculate the number of children of each fragment
}


bool Frag_Matrix::activate_next_fragment( float& next_z, unsigned long long int& max_search, Frag_Set* fs )
{

  unsigned long long int num_alis_to_search = get_number_of_alis_to_search( fs );
  
  if( num_alis_to_search >= max_search ) { return false; }

  cerr << "Search space: " << num_alis_to_search << "\t";

  if( fs->an_available_frag_exists() ) 
    {
      next_z = fill_Frag_Set_by_zscore( fs );
    }
  else
    {
      cerr << endl;
      return false;
    }
  
  cerr << "New frag z-score: " << next_z << endl;
  
  return true;
}

/*
void Frag_Matrix::add_in_between_frags( Frag_Set* fs, float zsc_thresh )
{
  int max_shift_diff = max_in_betw_shift;

  // first, find all pairs of Ali_Frags that are within max_shift_diff from each other ...
  // and are in order, store them in prev_frags and next_frags
  vector<Frag_ID> prev_frags, next_frags;

  Frag_ID prev_frag;
  Frag_ID next_frag;

  for( int i=1; i<=num_sses-2; i++ ) {
    for( int j=0; j<(int)fs->get_col(i)->Active.size(); j++ ) {

      for( int m=i+2; m<=num_sses; m++ ) {
	for( int n=0; n<(int)fs->get_col(m)->Active.size(); n++ ) {

	  prev_frag = make_Frag_ID(i,j);
	  next_frag = make_Frag_ID(m,n);

	  if( abs( get_frag( fs, prev_frag )->qt_shift - get_frag( fs, next_frag )->qt_shift ) <=
	      max_shift_diff ) {

	    prev_frags.push_back( prev_frag );
	    next_frags.push_back( next_frag );

	  }

	}
      }

    }
  }
  
  // next, examine the Available frags between each pair in prev_frags and next_frags ...
  // to be added, an in-between frag must:
  // 1) have a qt_shift within max_shift_diff of both prev_frag and next_frag
  // 2) be able to connect to both prev_frag and next_frag
  // 3) have a zcore such that the average zscore of prev_frag, next_frag and the new frag ...
  //    is greater than or equal to zsc_thresh

  // if several fragments in an SSE pass all these criteria, only activate the highest scoring one

  if( prev_frags.size() != next_frags.size() ) { // error checking
    cerr << "prev_frags and next_frags have different sizes.  Exiting." << endl;
    exit(-1);
  }

  for( int k=0; k<(int)prev_frags.size(); k++ )
    {
      prev_frag = prev_frags[k];
      next_frag = next_frags[k];

      for( int i=prev_frag.sse_idx+1; i<next_frag.sse_idx; i++ ) { // search in-between SSEs

	vector<Ali_Frag*> poss_frags; // possible frags that pass the above criteria

	for( int j=0; j<(int)fs->get_col(i)->Available.size(); j++ ) {
	  
	  Ali_Frag* curr_frag = &fs->Frag_Columns[i].Available[j];
	  
	  if( ( abs( curr_frag->qt_shift - get_frag( fs, prev_frag )->qt_shift ) > max_shift_diff ) ||
	      ( abs( curr_frag->qt_shift - get_frag( fs, next_frag )->qt_shift ) > max_shift_diff ) )
	    { // frag is not within range qt_shift-wise
	      continue;
	    }

	  if( !connection_is_valid( get_frag( fs, prev_frag ), curr_frag ) ||
	      !connection_is_valid( curr_frag, get_frag( fs, next_frag) ) )
	    {
	      continue;
	    }
	  
	  if( ( ( get_frag( fs, prev_frag )->z_score + 
		  curr_frag->z_score + 
		  get_frag( fs, next_frag )->z_score ) / 3.f ) < zsc_thresh )
	    {
	      continue; 
	    }

	  poss_frags.push_back( curr_frag );
	}

	// poss_frags now holds pointers to frags in Available that pass all the criteria
	// also, all frags in poss_frags are in the same SSE (Column[i]) 

	// since Available is ranked by score, the highest-scoring frag in poss_frags is on top
	// NOTE: if the above is true, do I need to keep looking once I've found one poss_frag?

	if( !poss_frags.empty() ) {

	  fs->Frag_Columns[i].activate_frag( poss_frags.front()->qt_shift );
	  cerr << "activated in-between frag!!!" << endl;
	}
 
      }
    }

}
*/

void Frag_Matrix::add_in_between_frags2( Frag_Set* main_fs )
{

  // first, find all pairs of Ali_Frags that are within max_shift_diff from each other, ...
  // are in order, have z-scores above 2, and an average z-score above 3
  // store them in prev_frags and next_frags

  vector<Frag_ID> prev_frags, next_frags;

  get_in_between_frag_pairs( main_fs, prev_frags, next_frags );

  exit(-1);

  for( int i=0; i<(int)prev_frags.size(); i++ ) {
    
    // focus on one pair
    if( ( prev_frags[i].sse_idx == 6 ) && ( next_frags[i].sse_idx == 10 ) ) {
      
      Frag_Set tmp_fs = *Main_FS; // create a temporary copy of the main Frag_Set

      // print the bounding frags
      cerr << "-----------------" << endl;
      cerr << "F1:" << endl;
      tmp_fs.get_frag( prev_frags[i] )->print( query_seq, templ_seq );
      cerr << endl;
      cerr << "F2:" << endl;
      tmp_fs.get_frag( next_frags[i] )->print( query_seq, templ_seq );
      cerr << endl;
      
      // determine constraints of region (defines a box in which to look for new fragments)
      int sse_left = prev_frags[i].sse_idx; // sse boundaries
      int sse_right = next_frags[i].sse_idx;
      int q_top = tmp_fs.get_frag( prev_frags[i] )->core_q1(); // query residue boundaries
      int q_bottom = tmp_fs.get_frag( next_frags[i] )->core_q0();
    
      // get active and available frags in R
      vector<Ali_Frag*> active_frags = tmp_fs.get_active_frags_in_region( sse_left, sse_right, q_top, q_bottom );
      vector<Ali_Frag*> available_frags = tmp_fs.get_available_frags_in_region( sse_left, sse_right, q_top, q_bottom );
      

      cerr << "Constraints on R:" << endl;
      cerr << "Between SSEs: " << sse_left << " and " << sse_right << endl;
      cerr << "Between query res: " << q_top << " and " << q_bottom << endl;
      cerr << "# active frags in R: " << active_frags.size() << endl;
      cerr << "# available frags in R: " << available_frags.size() << endl;
      cerr << endl;

      cerr << "active frags in R: " << endl;
      for( int j=0; j<(int)active_frags.size(); j++ ) {
	active_frags[j]->print( query_seq, templ_seq );
	cerr << endl;
      }
      cerr << endl << endl;
      cerr << "available frags in R: " << endl;
      for( int j=0; j<(int)available_frags.size(); j++ ) {
	available_frags[j]->print( query_seq, templ_seq );
	cerr << endl;
      }
      
    }
    
  }
  
  exit(-1);  
  
}


void Frag_Matrix::get_in_between_frag_pairs( Frag_Set* fs, 
					     vector<Frag_ID>& prev_frags,
					     vector<Frag_ID>& next_frags )
{
  int max_shift_diff = 3; // arbitrary!
  int max_sse_separation = 5; // arbitrary!

  cerr << "max_shift_diff: " << max_shift_diff << endl;
  pause();

  prev_frags.clear();
  next_frags.clear();

  Frag_ID prev_frag_id;
  Frag_ID next_frag_id;

  Ali_Frag *prev_frag, *next_frag;
  
  for( int i=1; i<=num_sses-2; i++ ) {
    
    vector<Ali_Frag*> active_frags_i = fs->get_active_frags(i);
    
    for( int j=0; j < (int)active_frags_i.size(); j++ ) {
      
      prev_frag = active_frags_i[j];
      prev_frag_id = prev_frag->get_id();

      for( int m=i+2; (m<=i+max_sse_separation && m<=num_sses); m++ ) {
	  
	vector<Ali_Frag*> active_frags_m = fs->get_active_frags(m);
	
	for( int n=0; n < (int)active_frags_m.size(); n++ ) {
	  
	  next_frag = active_frags_m[n];
	  next_frag_id = next_frag->get_id();

	  if( fs->frags_in_order( prev_frag, next_frag ) ) {

	    if( abs( prev_frag->qt() - next_frag->qt() ) <= max_shift_diff ) {
	      
	      if( ( ( prev_frag->zs() > 2.0f ) &&
		    ( next_frag->zs() > 2.0f ) ) &&
		  ( ( ( prev_frag->zs() + next_frag->zs() ) / 2.f ) > 3.0f ) ) {
		  
		prev_frags.push_back( prev_frag_id );
		next_frags.push_back( next_frag_id );
		
		cerr << "Will look for frags between: " << endl;
		prev_frag->print( query_seq, templ_seq );
		cerr << endl;
		next_frag->print( query_seq, templ_seq );
		cerr << endl << endl;
		pause();		  

	      }
	    }
	  }
	  
	}
      }
      
    }
  }
  
  cerr << "size of prev_frags: " << prev_frags.size() << endl;
  cerr << "size of next_frags: " << next_frags.size() << endl;

}


// access
Ali_Frag* Frag_Matrix::get_frag( Frag_Set* fs, Frag_ID f )
{
  return( fs->get_frag( f ) );
}


Ali_Frag* Frag_Matrix::get_frag( Frag_Set* fs, int sse_id, int frag_id )
{
  return( fs->get_frag( sse_id, frag_id ) );
}


Frag_ID Frag_Matrix::make_Frag_ID( int sse_id, int frag_id )
{
  Frag_ID res;
  res.sse_idx = sse_id;
  res.frag_idx = frag_id;
  return res;
}


// debugging
void Frag_Matrix::report_frag_quality( Frag_Set* fs )
{
  if( Compare_to_Native == NULL ) { return; } // cannot do this when not in tracking mode

  for( int i=1; i<=num_sses; i++ ) {

    cerr << "------SSE INFO----------" << endl;

    fs->get_col(i)->print_sse_info( templ_seq );

    int t_beg = fs->get_col(i)->t0;
    int t_end = fs->get_col(i)->t1;

    if( sse_is_native( t_beg, t_end ) ) {
      cerr << "NATIVE" << endl;

      float local_qt_shift = Compare_to_Native->get_local_qt_shift( t_beg, t_end );

      cerr << "Native shift: " << local_qt_shift << endl;
      cerr << "# Active frags:" << fs->get_col(i)->get_num_active_frags() << endl;
      cerr << "Top 5 (or less) closest frags:" << endl;
      cerr << "QT-shift (distance to native): " << endl;

      vector<Ali_Frag*> closest_frags = fs->get_col(i)->find_shift_neighbors( local_qt_shift, 5 );

      for( int i=0; i<(int)closest_frags.size(); i++ ) {
	Ali_Frag* tmp = closest_frags[i];
	cerr << tmp->qt() << "(" << abs( (float)tmp->qt() - local_qt_shift ) << ")\t";
      }
      cerr << endl;

    }
    else {
      cerr << "Not native." << endl;
    }

    cerr << endl;

    cerr << "------SSE INFO----------" << endl;
  }

}


void Frag_Matrix::report_full_sse_frag_set_info( Frag_Set* fs )
{
  if( Compare_to_Native == NULL ) { return; } // cannot do this when not in tracking mode

  for( int i=1; i<=num_sses; i++ ) {

    cerr << "------SSE FRAG SET----------" << endl;

    fs->get_col(i)->print_sse_info( templ_seq );

    int t_beg = fs->get_col(i)->t0;
    int t_end = fs->get_col(i)->t1;

    if( sse_is_native( t_beg, t_end ) ) {
      cerr << "NATIVE" << endl;

      float local_qt_shift = Compare_to_Native->get_local_qt_shift( t_beg, t_end );

      vector<Ali_Frag*> sse_frags = fs->get_col(i)->get_all_frags_qt_sorted();

      for( int j=0; j<(int)sse_frags.size(); j++ ) {

	sse_frags[j]->print_one_line( templ_seq, query_seq, cerr );
	cerr << ", " << (float)sse_frags[j]->qt() - local_qt_shift;

	if( fs->get_col(i)->get_frag_status( sse_frags[j] ) == 1 ) {
	  cerr << " -- ACTIVE ";
	}

	if( fs->get_col(i)->get_frag_status( sse_frags[j] ) == -1 ) {
	  cerr << " -- REDUNDANT";
	}

	if( fs->get_col(i)->get_frag_status( sse_frags[j] ) == -2 ) {
	  cerr << "Frag status undefined.  Frag not found in sse_frag_set. Exiting." << endl;
	  exit(-1);
	}

	cerr << endl;

      }

    }

    cerr << "------SSE FRAG SET----------" << endl;
  }

}


void Frag_Matrix::pause()
{
  cerr << "Paused. Enter a character and press Enter to continue." << endl;
  string pause;
  cin >> pause;
}
