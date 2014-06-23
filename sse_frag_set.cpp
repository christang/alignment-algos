#include "sse_frag_set.h"

using namespace std;

/*********************************************
* implementation of SSE_Frag_Set class *
*********************************************/

// constructor 1
SSE_Frag_Set::SSE_Frag_Set( int id, int t_beg, int t_end, int qt_beg, int qt_end, 
			    int q_len, int t_len, int type, 
			    vector<Ali_Frag> unused, int min_cov_res )
{
  // loads all fragments into 'Frags', makes them all available, sets each frag's ...
  // frag_id to equal its idx in 'Frags'
  // used in Frag_Matrix::create_all_fragments for frags in template SSEs

  active_flag = 1;     // setup flag values (hard-coded) -- define elsewhere?
  available_flag = 0;
  redundant_flag = -1;

  Frags.reserve( qt_end - qt_beg + 1 );
  Frags = unused; // store 'unused' in 'Frags'

  sse_id = id; // setup SSE-defining variables
  t0 = t_beg;
  t1 = t_end;
  sse_len = ( t1 - t0 ) + 1;

  qt_shift_lo = qt_beg;
  qt_shift_hi = qt_end;

  query_len = q_len;
  templ_len = t_len;

  ss_type = type;

  min_aligned_residues = min_cov_res;

  for( int i=0; i<(int)Frags.size(); i++ ) {
    Frags[i].sse_id = sse_id;
    Frags[i].frag_id = i; // the frag_id in Ali_Frag matches the index in Frags
    Frags[i].make_available(); // make all fragments available initially
  }

}


SSE_Frag_Set::~SSE_Frag_Set()
{}


// OPERATORS

// GENERAL

void SSE_Frag_Set::find_biggest_gap( int& max_gap, int& qt_max_gap_beg, int& qt_max_gap_end )
{
  // sets the input variables to reflect the size of the biggest gap and its location ...
  // in terms of qt_shift's

  // copy frags to a new vector 'ordered_frags' in order by sse_q0
  vector<Ali_Frag*> ordered_frags = get_ordered_frags();

  // now that the frags are ordered by qt_shift, find the biggest gap
  max_gap = -1;
  int cur_gap;

  for( int i=1; i<(int)ordered_frags.size(); i++ ) {

    cur_gap = ( ordered_frags[i]->qt() - ordered_frags[i-1]->qt() ) - 1;

    if( cur_gap > max_gap ) {
      max_gap = cur_gap;
      qt_max_gap_beg = ordered_frags[i-1]->qt() + 1;
      qt_max_gap_end = ordered_frags[i]->qt() - 1;
    }
  }

  // last, check to see if the gaps between the first valid query residue 
  // and the first frag in the list or the last valid query residue and the
  // last frag in the list are bigger than the biggest inter-frag gap found above

  if( ( ordered_frags.front()->qt() - qt_shift_lo ) > max_gap ) {
    max_gap = ordered_frags.front()->qt() - qt_shift_lo;
    qt_max_gap_beg = qt_shift_lo;
    qt_max_gap_end = ordered_frags.front()->qt() - 1;
  }

  if( ( qt_shift_hi - ordered_frags.back()->qt() ) > max_gap ) {
    max_gap = qt_shift_hi - ordered_frags.back()->qt();
    qt_max_gap_beg = ordered_frags.back()->qt() + 1;
    qt_max_gap_end = qt_shift_hi;
  }

}


void SSE_Frag_Set::fill_gap( int gap_beg, int gap_end )
{
  // establish ranges for frag search
  int q_top_range, q_bot_range;

  if( ( gap_end - gap_beg + 1 ) > 5 ) {
    q_top_range = gap_end - (int)( (float)( gap_end - gap_beg ) / 3.0f );
    q_bot_range = gap_beg + (int)( (float)( gap_end - gap_beg ) / 3.0f );
  }
  else {
    q_top_range = gap_end;
    q_bot_range = gap_beg;
  }

  vector<Ali_Frag*> available_frags = get_available_frags();

  // search through 'available_frags' for the highest-scoring frag in the range
  // NOTE: it is assumed that 'available_frags' is in order with the highest ...
  // scoring frag first since 'Frags' is in that order

  vector<Ali_Frag*>::iterator it;

  for( it = available_frags.begin(); it != available_frags.end(); it++ ) {
    if( ( (*it)->qt() >= q_bot_range ) && ( (*it)->qt() <= q_top_range ) ) {
      break;
    }
  }

  // exit if something goes wrong and no frag exists in the range  
  if( it == available_frags.end() ) {
    cerr << "never found a fragment in the range: sse_id " << sse_id << " - "
	 << gap_beg << " to " << gap_end << endl;
    exit(-1);
  }

  // old way
  /*
  Active.push_back( *it ); // add the frag to Active
  Active.back().frag_id = (int)Active.size() - 1; // give the new frag the proper frag_id
  Available.erase( it );    // delete it from Available
  */

  // new way
  (*it)->make_active();

}


vector<Ali_Frag*> SSE_Frag_Set::get_ordered_frags()
{
  // order all active frags by their qt shift (in ascending order)

  vector<Ali_Frag*> res = get_active_frags();

  // bubble sort 'res' according to qt_shift
  for( int i=0; i<(int)res.size()-1; i++ ) {
    for( int j=i+1; j<(int)res.size(); j++ ) {

      if( res[i]->qt() > res[j]->qt() ) { // swap these
	Ali_Frag* tmp = res[i];
	res[i] = res[j];
	res[j] = tmp;
      }
    }
  }

  return res;
}


vector<Ali_Frag*> SSE_Frag_Set::get_all_frags_qt_sorted()
{
  vector<Ali_Frag*> res;

  /*
  for( int i=0; i<(int)Active.size(); i++ ) { res.push_back( &Active[i] ); }
  for( int i=0; i<(int)Available.size(); i++ ) { res.push_back( &Available[i] ); }
  for( int i=0; i<(int)Redundant.size(); i++ ) { res.push_back( &Redundant[i] ); }
  */
  for( int i=0; i<(int)Frags.size(); i++ ) { res.push_back( &Frags[i] ); }

  for( int i=0; i<(int)res.size()-1; i++ ) {
    for( int j=i+1; j<(int)res.size(); j++ ) {

      if( res[j]->qt() < res[i]->qt() ) { // swap
	Ali_Frag* tmp = res[j];
	res[j] = res[i];
	res[i] = tmp;
      }

    }
  }

  return res;
}


vector<Ali_Frag*> SSE_Frag_Set::get_active_frags()
{
  // returns a vector<Ali_Frag*> representing all active frags

  vector<Ali_Frag*> res;

  for( int i=0; i<(int)Frags.size(); i++ ) {
    if( Frags[i].is_active() ) {
      res.push_back( &Frags[i] );
    }
  }

  return res;
}


vector<Ali_Frag*> SSE_Frag_Set::get_available_frags()
{
  // returns a vector<Ali_Frag*> representing all active frags

  vector<Ali_Frag*> res;

  for( int i=0; i<(int)Frags.size(); i++ ) {
    if( Frags[i].is_available() ) {
      res.push_back( &Frags[i] );
    }
  }

  return res;
}


float SSE_Frag_Set::get_highest_available_frag_zscore()
{
  // returns the zscore of the highest-scoring fragment in this column
  // assumes that 'Frags' (and thus 'available_frags' is sorted with the ...
  // highest-scoring frag first

  vector<Ali_Frag*> available_frags = get_available_frags();

  return available_frags.front()->zs();
}


void SSE_Frag_Set::activate_top_available_frag()
{
  // since we assume 'Frags' to be sorted by score (highest first), activate ...
  // the first available frag

  for( int i=0; i<(int)Frags.size(); i++ ) {
    if( Frags[i].is_available() ) {
      activate_frag( Frags[i].frag_id );
      return;
    }
  }

  // error-checking, should never get here
  cerr << "Could not find an available frag. Exiting." << endl;
  exit(-1);

}


void SSE_Frag_Set::set_frag_zscores()
{
  //  vector<Ali_Frag> all_frags; // load all frags into a single vector

  /*
  for( int i=0; i<(int)Active.size(); i++ ) {
    all_frags.push_back( Active[i] );
  }

  for( int i=0; i<(int)Available.size(); i++ ) {
    all_frags.push_back( Available[i] );
  }

  for( int i=0; i<(int)Redundant.size(); i++ ) {
    all_frags.push_back( Redundant[i] );
  }
  */

  // calc the average score for frags in this column
  float sum(0.0f), average;

  for( int i=0; i<(int)Frags.size(); i++ ) {
    sum += Frags[i].ss();
  }

  average = sum / (float)Frags.size();

  // calc the standard deviation for frags in this column
  float standard_deviation;
  sum = 0;

  for( int i=0; i<(int)Frags.size(); i++ ) {
    sum += pow( ( Frags[i].ss() - average ), 2 );
  }

  standard_deviation = sqrt( ( 1 / (float)Frags.size() ) * sum );

  // set frag z-scores
  /*
  for( int i=0; i<(int)Active.size(); i++ ) {
    Active[i].z_score = ( Active[i].ss() - average ) / standard_deviation;
  }

  for( int i=0; i<(int)Available.size(); i++ ) {
    Available[i].z_score = ( Available[i].ss() - average ) / standard_deviation;
  }

  for( int i=0; i<(int)Redundant.size(); i++ ) {
    Redundant[i].z_score = ( Redundant[i].ss() - average ) / standard_deviation;
  }
  */
  for( int i=0; i<(int)Frags.size(); i++ ) {
    Frags[i].z_score = ( Frags[i].ss() - average ) / standard_deviation;
  }

}


void SSE_Frag_Set::activate_frag( int selected_frag_id )
{
  int width;

  if( ss_type == 329 ) { width = 2; } // helix
  else if (ss_type == 330 ) { width = 0; } // strand
  else { cerr << "Invalid SSE type in SSE " << sse_id << ".  Exiting." << endl; exit(-1); }

  // get the qt_shifts of all frags within 'width' qt shifts of selected_frag_qt (from Available)
  vector<int> Neighbors = find_available_neighbors( selected_frag_id, width );

  Frags[ selected_frag_id ].make_active();

  // make neighboring frags redundant
  for( int i=0; i<(int)Neighbors.size(); i++ ) {
    Frags[ Neighbors[i] ].make_redundant();
  }

}


vector<int> SSE_Frag_Set::find_available_neighbors( int selected_frag_id, int neighborhood_width )
{
  // return a vector of ints representing the frag_id's of fragments within 'neigborhood_width' ...
  // of the qt_shift of selected_frag_id

  int center_frag_qt_shift = get_frag( selected_frag_id )->qt();

  vector<int> res;

  for( int i=0; i<(int)Frags.size(); i++ ) {

    if( Frags[i].is_available() ) { // only interested in available neighbors

      int qt_separation = abs( Frags[i].qt() - center_frag_qt_shift );

      if( ( qt_separation <= neighborhood_width ) && ( qt_separation != 0 ) ) {
	res.push_back( Frags[i].frag_id );
      }

    }

  }

  return res;
}


bool SSE_Frag_Set::an_available_frag_exists() const
{
  // determine whether an available frag remains in this SSE_Frag_Set

  for( int i=0; i<(int)Frags.size(); i++ ) {
    if( Frags[i].is_available() ) { return true; }
  }

  return false;
}


vector<Ali_Frag*> SSE_Frag_Set::find_shift_neighbors( float qt_target, int num_frags )
{
  // returns a vector of Ali_Frag*'s (of length 'num_frags') containing the active ...
  // fragments that are closest in terms of qt_shift to 'qt_target'

  vector<Ali_Frag*> res;

  // fill res with pointers to all active fragments
  for( int i=0; i<(int)Frags.size(); i++ ) {
    if( Frags[i].is_active() ) { res.push_back( &Frags[i] ); }
  }

  // sort res according to each frag's absolute qt_shift distance from qt_target
  for( int i=0; i<(int)res.size() - 1; i++ ) {
    for( int j=i+1; j<(int)res.size(); j++ ) {

      float qt_diff_i = abs( (float)res[i]->qt() - qt_target );
      float qt_diff_j = abs( (float)res[j]->qt() - qt_target );

      if( qt_diff_j < qt_diff_i ) { // swap i and j
	Ali_Frag* tmp = res[i];
	res[i] = res[j];
	res[j] = tmp;
      }
    }
  }

  // shorten res to num_frags elements
  while( (int)res.size() > num_frags ) {
    res.pop_back();
  }

  return res;
}


int SSE_Frag_Set::get_frag_status( Ali_Frag* af )
{
  return( af->get_status() );
}


int SSE_Frag_Set::get_num_active_frags()
{
  // return the number of active frags in this SSE_Frag_Set

  vector<Ali_Frag*> tmp = get_active_frags();
  return tmp.size();
}


// PRINTING
void SSE_Frag_Set::print_frag_scores( ofstream& out ) const
{

  out << "ACTIVE" << endl;

  for( int i=0; i<(int)Frags.size(); i++ ) {
    if( Frags[i].is_active() ) {
      out << Frags[i].frag_id << "\t";
      out << Frags[i].core_len() << "\t";
      out << Frags[i].qt() << "\t";
      out << Frags[i].ss() << "\t";
      out << Frags[i].ss() / (float)Frags[i].core_len() << "\t";
      out << Frags[i].z_score << endl;
    }
  }

  out << "AVAILABLE" << endl;

  for( int i=0; i<(int)Frags.size(); i++ ) {
    if( Frags[i].is_available() ) {
      out << Frags[i].frag_id << "\t";
      out << Frags[i].core_len() << "\t";
      out << Frags[i].qt() << "\t";
      out << Frags[i].ss() << "\t";
      out << Frags[i].ss() / (float)Frags[i].core_len() << "\t";
      out << Frags[i].z_score << endl;
    }
  }

  out << "REDUNDANT" << endl;

  for( int i=0; i<(int)Frags.size(); i++ ) {
    if( Frags[i].is_redundant() ) {
      out << Frags[i].frag_id << "\t";
      out << Frags[i].core_len() << "\t";
      out << Frags[i].qt() << "\t";
      out << Frags[i].ss() << "\t";
      out << Frags[i].ss() / (float)Frags[i].core_len() << "\t";
      out << Frags[i].z_score << endl;
    }
  }

  out << endl;

}

void SSE_Frag_Set::print_frag_scores( ostream& out ) const
{

  out << "ACTIVE" << endl;

  for( int i=0; i<(int)Frags.size(); i++ ) {
    if( Frags[i].is_active() ) {
      out << Frags[i].frag_id << "\t";
      out << Frags[i].core_len() << "\t";
      out << Frags[i].qt() << "\t";
      out << Frags[i].ss() << "\t";
      out << Frags[i].ss() / (float)Frags[i].core_len() << "\t";
      out << Frags[i].z_score << endl;
    }
  }

  out << "AVAILABLE" << endl;

  for( int i=0; i<(int)Frags.size(); i++ ) {
    if( Frags[i].is_available() ) {
      out << Frags[i].frag_id << "\t";
      out << Frags[i].core_len() << "\t";
      out << Frags[i].qt() << "\t";
      out << Frags[i].ss() << "\t";
      out << Frags[i].ss() / (float)Frags[i].core_len() << "\t";
      out << Frags[i].z_score << endl;
    }
  }

  out << "REDUNDANT" << endl;

  for( int i=0; i<(int)Frags.size(); i++ ) {
    if( Frags[i].is_redundant() ) {
      out << Frags[i].frag_id << "\t";
      out << Frags[i].core_len() << "\t";
      out << Frags[i].qt() << "\t";
      out << Frags[i].ss() << "\t";
      out << Frags[i].ss() / (float)Frags[i].core_len() << "\t";
      out << Frags[i].z_score << endl;
    }
  }

  out << endl;

}


void SSE_Frag_Set::print_sse_info( ostream& out ) const
{
  out << "SSE id: " << sse_id << endl;

  out << "Type: ";
  if( ss_type == 329 ) { out << "Helix" << endl; }
  else if( ss_type == 330 ) { out << "Strand" << endl; }
  else { out << "Undefined" << endl; }

  out << "T: " << t0 << " - " << t1 << endl;
  out << "QT: " << qt_shift_lo << " - " << qt_shift_hi << endl;
}

void SSE_Frag_Set::print_sse_info( string templ_seq, ostream& out ) const
{
  print_sse_info( out );

  out << "Seq: ";

  for( int t=t0; t<=t1; t++ ) {
    out << templ_seq[t];
  }

  out << endl;
}

// DEBUGGING
void SSE_Frag_Set::get_statistics( int& num_frags, int& ali_len, float& avg_gap, int& max_gap )
{
  int tmp1, tmp2; // throw away variables

  num_frags = get_num_active_frags();
  ali_len = ( query_len - sse_len ) - 1;
  avg_gap = (float)( ali_len - num_frags ) / (float)( num_frags + 1 );
  find_biggest_gap( max_gap, tmp1, tmp2 );
}
