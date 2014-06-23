#include "frag_set.h"

using namespace std;

/*********************************************
* implementation of Frag_Set class *
*********************************************/

// constructor 1
Frag_Set::Frag_Set()
{
  set_verbose( false );
}

Frag_Set::~Frag_Set()
{
}

// OPERATORS
vector<Ali_Frag*> Frag_Set::operator-( Frag_Set& fs )
{
  // return a list of fragments that are active in 'this', but not in 'fs'

  vector<Ali_Frag*> res;

  vector<Ali_Frag*> active_frags;

  // load 'active_frags' with all active frags from all columns
  for( int i=1; i<=num_sses; i++ ) {
    vector<Ali_Frag*> tmp = get_active_frags(i);
    active_frags.insert( active_frags.end(), tmp.begin(), tmp.end() );
  }

  for( int i=0; i<(int)active_frags.size(); i++ ) {

    if( !fs.get_frag( active_frags[i]->get_id() )->is_active() ) {
      res.push_back( active_frags[i] );
    }

  }

  return res;
}

// GENERAL

void Frag_Set::clear_all()
{
  Frag_Columns.clear();
}

void Frag_Set::add_column( SSE_Frag_Set set )
{
  //  cerr << "add_column: top" << endl;
  Frag_Columns.push_back( set );
  //  cerr << "add_column: end" << endl;
}


void Frag_Set::set_num_sses()
{
  // determines the number of sses in the Frag_Set, only call once all columns are added
  // assumes the first and last column are the N- and C-terminal caps

  num_sses = Frag_Columns.size() - 2;
}


void Frag_Set::activate_terminal_caps()
{

  // activate N-terminal cap
  Frag_Columns.front().Frags.front().make_active();

  // activate C-terminal cap
  Frag_Columns.back().Frags.front().make_active();

  // NOTE:  N- and C-terminal columns are assumed to only hold one fragment (the cap)

  set_num_sses(); // NOTE: is there a more logical place to do this?
}

void Frag_Set::initialize_all_zscores()
{
  for( int i=1; i<(int)Frag_Columns.size()-1; i++ ) {
    Frag_Columns[i].set_frag_zscores();
  }
}

void Frag_Set::seed_all_columns( int num_frags )
{
  cerr << "seed_all_columns: num_sses: " << num_sses << endl;

  for( int i=1; i<=num_sses; i++ ) {
    Frag_Columns[i].activate_top_available_frag();
  }

}


void Frag_Set::count_frag_children()
{
  for( int i=num_sses; i>=0; i-- ) {
    for( int j=0; j<num_frags_in_sse(i); j++ ) {

      Ali_Frag* curr_frag = get_frag( i, j );

      unsigned long long int sum_children = 0;

      for( int k=0; k<(int)curr_frag->num_next(); k++ ) {
	Ali_Frag* next_frag = get_frag( curr_frag->get_next(k).next_frag );
	sum_children += ( 1 + next_frag->get_num_children() );
      }
      curr_frag->set_num_children( sum_children );
    }
  }
}


float Frag_Set::activate_next_best_available_frag()
{
  // activate the highest scoring availabe frag in this Frag_Set, return its z-score

  float max_zscore(-9999.f); // arbitrary low_number
  int max_sse_id(-1); // flag to catch error

  for( int i=1; i<=num_sses; i++ ) {

    float cur_zscore = Frag_Columns[i].get_highest_available_frag_zscore();

    if( max_zscore < cur_zscore ) {
      max_zscore = cur_zscore;
      max_sse_id = Frag_Columns[i].sse_id;
    }
  }

  if( max_sse_id == -1 ) {
    cerr << "Could not find a highest-scoring available frag. Exiting." << endl;
    exit(-1);
  }

  // activate the highest-scoring frag in Frag_Columns[ max_sse_id ]
  Frag_Columns[ max_sse_id ].activate_top_available_frag();

  return max_zscore;
}


bool Frag_Set::an_available_frag_exists()
{
  // determine whether an available frag remains in this Frag_Set

  for( int i=1; i<=num_sses; i++ ) {
    if( Frag_Columns[i].an_available_frag_exists() ) { return true; }
  }

  return false;

}


vector<Ali_Frag*> Frag_Set::get_active_frags( int i )
{
  return( Frag_Columns[i].get_active_frags() );
}


vector<Ali_Frag*> Frag_Set::get_available_frags( int i )
{
  return( Frag_Columns[i].get_available_frags() );
}


vector<Ali_Frag*> Frag_Set::export_all_frags()
{
  // exports all frags in this Frag_Set (except the terminal caps)

  vector<Ali_Frag*> res;

  for( int i=1; i<=num_sses; i++ ) {
    vector<Ali_Frag*> tmp = Frag_Columns[i].get_active_frags();
    res.insert( res.end(), tmp.begin(), tmp.end() );
    tmp = Frag_Columns[i].get_available_frags();
    res.insert( res.end(), tmp.begin(), tmp.end() );
  }

  return res;
}


bool Frag_Set::frags_in_order( int t_prev_end, int q_prev_end, int t_next_beg, int q_next_beg )
{
  return( ( q_next_beg > ( q_prev_end + 1 ) ) && ( t_next_beg > ( t_prev_end + 1 ) ) );
}


bool Frag_Set::frags_in_order( Ali_Frag* af1, Ali_Frag* af2 )
{
  return( ( af1->core_t1()+1 < af2->core_t0() ) && ( af1->core_q1()+1 < af2->core_q0() ) );
}


// REGION-FOCUSED
vector<Ali_Frag*> Frag_Set::get_active_frags_in_region( int sse_l, int sse_r, int q_top, int q_bot )
{
  vector<Ali_Frag*> res;

  for( int i=sse_l+1; i<sse_r; i++ ) {
    vector<Ali_Frag*> tmp = get_active_frags(i);

    for( int j=0; j<(int)tmp.size(); j++ ) {
      if( ( tmp[j]->is_active() ) && 
	  ( tmp[j]->sse_id > sse_l ) &&
	  ( tmp[j]->sse_id < sse_r ) &&
	  ( tmp[j]->core_q0() > q_top ) &&
	  ( tmp[j]->core_q1() < q_bot ) ) {
	res.push_back( tmp[j] );
      }
      
    }
  }

  return res;
}


vector<Ali_Frag*> Frag_Set::get_available_frags_in_region( int sse_l, int sse_r, int q_top, int q_bot )
{
  vector<Ali_Frag*> res;

  for( int i=sse_l+1; i<sse_r; i++ ) {
    vector<Ali_Frag*> tmp = get_available_frags(i);
    
    for( int j=0; j<(int)tmp.size(); j++ ) {
      if( ( tmp[j]->is_available() ) && 
	  ( tmp[j]->sse_id > sse_l ) &&
	  ( tmp[j]->sse_id < sse_r ) &&
	  ( tmp[j]->core_q0() > q_top ) &&
	  ( tmp[j]->core_q1() < q_bot ) ) {
	res.push_back( tmp[j] );
      }

    }
  }
    
  return res;
}


// ACCESS
int Frag_Set::num_frags_in_sse( int sse )
{
  return( Frag_Columns[sse].get_num_active_frags() );
}


Ali_Frag* Frag_Set::get_frag( Frag_ID f )
{
  return( Frag_Columns[f.sse_idx].get_frag(f.frag_idx) );
}

Ali_Frag* Frag_Set::get_frag( int sse_id, int frag_id )
{
  return( Frag_Columns[sse_id].get_frag(frag_id) );
}


// DEBUGGING
void Frag_Set::print_all_frags()
{
  for( int i=0; i<(int)Frag_Columns.size(); i++ ) {
    Frag_Columns[i].print_frag_scores( cerr );
    cerr << endl;
  }
}

void Frag_Set::print_all_sses()
{
  for( int i=0; i<(int)Frag_Columns.size(); i++ ) {
    Frag_Columns[i].print_sse_info( cerr );
    cerr << endl;
  }
}


void Frag_Set::print_all_frag_connections()
{
  for( int i=1; i<=num_sses; i++ ) {

    cerr << "SSE_ID: " << i << endl << endl;

    vector<Ali_Frag*> active_frags = get_active_frags(i);

    for( int j=0; j<(int)active_frags.size(); j++ ) {

      int tmp_sse_id = get_frag( i, j )->sse_id;
      int tmp_frag_id = get_frag( i, j )->frag_id;
      
      cerr << "Frag: (" << tmp_sse_id << "," << tmp_frag_id << "):" << endl;
      get_frag( i, j )->print(query_seq, templ_seq);      
      
      cerr << "Connections: ";
      
      for( int k=0; k<(int)active_frags[j]->num_next(); k++ ) {
	Frag_ID f = active_frags[j]->get_next(k).next_frag;
	cerr << "(" << f.sse_idx << "," << f.frag_idx << ") ";
      }
      
      cerr << endl << endl;
    }
    
    cerr << endl;
    cerr << "--------------------------------" << endl;
  }

}


void Frag_Set::report_statistics()
{
  int num_frags, ali_len, max_gap;
  float avg_gap;
  int sum_frags(0), sum_ali_len(0), max_gap_overall(-1);

  for( int i=1; i<=num_sses; i++ ) {
    get_col(i)->get_statistics( num_frags, ali_len, avg_gap, max_gap );

    cerr << "Stats: sse: " << i << " - #frags: " << num_frags
	 << ", avg gap: " << avg_gap
	 << ", max gap: " << max_gap << endl;

    sum_frags += num_frags;
    sum_ali_len += ali_len;
    if( max_gap > max_gap_overall ) { max_gap_overall = max_gap; }
  }

  cerr << "Stats (overall): #frags: " << sum_frags
       << ", avg gap: " << (float)( sum_ali_len - sum_frags ) / (float)( sum_frags + num_sses )
       << ", max gap: " << max_gap_overall << endl;

}
