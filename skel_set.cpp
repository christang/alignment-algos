#include "skel_set.h"

using namespace std;

/*********************************************
* implementation of Skel_Set class *
*********************************************/

// constructor 1
Skel_Set::Skel_Set( int min_ali, float min_CO,
		    int max_kept, float max_size,
		    Frag_Set* fs, Ali_Str_Info* asi, Alignment_Strand_Evaluator* ase, Ali_Dist* adm )
{

  cerr << "Skel_Set constructor: top" << endl;

  // load paramters
  Frags = fs;
  Str_Data = asi;
  Strand_Eval = ase;
  Measurer = adm;

  min_aligned_residues = min_ali;
  min_SSE_CO_fraction = min_CO;

  max_alis = max_kept;
  max_cluster_size = max_size;

  max_bad_alis = 100;

  templ_seq = Str_Data->get_templ_seq();
  query_seq = Str_Data->get_query_seq();

  num_culled_by_coverage = 0;
  num_culled_by_contact_order = 0;
  num_culled_by_strand_rules = 0;
  num_culled_by_score = 0;

  top_constrained_skel = NULL;

  tracking_mode = ( adm != NULL );

  if( tracking_mode ) {
    ofs_cov.open( "track_low_coverage.txt" );
    ofs_sse_co.open( "track_low_CO.txt" );
    ofs_strand_rules.open( "track_bad_strands.txt" );
    ofs_score.open( "track_low_score.txt" );
  }

  skels_created = 0; // keep track of skels in memory to check for memory leaks
  skels_deleted = 0;

  for( int i=0; i<(int)get_frag(0,0)->num_next(); i++ ) { // get starting skeletons
    Start_Skels.push_back( new Skel_Ali( get_frag(0,0)->get_next(i),
					 Str_Data, Frags,
					 0 ) );
    skels_created++;
  }

  // set minimum SSE contact order threshold for alignments
  float template_SSE_CO = find_template_SSE_CO();

  cerr << "Template SSE_CO: " << template_SSE_CO << endl;
  min_SSE_CO = min_SSE_CO_fraction * template_SSE_CO;
  cerr << "Minimum SSE_CO: " << min_SSE_CO << endl; 

  cerr << "Skel_Set constructor: end" << endl;
}

Skel_Set::~Skel_Set()
{
  // delete all skel alis in memory

  while( !Top_Skels.empty() ) {
    delete Top_Skels.back();
    Top_Skels.pop_back();
    skels_deleted++;
  }

  while( !Low_Coverage.empty() ) {
    delete Low_Coverage.back();
    Low_Coverage.pop_back();
    skels_deleted++;
  }

  while( !Low_SSE_CO.empty() ) {
    delete Low_SSE_CO.back();
    Low_SSE_CO.pop_back();
    skels_deleted++;
  }

  while( !Bad_Strands.empty() ) {
    delete Bad_Strands.back();
    Bad_Strands.pop_back();
    skels_deleted++;
  }

  while( !Low_Score.empty() ) {
    delete Low_Score.back();
    Low_Score.pop_back();
    skels_deleted++;
  }

  cerr << "Skel_set destructor: skels_created: " << skels_created << ", skels_deleted: " << skels_deleted << endl;
}

// GENERAL


void Skel_Set::find_top_skeletons()
{
  for( int i=0; i<(int)Start_Skels.size(); i++ ) {
    //    cerr << "fts: i: " << i << endl;
    grow_skel( Start_Skels[i] );
  }

  cerr << "Num culled by coverage: " << num_culled_by_coverage << endl;
  cerr << "Num culled by contact order: " << num_culled_by_contact_order << endl;
  cerr << "Num culled by strand rules: " << num_culled_by_strand_rules << endl;
  cerr << "Num culled by score: " << num_culled_by_score << endl;

  num_culled_by_coverage = 0;
  num_culled_by_contact_order = 0;
  num_culled_by_strand_rules = 0;
  num_culled_by_score = 0;

}


void Skel_Set::find_top_constrained_skel( Skel_Ali* orig )
{
  top_constrained_skel = NULL;

  Frag_ID orig_first = orig->get_connection(0).next_frag;

  for( int i=0; i<get_frag(0,0)->num_next(); i++ ) {

    Frag_Connection tmp_fc = get_frag(0,0)->get_next(i);

    if( ( tmp_fc.next_frag.sse_idx < orig_first.sse_idx ) ||
	( ( tmp_fc.next_frag.sse_idx == orig_first.sse_idx ) &&
	  ( tmp_fc.next_frag.frag_idx == orig_first.frag_idx ) ) ) {

      Skel_Ali* sa = new Skel_Ali( tmp_fc, orig->get_str_info(), orig->get_frag_set(), 0 );
      skels_created++;
      grow_constrained_skel( sa, orig, 1 );

    }

  }

  if( top_constrained_skel == NULL ) {
    //    cerr << "After find_top_constrained_skel, pointer is still NULL." << endl;
    //    cerr << "Returning original skel" << endl;
    //    pause();
    //    return orig;

    cerr << "Something's strange: grow_constrained_skel did not find the original skel. Exiting." << endl;
    exit(-1);
  }

}


void Skel_Set::grow_skel( Skel_Ali* sa )
{
  // Adds an extra C-terminal fragment to 'sa' from the list of the current last fragment's
  // list of connections. If sa's last fragment is the C-terminus, stop recursing and evaluate
  // the alignment.

  // LOW COVERAGE ELIMINATION CASE
  if( sa->get_num_aligned() + Str_Data->get_TSR_to_C( sa->get_last_templ_res_idx() )
      < min_aligned_residues ) {
    // if sa's coverage is so low that it cannot surpass the minimum even if it aligns to every 
    // upcoming template SSE residue, then stop this branch now
    pre_empt_low_coverage( sa );
    return;
  }

  // TERMINATION CASE (last frag is C-terminus)
  if( sa->last_frag_is_C_terminal() ) {
    handle_completed_skel( sa );
    return; // Done.  This skel does not need to be looked at or recursed further.
  }
  else{ // RECURSIVE CASE (there are more frags to add)

    Frag_Connection curr = sa->get_last_connection(); // the last frag thus far in the ali

    for( int i=0; i<(int)get_frag( curr.next_frag )->num_next(); i++ ) {

      Frag_Connection tmp_fc = get_frag( curr.next_frag )->get_next(i);
      Skel_Ali* new_skel_ali = new Skel_Ali( sa ); // start a new skeleton
      new_skel_ali->add_connection( tmp_fc ); // elongate the skeleton
      skels_created++;
      
      grow_skel( new_skel_ali ); // recurse
    }
    
  }

  delete sa; // The skel has spawned all its children, including one that ends right here.
             // (i.e., sa does not end in a C-terminal frag.)  sa has thus run its course,
             // so delete it.
  skels_deleted++;
}


void Skel_Set::grow_constrained_skel( Skel_Ali* sa, Skel_Ali* orig, int post_idx )
{
  // similar to 'grow_skel' above, but 'sa' must grow only between the fragments already present
  // in 'orig', which is a complete, capped skel ali.  'post' is the next fragment that 'sa' must
  // grow toward without skipping

  // NOTE: don't need low coverage elimination here because 'orig' already passed coverage test

  // TERMINATION CASE (last frag is C-terminus)
  if( sa->last_frag_is_C_terminal() ) {
    handle_completed_constrained_skel( sa );
    return; // Done.  This skel does not need to be looked at or recursed further.
  }
  else{ // RECURSIVE CASE (there are more frags to add)

    Ali_Frag* post = get_frag( orig->get_connection( post_idx ).next_frag );
    
    Frag_Connection curr = sa->get_last_connection(); // the last frag thus far in the ali

    for( int i=0; i<(int)get_frag( curr.next_frag )->num_next(); i++ ) {

      Frag_Connection tmp_fc = get_frag( curr.next_frag )->get_next(i);

      // NOTE: could merge the first two conditions with a '>' function for Ali_Frag

      if( tmp_fc.next_frag.sse_idx > post->sse_id ) {
	// since next_frags are in order according to sse_id, once post->sse_id is passed ...
	// all the following next_frags are past too.  so just break.
	break;
      }

      if( ( tmp_fc.next_frag.sse_idx == post->sse_id ) &&
	  ( tmp_fc.next_frag.frag_idx > post->frag_id ) ) {
	// since next_frags are also in order according to frag_id (within an sse_id), once ...
	// post->frag_id is passed, all the following next_frags are bad too.  so break.
	break;
      }

      if( ( tmp_fc.next_frag.sse_idx == post->sse_id ) &&
	  ( tmp_fc.next_frag.frag_idx < post->frag_id ) ) {
	// next_frag is in the same SSE as post, but is a different frag
	continue;
      }

      // NOTE: could simplify this with an '==' function for Ali_Frag

      if( !( ( tmp_fc.next_frag.sse_idx == post->sse_id ) &&
	     ( tmp_fc.next_frag.frag_idx == post->frag_id ) ) && 
	  ( !Frags->frags_in_order( get_frag( tmp_fc.next_frag ), post ) ) ) {
	// even if next_frag has a lower sse_id than post, if it does not precede it ...
	// in the query sequence, post will never be reached later.  so skip it. ...
	// but first must check to make sure it is not the same frag in the same SSE.
	continue;
      }

      Skel_Ali* new_skel_ali = new Skel_Ali( sa ); // start a new skeleton
      new_skel_ali->add_connection( tmp_fc ); // elongate the skeleton
      skels_created++;
      
      int next_post_idx = find_next_post( new_skel_ali, orig, post_idx );
      
      grow_constrained_skel( new_skel_ali, orig, next_post_idx ); // recurse
    }
  }
  
  delete sa; // The skel has spawned all its children, including one that ends right here.
             // (i.e., sa does not end in a C-terminal frag.)  sa has thus run its course,
             // so delete it.
  skels_deleted++;
  
}


void Skel_Set::pre_empt_low_coverage( Skel_Ali* sa )
{

  if( tracking_mode && 
      ( (float)sa->get_num_aligned() > ( 0.75f * (float)min_aligned_residues ) ) ) {
    // only bother tracking 'sa' if 'tracking_mode' is on and coverage is above 75% of cutoff

    // pre-emptively cap sa with a C-term-frag, if not already capped
    if( !sa->last_frag_is_C_terminal() ) {
      Frag_Connection last_fc = sa->get_last_connection(); // the last frag thus far in the ali
      Frag_Connection cap_fc = get_frag( last_fc.next_frag )->get_last_next();
      sa->add_connection( cap_fc ); // elongate the skeleton
    }
    
    // sa is now a full alignment with a C-term cap, ready to be measured against the native
    handle_culled_skel_ali( sa, 1 ); // 1=low coverage
  }
  else { // don't bother tracking, just delete sa
    delete sa;
    skels_deleted++;
  }

}


void Skel_Set::handle_completed_skel( Skel_Ali* sa )
{
  // once a skel has been completed (i.e., a C-terminal cap is at the end),
  // this function checks all the filters and moves to 'sort_top_skels' or
  // 'handle_culled_skel_ali' depending on whether it passes and whether
  // tracking_mode is on or off

  sa->calc_skel_SSE_CO();

  // check critera for accepting alignment
  int reason(-1);
  
  if( passes_all_filters( sa, reason ) ) {

    find_top_constrained_skel( sa ); // fill in loops, maybe

    if( *top_constrained_skel != *sa ) { // avoid duplicates, top_constrained_skel will be found later 

      delete sa;
      delete top_constrained_skel;
      skels_deleted += 2;
    }
    else {

      delete top_constrained_skel;
      skels_deleted++;

      sa->set_param( sa->get_score() );
      sort_top_skels( sa );
    }
  }
  else { // did not pass all rules
    
    if( tracking_mode ) {
      handle_culled_skel_ali( sa, reason );
    }
    else { // no tracking, so just delete skel
      delete sa;
      skels_deleted++;
    }
  }

}


void Skel_Set::handle_completed_constrained_skel( Skel_Ali* sa )
{
  // like 'handle_completed_skel' above, but doesn't bother with tracking_mode
  // or sort_top_skels.  This function is only interested in keepint the top
  // constrainted skel ali.

  sa->calc_skel_SSE_CO();

  // check critera for accepting alignment
  int reason(-1);
  
  if( passes_all_filters( sa, reason ) ) {
    sa->set_param( sa->get_score() );

    if( top_constrained_skel != NULL ) {

      if( sa->get_score() > top_constrained_skel->get_score() ) {

	// copy 'sa' onto 'top_constrained_skel'
	delete top_constrained_skel;
	skels_deleted++;

	top_constrained_skel = sa;
      }
      else { // 'sa' is lower-scoring, so delete it
	delete sa;
	skels_deleted++;
      }

    }
    else { // top_constrained_skel is empty, so fill it with sa
      top_constrained_skel = sa;
    }

  }
  else { // did not pass all rules
    
    delete sa;
    skels_deleted++;
  }
  

}
 

int Skel_Set::find_next_post( Skel_Ali* curr, Skel_Ali* orig, int old_post_idx )
{
  // return an int representing the next connection in 'orig' whose next_frag is the ...
  // first one  not in 'curr'.
  // assume that 'curr' has not skipped any fragments in 'orig' up to now

  Ali_Frag* curr_last_frag = get_frag( curr->get_last_connection().next_frag );
  Ali_Frag* old_post_frag = get_frag( orig->get_connection( old_post_idx ).next_frag );

  if( curr_last_frag->sse_id < old_post_frag->sse_id ) {
    return old_post_idx;
  }

  if( curr_last_frag->sse_id == old_post_frag->sse_id ) { // last frag and post in same SSE

    if( curr_last_frag->frag_id == old_post_frag->frag_id ) { // last frag is post
      return( old_post_idx + 1 );
    }
    else{ // error! last frag and post cannot be in same SSE and not the same frag
      cerr << "Error: Frag in skel ali is in same SSE, but different frag than post. Exiting." << endl;
      exit(-1);
    }

  }
  else { // error! last frag's sse_id must be greater than post's
    cerr << "Error: Frag in skel ali has passed that in post. Exiting." << endl;
    exit(-1);
  }

}


bool Skel_Set::passes_all_filters( Skel_Ali* sa, int& reason )
{
  // return true if 'sa' passes all the criteria.  If false, set 'reason' to reflect why.

    if( sa->get_num_aligned() < min_aligned_residues ) {
      reason = 1;
      return false;
    }

    if( sa->get_contact_order() < min_SSE_CO ) {
      reason = 2;
      return false;
    }

    if( Strand_Eval->ali_passes_rules( sa->get_sse_id_list() ) ) {
      reason = 3;
      return false;
    }

    return true;
}


void Skel_Set::sort_top_skels( Skel_Ali* sa )
{

  list<Skel_Ali*>::reverse_iterator rit = Top_Skels.rbegin();

  while( ( rit !=Top_Skels.rend() ) && ( (*rit)->get_param() < sa->get_param() ) ) {
    rit++;
  }
  Top_Skels.insert( rit.base(), sa );

  if( (int)Top_Skels.size() > max_alis ) {

    Skel_Ali* last_skel = Top_Skels.back();

    if( tracking_mode ) {
      // last_skel is being culled for the first time by being pushed out of Top_Skels for its low score
      handle_culled_skel_ali( last_skel, 4 ); // sort it into list of low-scoring skels
    }
    else { // we're not tracking so just delete it
      delete last_skel; // free up the memory from this alignment
      skels_deleted++;
    }

    Top_Skels.pop_back();  // erase from the list of top alignments, whether tracking or not
  }

}

void Skel_Set::sort_culled_skels( Skel_Ali* sa, list<Skel_Ali*>& skel_list, int max )
{

  list<Skel_Ali*>::reverse_iterator rit = skel_list.rbegin();

  while( ( rit != skel_list.rend() ) && ( (*rit)->get_param() > sa->get_param() ) ) {
    rit++;
  }
  skel_list.insert( rit.base(), sa );

  if( (int)skel_list.size() > max ) {

    // can assume tracking_mode is true (wouldn't call this function otherwise)
    delete skel_list.back(); // free up the memory from this alignment
    skel_list.pop_back();    // erase from the list of top alignments
    skels_deleted++;

  }

}


void Skel_Set::handle_culled_skel_ali( Skel_Ali* sa, int reason )
{
  Measurer->load_test( sa->export_vrp() );

  float dist = Measurer->get_dist_between_main_and_test();

  sa->set_shift( dist );
  sa->set_param( sa->get_shift() );

  switch( reason ) {
  case 1:
    sort_culled_skels( sa, Low_Coverage, max_bad_alis );
    num_culled_by_coverage++;
    break;
  case 2:
    sort_culled_skels( sa, Low_SSE_CO, max_bad_alis );
    num_culled_by_contact_order++;
    break;
  case 3:
    sort_culled_skels( sa, Bad_Strands, max_bad_alis );
    num_culled_by_strand_rules++;
    break;
  case 4:
    sort_culled_skels( sa, Low_Score, max_bad_alis );
    num_culled_by_score++;
    break;
  default:
    cerr << "No valid reason given for this ali's exclusion." << endl;
  }

}


float Skel_Set::find_template_SSE_CO()
{
  vector<bool> contacting_residues;
  contacting_residues.resize( templ_seq.size(), false );

  int num_residues_in_contact(0);

  for( int i=0; i<Str_Data->get_num_templ_SSEs(); i++ ) {
    for( int j=Str_Data->get_SSE_Data(i).beg_id; j<Str_Data->get_SSE_Data(i).end_id; j++ ) {

      for( int m=0; m<Str_Data->get_num_templ_SSEs(); m++ ) {

	if ( m == i ) { continue; } // ignore intra-SSE contacts

	for( int n=Str_Data->get_SSE_Data(m).beg_id; n<Str_Data->get_SSE_Data(m).end_id; n++ ) {

	  if( j == n ) { continue; }

	  if( Str_Data->get_contact( j, n ) ) {

	    if( !contacting_residues[j] ) {
	      contacting_residues[j] = true;
	      num_residues_in_contact++;
	    }

	    if( !contacting_residues[n] ) {
	      contacting_residues[n] = true;
	      num_residues_in_contact++;
	    }

	  }
	}
      }
    }
  }

  int num_SSE_residues(0);

  for( int i=0; i<Str_Data->get_num_templ_SSEs(); i++ ) {
    num_SSE_residues += ( Str_Data->get_SSE_Data(i).end_id - Str_Data->get_SSE_Data(i).beg_id + 1 );
  }

  return( (float)num_residues_in_contact / (float)num_SSE_residues );
}


void Skel_Set::send_culled_alis_to_files()
{

  cerr << "Low_Coverage" << endl;
  int idx(1);
  for( list<Skel_Ali*>::iterator it=Low_Coverage.begin(); it!=Low_Coverage.end(); it++ ) {
    (*it)->print( query_seq, templ_seq, min_aligned_residues, ofs_cov );
    cerr << "shift: " << (*it)->get_shift() 
	 << ", coverage: " << (*it)->get_num_aligned() << " of " << min_aligned_residues << endl;
    idx++;
  }
  cerr << endl << endl;

  cerr << "Low_SSE_CO" << endl;
  idx = 0;
  for( list<Skel_Ali*>::iterator it=Low_SSE_CO.begin(); it!=Low_SSE_CO.end(); it++ ) {
    (*it)->print( query_seq, templ_seq, min_aligned_residues, ofs_sse_co );
    cerr << "shift: " << (*it)->get_shift() 
	 << ", SSE_CO: " << (*it)->get_contact_order() << " of " << min_SSE_CO << endl;
    idx++;
  }
  cerr << endl << endl;

  cerr << "Bad_Strands" << endl;
  idx = 0;
  for( list<Skel_Ali*>::iterator it=Bad_Strands.begin(); it!=Bad_Strands.end(); it++ ) {
    (*it)->print( query_seq, templ_seq, min_aligned_residues, ofs_strand_rules );
    cerr << "shift: " << (*it)->get_shift() << endl;
    idx++;
  }
  cerr << endl << endl;

  cerr << "Low_Score" << endl;
  idx = 0;
  for( list<Skel_Ali*>::iterator it=Low_Score.begin(); it!=Low_Score.end(); it++ ) {
    (*it)->print( query_seq, templ_seq, min_aligned_residues, ofs_score );
    cerr << "shift: " << (*it)->get_shift()
	 << ", score: " << (*it)->get_score() << endl;
    idx++;
  }
  cerr << endl << endl;

}


void Skel_Set::cluster_alignments()
{
/*
  // Purpose: perform UPGMA clustering on all alignments in Top_Skels

  // transfer the Top_Skels list to the Skel_Alis vector
  vector<Skel_Ali*> Skel_Alis;

  while( !Top_Skels.empty() ) {
    Skel_Alis.push_back( Top_Skels.front() );
    Top_Skels.pop_front();
  }

  get_exact_inter_ali_areas( Skel_Alis );

  clusterer = new UPGMA_Clusterer( inter_ali_area, Skel_Alis.size() );

  cerr << "setup new clusterer" << endl;

  clusterer->cluster();

  cerr << "performed clustering" << endl;

  cerr << "cluster_alignments: max_cluster_size: " << max_cluster_size << endl;

  clusterer->find_clusters_under_threshold( max_cluster_size );
  //  clusterer->print_clusters_under_threshold();

  cerr << "cluster_alignments: # clusters found: " << clusterer->get_num_clusters() << endl;

  // arbitrarily selects the first member of each cluster to represent the cluster
  // NOTE: should select the centroid or highest-scoring member of each cluster
  for( int i=0; i<clusterer->get_num_clusters(); i++ ) {
    Top_Skels.push_back( Skel_Alis[ clusterer->get_member_index( i, 0 ) ] );
  }

  // sort Top_Skels
  list<Skel_Ali*> Sorted_Skels;

  while( !Top_Skels.empty() ) {

    Skel_Ali* tmp = Top_Skels.front();
    list<Skel_Ali*>::iterator jt = Sorted_Skels.begin();

    while( ( jt != Sorted_Skels.end() ) && ( tmp->get_score() < (*jt)->get_score() ) ) {
      jt++;
    }
    // jt now points to the position tmp should be inserted before

    Sorted_Skels.insert( jt, tmp ); // insert the front of Top_Skels into the
                                    // correct spot in Sorted_Skels

    Top_Skels.pop_front(); // remove the inserted element and shorten Top_Skels
  }

  // Sorted_Skels should now be a sorted (descending) version of Top_Skels, so swap them
  Top_Skels.swap( Sorted_Skels );
*/
}


void Skel_Set::get_exact_inter_ali_areas( const vector<Skel_Ali*>& Alis )
{

  // determine the distances between skel alis and fill in inter_ali_areas
  inter_ali_area = new float* [ (int)Alis.size() ];

  for( int i=0; i<(int)Alis.size(); i++ ) {
    inter_ali_area[i] = new float [ i+1 ];
  }

  //////////
  Ali_Dist tmp;
  int a1(21), a2(65);
  cerr << "Loading main: ali " << a1 << endl;
  tmp.load_main( Alis[a1]->export_vrp() );
  cerr << "Loading test: ali " << a2 << endl;
  tmp.load_test( Alis[a2]->export_vrp() );

  cerr << "Distance between " << a1 << " and " << a2 << ": "
       << tmp.get_area_between_main_and_test() << endl;

  pause();
  ////////

  //////////
  Ali_Dist tmp2;
  a1 = 65;
  a2 = 21;
  cerr << "Loading main: ali " << a1 << endl;
  tmp2.load_main( Alis[a1]->export_vrp() );
  cerr << "Loading test: ali " << a2 << endl;
  tmp2.load_test( Alis[a2]->export_vrp() );

  cerr << "Distance between " << a1 << " and " << a2 << ": "
       << tmp2.get_area_between_main_and_test() << endl;

  pause();
  ////////


  float last_percentage(-0.01f), curr_percentage;

  for( int i=0; i<(int)Alis.size(); i++ ) {

    curr_percentage = pow( ( (float)i / (float)Alis.size() ), 2 );

    if( ( curr_percentage - last_percentage ) > 0.01f ) {
      last_percentage = curr_percentage;
      cerr << "finding distances: " << curr_percentage*100 << "% done." << endl;
    }

    Ali_Dist X;

    X.load_main( Alis[i]->export_vrp() );

    for( int j=0; j<i; j++ ) {

      X.load_test( Alis[j]->export_vrp() );

      inter_ali_area[i][j] = ( X.get_area_between_main_and_test() );

      if( inter_ali_area[i][j] < 0.f ) {
	cerr << "invalid area measurement between alis " 
	     << i << " and " << j << ": " << inter_ali_area[i][j] << endl;

	pause();
       }
    }

    inter_ali_area[i][i] = 0.f; // don't bother calculating this one
  }

  cerr << "get_exact: end: done calculating distances." << endl;
}


// ACCESS

Ali_Frag* Skel_Set::get_frag( Frag_ID f ) const
{
  //return source_fm->get_frag( f );
  return( Frags->get_frag( f ) );
}

Ali_Frag* Skel_Set::get_frag( int s, int f ) const
{
  //  return source_fm->get_frag( s, f );
  return( Frags->get_frag( s, f ) );
}


// DEBUGGING

void Skel_Set::print_Top_Skels()
{
  cerr << "Printing " << Top_Skels.size() << " skels of Top_Skels:" << endl;

  int index=0;

  list<Skel_Ali*>::iterator it = Top_Skels.begin();

  while( it != Top_Skels.end() ) {
    cerr << "Top_Skel[" << index << "]" << endl;
    (*it)->print( query_seq, templ_seq, min_aligned_residues );
    it++;
    index++;
  }

}


void Skel_Set::pause()
{
  cerr << "Paused. Enter a character and press Enter to continue." << endl;
  string pause;
  cin >> pause;
}
