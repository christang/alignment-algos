/**
 *  Package HMAP2.1
 *  File: ssss.h
 *  Desc: Search for sub-optimal alignments based on 
 *        sampling shifts in alignments to template SSEs
 *
 *  Created on 12/18/06
 *  Author: kuziemko @ honig lab
 *
 */

#ifndef _HMAP2_SSSS
#define _HMAP2_SSSS

#include <string>
#include <fstream>
#include <sstream>

#include "alignment.h"
#include "clusterset.h"
#include "enumerator.h"
#include "dpmatrix.h"
#include "kmedoidclusterer.h"
#include "optimal_subali.h"
#include "noalib.h"
#include "sflags.h"

#include "ali_str_info.h"
#include "frag_set.h"
#include "frag_matrix.h"
#include "skel_set.h"
#include "skel_ali.h"
#include "ali_frag.h"
#include "ali_dist.h"
#include "ssss_shared_defs.h"

using namespace std;


template <class S1, class S2, class Etype>
class SSSS : public Enumerator<S1,S2,Etype> {
  
  // SSSS = Sample Shifts in Secondary Structures
  // This class is a sub-optimal alignment enumerator.

  // It has two phases.  First, it scans through a DP matrix
  // looking for all possible alignments to each template secondary
  // structure element tht would lead to decent (above a threshold)
  // global alignments.  The results is a list of alignment fragments
  // for each template SSE.  The second phase uses a different scoring
  // function to recombine these fragments into full alignments and
  // returns the top N alignments.


  // The structs below are a bit confusing.  An 'SSE_Frag' is made
  // for each alignment fragment found for a particular SSE.

 public:
  
  typedef AlignedPairList<S1,S2> SingleAlignment;
  typedef AlignedPair<S1,S2>     SinglePair;

  SSSS (const NOaliParams&, const Evaluator<S1,S2,Etype>&, DPMatrix<S1,S2,Etype>* dpm_tmp,
	int, int, float, float, int, int, float, int, string );

  ~SSSS();

  inline int estimateSize () const { return params->number_suboptimal; } ;
  vector<Ali_Frag*> GetAlignedFrags();
  vector<Ali_Frag*> get_new_frags();
  void setNFill(int m) { nfill=m; }
  void enumerate (DPMatrix<S1,S2,Etype>& dpm_fwd,
		  AlignmentSet<S1,S2,Etype>& as);

  void build_alignments();

  //  void choose_fragments_for_ali();

 private:
  
  // reads flags file to get locations of template sses
  vector<SSE_Data> find_template_sses( int );

  void fill_frag_matrix();
  void fill_frag_matrix2();
  void setup_data_structures();

  int get_user_selection();
  void pause();

  // takes in a skeleton_ali and returns a full alignment for final printing
  SingleAlignment convert_skel_to_APL( Skel_Ali* );
  void output_pir_ali( Skel_Ali, int, DPMatrix<S1,S2,Etype>&,AlignmentSet<S1,S2,Etype>& as , ostream& os=cout);

  const NOaliParams* params;
  const Evaluator<S1,S2,Etype> evaluator;
  DPMatrix<S1,S2,Etype>* dpm;

  //  Alignment_Data* Info;
  Ali_Str_Info* Str_Data;
  Frag_Set *All_Frags, Old_Frags;
  Alignment_Strand_Evaluator* Strand_Eval;
  Frag_Matrix* Main_Frag_Selector;
  Skel_Set* Alignment_Builder;
  Ali_Dist* Dist_Measurer;

  unsigned int max_subopt; // # of alignments that user wants returned from 'sample_shifts'
  unsigned long long int max_alis_to_search;

  float min_coverage;
  float min_SSE_CO;

  int max_in_betw_shift;

  unsigned int ali_mode; // 0 = align fragments to entire template SSE
                         // 1 = align to part of template SSE
  
 
  unsigned int nfill;  // Number of fragments to add each time enumerate is called.

  float max_avg_shift;  // when clustering, the maximum number of residues one alignment
                        // of a cluster can be from another member

  int tracking;
  string native_ali; // filename of native alignment in fasta format

  bool tracking_mode;

  const S1* query;
  const S2* templ;

  string query_seq;
  string templ_seq;
  int query_len, templ_len; // sequence lengths

  //  int mode;

  //  float shorter_seq_len;

  float min_ali_residues;
  float max_contact_dist;

  vector<SSE_Data> SSEs;
  
  float** sims;
  float** cb_dists;
  bool** templ_contacts;

  bool* query_predicted_loops;

  list<Skel_Ali> Returned_Skel_Alis;

  int* TSR_to_Nterm;
  int* TSR_to_Cterm;

  bool** Strand_Pairings;

  float templ_SSE_CO;
  
  map<string, SingleAlignment> loops;  // to keep track of loop alignments using string containing endpoints as a key

  // debugging
  int ali_counter;

};

//constructor
template <class S1, class S2, class Etype>
SSSS<S1,S2,Etype>::
SSSS (const NOaliParams& p, const Evaluator<S1,S2,Etype>& eval, DPMatrix<S1,S2,Etype>* dpm_tmp,
      int num_alis_kept, int max_alis,
      float min_cov, float min_CO, 
      int max_frag_shift, int ali_how, 
      float max_cluster_shift,
      int track, string nat_ali )
  : params(&p), evaluator(eval), dpm(dpm_tmp),
     max_subopt(num_alis_kept), max_alis_to_search(max_alis),
     min_coverage( min_cov ), min_SSE_CO( min_CO ),
     max_in_betw_shift( max_frag_shift ),
     ali_mode( ali_how ),
     max_avg_shift( max_cluster_shift ),
     tracking( track ), native_ali( nat_ali )
{

  // determine lengths of query and template sequences
  query_len = dpm->getQuerySize() - 1;
  templ_len = dpm->getTemplateSize() - 1;

  // establish min_ali_residues (for low coverage filter)
  min_ali_residues = min_coverage * ( query_len - 1 );

  // setup query and template sequence strings
  query_seq = *dpm->getQuerySequence()->getString();
  templ_seq = *dpm->getTemplateSequence()->getString();

  // setup query and template Sequence variables
  query = dpm->getQuerySequence();
  templ = dpm->getTemplateSequence();

  max_contact_dist = 6.f; // default/hardcoded value defining the maximum distance between 2 contacting c-betas

  tracking_mode = ( tracking == 1 ); // turn tracking_mode on or off

  nfill=0;

  // initialize large class variables to NULL
  Str_Data = NULL;
  All_Frags = NULL;
  Strand_Eval = NULL;
  Main_Frag_Selector = NULL;
  Alignment_Builder = NULL;
  Dist_Measurer = NULL;

  Str_Data = new Ali_Str_Info();

  // fill in all arrays, matrices of needed data
  setup_data_structures();

  // setup Dist_Measurer for tracking culled alis and clustering, only if tracking is on (and native ali is given)
  if( native_ali != "" ) {
    Dist_Measurer = new Ali_Dist();
    Dist_Measurer->load_main( native_ali );
  }

  // setup the Frag_Set
  All_Frags = new Frag_Set();

  cerr << "Made empty Frag_Set" << endl;

  // setup the strand rule evaluator
  Strand_Eval = new Alignment_Strand_Evaluator();

  cerr << "Made empty Strand rule evaluator" << endl;

  // setup the main Frag_Matrix
  Main_Frag_Selector = new Frag_Matrix(	(int)min_ali_residues,
					All_Frags, Str_Data,
					max_in_betw_shift, ali_mode, Dist_Measurer );

  cerr << "Made new Frag_Matrix" << endl;

  Main_Frag_Selector->create_all_fragments( All_Frags ); // find top scoring frag in each shift in each SSE, place each frag in Available

  cerr << "created all fragments" << endl;

  All_Frags->initialize_all_zscores(); // calculate a z-scores for each fragment

  cerr << "initialized all z-scores" << endl;

  All_Frags->seed_all_columns( 1 ); // move the highest scoring fragment in each SSE to Active

  cerr << "seeded each column" << endl;

  Main_Frag_Selector->find_fragment_connections( All_Frags ); // establish connections between fragments

  cerr << "found initial fragment connections" << endl;

  All_Frags->count_frag_children(); // determine how many alignments proceed from each fragment

  cerr << "done counting frag_children" << endl;

  // establish strand pairing evaluator
  Strand_Eval->load_SSE_contacts( (int)SSEs.size()+2, Strand_Pairings );
  Strand_Eval->load_All_Strands( SSEs );
  //  Main_Frag_Selector->determine_strand_pairing_rules();
  Strand_Eval->determine_rules();

  cerr << "done determining strand rules" << endl;

  ali_counter = 0;

  // print paramter values
  cerr << "Max alis to return: ...... " << max_subopt << endl;
  cerr << "Max alis to search: ...... " << max_alis_to_search << endl;
  cerr << "Minimum frac coverage: ... " << min_coverage << endl;
  cerr << "Minimum frac SSE_CO: ..... " << min_SSE_CO << endl;
  cerr << "Full vs. Partial ali: .... " << ali_mode << endl;
  cerr << "Max cluster size: ........ " << max_avg_shift << endl;
  cerr << "Tracking mode ............ " << tracking_mode << endl;
  cerr << "Native alignment ......... " << native_ali << endl;
}

//destructor
template <class S1, class S2, class Etype>
SSSS<S1,S2,Etype>::~SSSS()
{

  if( Main_Frag_Selector!=NULL ) delete Main_Frag_Selector;
  else return;

  Main_Frag_Selector=NULL;

  // delete similarity matrix, sims
  for( int i=0; i<(int)query_seq.size(); i++ ) {
    delete sims[i];
  }
  delete sims;

  // delete the inter-residue distance matrix
  for( int i=0; i<(int)templ_seq.size(); i++ ) {
    delete cb_dists[i];
  }
  delete cb_dists;

  // delete the contact matrix
  for( int i=0; i<(int)templ_seq.size(); i++ ) {
    delete templ_contacts[i];
  }
  delete templ_contacts;

  // delete the query_predicted_loops array
  delete query_predicted_loops;

  // delete the TSR arrays
  delete TSR_to_Nterm;
  delete TSR_to_Cterm;

  // delete the strand pairings matrix
  for( int i=0; i<(int)SSEs.size()+2; i++ ) {
    delete Strand_Pairings[i];
  }
  delete Strand_Pairings;

  if( Main_Frag_Selector != NULL ) { delete Main_Frag_Selector; }
  if( Alignment_Builder != NULL ) { delete Alignment_Builder; }
  if( Dist_Measurer != NULL ) { delete Dist_Measurer; }

}



template <class S1, class S2, class Etype>
void SSSS<S1,S2,Etype>::
enumerate (DPMatrix<S1,S2,Etype>& dpm_fwd,
	   AlignmentSet<S1,S2,Etype>& as)
{

  cerr << "enumerate: top" << endl;

  // add fragments and determine connections, number of alis, etc.

  fill_frag_matrix();
  
  //fill_frag_matrix2();

  cerr << "filled frag matrix" << endl;

  // establish which fragments can start an alignment (should be set, but make sure)
  Main_Frag_Selector->find_N_terminal_connections( All_Frags, 0 );

  cerr << "found N terminal cnxns" << endl;

  // give information on fragments distance from native (only prints if native_ali is given)
  Main_Frag_Selector->report_frag_quality( All_Frags );
  Main_Frag_Selector->report_full_sse_frag_set_info( All_Frags );
  cerr << "Final number of alis to search: " << Main_Frag_Selector->get_number_of_alis_to_search( All_Frags ) << endl;

  // print info on fragment set
  cerr << "Fragment Set statistics: " << endl;
  All_Frags->report_statistics();
  cerr << endl;

  // search the fragment set for high-scoring alignments
  build_alignments();

  //  cerr << "Exiting..." << endl;
  //  exit(-1);

  // print some alignment data
  cerr << endl << endl;
  cerr << "Alignment info: " << endl;
  cerr << "Min aligned residues (coverage): " << (int)min_ali_residues << endl;
  cerr << "Number of alignments found: " << Returned_Skel_Alis.size() << endl;
  cerr << "Template size: " << templ_len - 1 << "/" << SSEs.size() << endl;

  if( !Returned_Skel_Alis.empty() ) {
    cerr << "High score: " << Returned_Skel_Alis.front().get_score() << endl;
    cerr << "Low score:  " << Returned_Skel_Alis.back().get_score() << endl;
  }

  // print alignments
  as.clear();
  int ali_id(1);
  
  while( !Returned_Skel_Alis.empty() ) {

    Skel_Ali tmp = Returned_Skel_Alis.front();
    Returned_Skel_Alis.pop_front();
    output_pir_ali( tmp, ali_id, dpm_fwd, as );
    ali_id++;
  }

}



template <class S1, class S2, class Etype>
void SSSS<S1,S2,Etype>::build_alignments()
{
  // setup a Skel_Set variable, find top alis, cluster (if size != 0), ...
  // store returned alis in Returned_Skel_Alis

  Returned_Skel_Alis.clear();  // resulting alignments saved here

  // setup the Skel_Set for searching through the fragment set
  Alignment_Builder = new Skel_Set( (int)min_ali_residues, min_SSE_CO,
				    (int)max_subopt, max_avg_shift * templ_len,
				    All_Frags, Str_Data, Strand_Eval, Dist_Measurer );

  cerr << "enumerate: about to find top skels" << endl;

  Alignment_Builder->find_top_skeletons();

  if( tracking_mode ) { Alignment_Builder->send_culled_alis_to_files(); }

  if( max_avg_shift > 0.0f ) {
    Alignment_Builder->cluster_alignments();
  }

  list<Skel_Ali*> tmp_skels = Alignment_Builder->get_Top_Skels();

  // save skels
  while( !tmp_skels.empty() ) {
    Skel_Ali* tmp = tmp_skels.front();
    tmp_skels.pop_front();
    Returned_Skel_Alis.push_back( *tmp ); // save a permanent version
  }

  delete Alignment_Builder; // clear this Skel_Set variable
  Alignment_Builder = NULL;
}

/*
template <class S1, class S2, class Etype>
void SSSS<S1,S2,Etype>::choose_fragments_for_ali()
{

  fill_frag_matrix(); // add all the fragments and determine connections, number of alis, etc.

  Main_Frag_Selector->print_all_frags();

  Skel_Set All_Skels( SSEs, (int)min_ali_residues, min_SSE_CO,
		      (int)max_subopt, cb_dists, templ_contacts, 
		      max_avg_shift * templ_len,
		      query_seq, templ_seq, TSR_to_Nterm, TSR_to_Cterm, Main_Frag_Selector );

  cerr << "Press any key and hit Enter to see the possible starting fragments." << endl;
  string pause;
  cin >> pause;
  cerr << endl;

  for( int i=0; i<(int)All_Skels.Start_Skels.size(); i++ ) {
    cerr << i+1 << ")" << endl;
    Ali_Frag* first_frag = Main_Frag_Selector->get_frag( All_Skels.Start_Skels[i]->members.front().next_frag );
    first_frag->print( query_seq, templ_seq );
    cerr << endl;
  }

  int choice = get_user_selection();

  Skel_Ali one_skel( All_Skels.Start_Skels[choice-1] );

  cerr << "You have chosen to start with: " << endl;
  one_skel.print( query_seq, templ_seq, (int)min_ali_residues );

  while( true ) {

    Frag_Connection curr = one_skel.members.back(); // the last frag thus far in the ali
   
    int num_next_frags = (int)Main_Frag_Selector->get_frag( curr.next_frag )->next.size();
 
    if( num_next_frags == 0 ) { break; } // only true for C-terminal fragment

    cerr << "Your next choices are: " << endl;

    for( int i=0; i<num_next_frags; i++ ) {
    
      Frag_Connection tmp_fc = Main_Frag_Selector->get_frag( curr.next_frag )->next[i];
      
      cerr << i+1 << ")" << endl;
      Main_Frag_Selector->get_frag( tmp_fc.next_frag )->print( query_seq, templ_seq );
      cerr << endl;
    }

    choice = get_user_selection();
    
    one_skel.add_member( Main_Frag_Selector->get_frag( curr.next_frag )->next[choice-1] ); // elongate the skeleton

    cerr << "You now have: " << endl;
    one_skel.print( query_seq, templ_seq, (int)min_ali_residues );

    cerr << "Press any key and hit Enter to see the fragments that can follow." << endl;
    string pause;
    cin >> pause;
    cerr << endl;

  }

  cerr << "Final skeleton alignment: " << endl;
  one_skel.print( query_seq, templ_seq, (int)min_ali_residues );


}
*/

template <class S1, class S2, class Etype>
vector<Ali_Frag*> SSSS<S1,S2,Etype>::get_new_frags()

{

return (*All_Frags - Old_Frags);

}

template <class S1, class S2, class Etype>
void SSSS<S1,S2,Etype>::fill_frag_matrix()
{
  int i;
  cerr << endl;
  cerr << "Adding fragments until search space exceeds maximum:" << endl;

  float curr_zscore(0);
  
  Old_Frags = *All_Frags; // copy state

  if(nfill>0) for(i=0;i<(int)nfill;i++) {
    if(!Main_Frag_Selector->activate_next_fragment( curr_zscore, max_alis_to_search, All_Frags ) )
      break;
    }
  else while(Main_Frag_Selector->activate_next_fragment( curr_zscore, max_alis_to_search, All_Frags ) ) {
    }

  cerr << "Last frag z-score: " << curr_zscore << endl;
  cerr << endl;
  
}


template <class S1, class S2, class Etype>
void SSSS<S1,S2,Etype>::fill_frag_matrix2()
{

  cerr << endl;
  cerr << "Adding fragments until search space exceeds half maximum:" << endl;

  unsigned long long int half_max = max_alis_to_search / 2;

  cerr << "max_alis_to_search: " << max_alis_to_search << endl;
  cerr << "half_max: " << half_max << endl;

  float curr_zscore(0);

  while( Main_Frag_Selector->activate_next_fragment( curr_zscore, half_max, All_Frags ) )
    { }

  cerr << "Search space size at halfway: "
       << Main_Frag_Selector->get_number_of_alis_to_search( All_Frags ) << endl;

  Main_Frag_Selector->add_in_between_frags2( All_Frags );

  cerr << "Last frag z-score: " << curr_zscore << endl;
  cerr << endl;
  
}


template <class S1, class S2, class Etype>
void SSSS<S1,S2,Etype>::
output_pir_ali( Skel_Ali sa, int ali_id, DPMatrix<S1,S2,Etype>& dpm,AlignmentSet<S1,S2,Etype>& as, ostream& os)
{
  string t_seq("^"), q_seq("^");
  int next_t_res(1), next_q_res(1);

  ali_counter++;

  string key;

  map<string,AlignedPairList<S1,S2> >::iterator map_it = loops.find( key );

  string tmp_templ_str;
  string tmp_query_str;

  int t_loop_beg, q_loop_beg, t_loop_end, q_loop_end;

  os << "#start" << endl;

  // step through all connections
  for( int i=1; i<(int)sa.num_connections(); i++ ) {

    Ali_Frag* prev_af = sa.get_frag( sa.get_connection(i-1).prev_frag );
    Ali_Frag* next_af = sa.get_frag( sa.get_connection(i-1).next_frag );

    // establish frag endpoints
    int t_beg = sa.get_connection(i-1).next_beg_res_idx;
    int t_end = sa.get_connection(i).prev_end_res_idx;
    
    int q_beg = next_af->q( t_beg );
    int q_end = next_af->q( t_end );
    
    // establish loop endpoints
    t_loop_beg = next_t_res;
    q_loop_beg = next_q_res;
    
    t_loop_end = t_beg - 1;
    q_loop_end = q_beg - 1;

    // build string holding loop endpoints
    stringstream key_sstr;
    key_sstr << t_loop_beg-1 << "\t" << q_loop_beg-1 << "\t" << t_loop_end+1 << "\t" << q_loop_end+1;
    key.clear();
    key = key_sstr.str();
    
    // search map with string to see if alignment for this loop has been previously determined
    map_it = loops.find( key );
    
    if( map_it == loops.end() ) {

      if( next_af->sse_id - prev_af->sse_id == 1 ) { // no template SSEs skipped, use nalign to align loop

	// get a sub-dpm for just this region containing the diagonal
	DPMatrix<S1,S2,Etype> tmp_sub_dpm( *dpm.getQuerySequence(),
					   *dpm.getTemplateSequence(),
					   *dpm.getEvaluator(),
					   q_loop_beg-1, t_loop_beg-1, q_loop_end+1, t_loop_end+1,
					   fwd );
	
	// get an optimal sub-alignment through this region
	Optimal_Subali<S1,S2,Etype> opt_subali( q_loop_beg-1, t_loop_beg-1, 
						q_loop_end+1, t_loop_end+1 );
	
	AlignmentSet<S1,S2,Etype> ali_set( tmp_sub_dpm, opt_subali );
	
	loops[key] = ali_set[0];

      }
      else { // loop skips at least one template SSE
	
	AlignedPairList<S1,S2> loop_ali;
	
	loop_ali.append( q_loop_beg-1, t_loop_beg-1 ); // add the pair representing the end of the previous frag
	
	int prev_frag_sse_id = prev_af->sse_id;
	int loop_frag_sse_id = prev_frag_sse_id + 1; // the sse id of the first skipped template SSE
	
	SSE_Data loop_frag = Str_Data->get_SSE_Data( loop_frag_sse_id );

	int num_query_loop_res = q_loop_end - q_loop_beg + 1;
	int num_templ_loop_res = loop_frag.beg_id - t_loop_beg;
	int min_loop_res = min( num_query_loop_res, num_templ_loop_res );

	for( int j=0; j<min_loop_res; j++ ) { // align the minimum number of shared loop residues
	  loop_ali.append( q_loop_beg+j, t_loop_beg+j );
	}

	loop_ali.append( q_loop_end+1, t_loop_end+1 ); // skip to the pair starting the next fragment

	loops[key] = loop_ali;
      }

    }

    map_it = loops.find( key ); // set map_it to point to the proper loop
    
    // get the template and query sequences
    tmp_templ_str = map_it->second.get_templ_string( templ_seq );
    tmp_query_str = map_it->second.get_query_string( query_seq );

    // remove the first and last residues of both sequences (shared with fragment endpoints)
    tmp_templ_str = tmp_templ_str.substr( 1, tmp_templ_str.size() - 2 );
    tmp_query_str = tmp_query_str.substr( 1, tmp_query_str.size() - 2 );

    // insert the loop into the alignment sequences
    t_seq.append( tmp_templ_str );
    q_seq.append( tmp_query_str );

    // insert matched pairs representing fragment
    for( int t=t_beg; t<=t_end; t++ ) {
      t_seq.append( templ_seq.substr( t, 1 ) );
      q_seq.append( query_seq.substr( next_af->q( t ), 1 ) );
    }

    // update next_t_res and next_t_res
    next_t_res = t_end + 1;
    next_q_res = q_end + 1;

  }


  // calculate and append C-terminal loop

  // establish loop endpoints
  t_loop_beg = next_t_res;
  q_loop_beg = next_q_res;

  t_loop_end = templ_seq.size() - 1;
  q_loop_end = query_seq.size() - 1;

  // build string holding loop endpoints
  stringstream key_sstr;
  key_sstr << t_loop_beg-1 << "\t" << q_loop_beg-1 << "\t" << t_loop_end+1 << "\t" << q_loop_end+1;
  key.clear();
  key = key_sstr.str();

  // search map with string to see if alignment for this loop has been previously determined
  map_it = loops.find( key );

  if( map_it == loops.end() ) {

    // get a sub-dpm for just this region containing the diagonal
    DPMatrix<S1,S2,Etype> tmp_sub_dpm( *dpm.getQuerySequence(),
				       *dpm.getTemplateSequence(),
				       *dpm.getEvaluator(),
				       q_loop_beg-1, t_loop_beg-1, q_loop_end, t_loop_end,
				       fwd );
    
    // get an optimal sub-alignment through this region
    Optimal_Subali<S1,S2,Etype> opt_subali( q_loop_beg-1, t_loop_beg-1, 
					    q_loop_end, t_loop_end );
      
    AlignmentSet<S1,S2,Etype> ali_set( tmp_sub_dpm, opt_subali );
    
    loops[key] = ali_set[0];
    
    map_it = loops.find( key ); // set it to point to new loop
  }

  // get the template and query sequences
  tmp_templ_str = map_it->second.get_templ_string( templ_seq );
  tmp_query_str = map_it->second.get_query_string( query_seq );
  
  // remove the first and last residues of both sequences (shared with fragment endpoints)
  tmp_templ_str = tmp_templ_str.substr( 1, tmp_templ_str.size() - 2 );
  tmp_query_str = tmp_query_str.substr( 1, tmp_query_str.size() - 2 );
  
  t_seq.append( tmp_templ_str );   // insert the loop into the alignment sequences
  q_seq.append( tmp_query_str );
  
  t_seq.append( "*" );   // append final '*' character for proper PIR files
  q_seq.append( "*" );

  // print PIR-formatted files to cout
  int max_line_length(60);

  os << ">P1;templ" << endl;
  os << "structure:" << endl;

  int last_break( 0 );
  int next_break( max_line_length );

  // print all but the last line of the template
  while( (int)t_seq.size() > next_break ) {
    os << t_seq.substr( last_break, max_line_length ) << endl;
    last_break = next_break;
    next_break += max_line_length;
  }

  // print the last line of the template
  os << t_seq.substr( last_break ) << endl;
  last_break = next_break;
  next_break += max_line_length;

  os << ">P1;query" << endl;
  os << "sequence:mdl_" << ali_id << endl;

  last_break = 0;
  next_break = max_line_length;

  // print all but the last line of the query
  while( (int)q_seq.size() > next_break ) {
    os << q_seq.substr( last_break, max_line_length ) << endl;
    last_break = next_break;
    next_break += max_line_length;
  }

  // print the last line of the query
  os << q_seq.substr( last_break ) << endl;
  last_break = next_break;
  next_break += max_line_length;

  //  os << endl;
  os << "#end" << endl;

  SingleAlignment al;
  int t_idx(1),q_idx(1);
  for(int i=1;i<(int)t_seq.size();i++) {
    char tc=t_seq[i];
    char qc=q_seq[i];
    if(tc=='-') {
      if(qc=='-') continue;
      else q_idx++;
      continue;
      }
    if(qc=='-') {
      if(tc=='-') continue;
      else t_idx++;
      continue;
      }
    al.append(q_idx++,t_idx++);
    }
  as.push_back(al);

}

template <class S1, class S2, class Etype>
void SSSS<S1,S2,Etype>::
setup_data_structures()
{

  // load sequence strings and lengths
  Str_Data->load_seq_lengths( (int)templ_seq.size(), (int)query_seq.size() );
  Str_Data->load_seq_strings( templ_seq, query_seq );

  // setup the similarity matrix, sims**
  sims = new float* [ query_seq.size() ];

  for( int i=0; i<(int)query_seq.size(); i++ ) {
    sims[i] = new float[ templ_seq.size() ];
  }

  for( int i=0; i<(int)query_seq.size(); i++ ) {
    for( int j=0; j<(int)templ_seq.size(); j++ ) {
      sims[i][j] = dpm->getSim(i,j);
    }
  }

  Str_Data->load_sims( sims );


  // setup the CB_distance matrix, cb_dists**
  cb_dists = new float* [ templ_seq.size() ];

  for( int i=0; i<(int)templ_seq.size(); i++ ) {
    cb_dists[i] = new float [ templ_seq.size() ];
  }

  for( int i=0; i<(int)templ_seq.size(); i++ ) {
    for( int j=0; j<(int)templ_seq.size(); j++ ) {
      cb_dists[i][j] = ( templ->at(i)->rdata.cb - templ->at(j)->rdata.cb ).norm();
    }
  }

  Str_Data->load_cb_dists( cb_dists );


  // setup the contact matrix, templ_contacts**
  templ_contacts = new bool* [ templ_seq.size() ];

  for( int i=0; i<(int)templ_seq.size(); i++ ) {
    templ_contacts[i] = new bool [ templ_seq.size() ];
  }

  for( int i=0; i<(int)templ_seq.size(); i++ ) {
    for( int j=0; j<(int)templ_seq.size(); j++ ) {
      templ_contacts[i][j] = ( cb_dists[i][j] < max_contact_dist );
    }
  }

  for( int i=0; i<(int)templ_seq.size(); i++ ) { // initialize contacts matrix
    templ_contacts[0][i] = false;
    templ_contacts[i][0] = false;
    templ_contacts[templ_seq.size()-1][i] = false;
    templ_contacts[i][templ_seq.size()-1] = false;
  }

  Str_Data->load_contacts( templ_contacts );


  // setup the query predicted loops array
  query_predicted_loops = new bool [ query_len+1 ];

  for( int i=0; i<query_len+1; i++ ) {

    // identify confidently-predicted loop residues
    query_predicted_loops[i] = ( ( query->at(i)->sse_values[2] == 1.f ) && 
				 ( query->at(i)->sse_confid > .85f ) );
  }

  Str_Data->load_query_predicted_loops( query_predicted_loops );


  // setup the template SSE data structures
  int templ_res_idx(0);
  int sse_id(1);
  int min_sse_len(3);
  while( templ_res_idx < (int)templ_seq.size() ) {

    // advance until next SSE is found
    while( ( templ_res_idx < (int)templ_seq.size() ) &&
	   ( templ->at(templ_res_idx )->rdata.isse == -1 ) ) {
      templ_res_idx++;
    }

    if( templ_res_idx >= (int)templ_seq.size() ) { break; }


    // templ_res_idx is now at the first residue of the next SSE
    SSE_Data tmp;
    tmp.beg_id = templ_res_idx; // save the SSE start
    tmp.ss_type = templ->at(templ_res_idx)->rdata.sse_type; // save the SSE type

    // advance to one past the end of the current SSE
    while( ( templ_res_idx < (int)templ_seq.size() ) &&
	   ( templ->at(templ_res_idx )->rdata.isse != -1 ) ) {
      templ_res_idx++;
    }

    tmp.end_id = templ_res_idx - 1; // save the SSE end

    // skip short SSEs
    if( ( tmp.end_id - tmp.beg_id ) + 1 < min_sse_len ) { continue; }

    tmp.sse_id = sse_id;
    sse_id++;

    SSEs.push_back( tmp );
  }

  Str_Data->load_SSE_Data( SSEs );


  // setup the max TSR arrays
  TSR_to_Nterm = new int [ templ_seq.size() ];
  TSR_to_Cterm = new int [ templ_seq.size() ];

  int idx = 0;
  while( (int)idx < SSEs[0].beg_id ) { TSR_to_Nterm[ idx ] = 0; idx++; }

  for( int i=0; i<(int)SSEs.size()-1; i++ ) {
    for( idx=SSEs[i].beg_id; (int)idx<=SSEs[i].end_id; idx++ ) {
      TSR_to_Nterm[ idx ] = TSR_to_Nterm[ idx-1 ] + 1;
    }

    while( (int)idx < SSEs[i+1].beg_id ) {
      TSR_to_Nterm[ idx ] = TSR_to_Nterm[ idx-1 ];
      idx++;
    }
  }

  for( idx=SSEs.back().beg_id; (int)idx<=SSEs.back().end_id; idx++ ) {
    TSR_to_Nterm[ idx ] = TSR_to_Nterm[ idx-1 ] + 1;
  }

  while( idx < (int)templ_seq.size() ) { 
    TSR_to_Nterm[ idx ] = TSR_to_Nterm[ SSEs.back().end_id ];
    idx++;
  }

  int total_TSR = TSR_to_Nterm[ templ_seq.size() - 1 ];

  for( int i=0; i<(int)templ_seq.size(); i++ ) {
    TSR_to_Cterm[i] = total_TSR - TSR_to_Nterm[i];
  }

  for( int i=0; i<(int)SSEs.size(); i++ ) {
    for( idx=SSEs[i].beg_id; (int)idx<=SSEs[i].end_id; idx++ ) {
      TSR_to_Cterm[ idx ] = ( total_TSR + 1 ) - TSR_to_Nterm[ idx ];
    }
  }

  Str_Data->load_TSRs( TSR_to_Nterm, TSR_to_Cterm );


  // set up array (Strand_Pairings**) to record beta-strand pairings based on H-bond connectivities
  Strand_Pairings = new bool* [ SSEs.size()+2 ];

  for( int i=0; i<(int)SSEs.size()+2; i++ ) { // create space for Strand_Pairings
    Strand_Pairings[i] = new bool [ i+1 ];
  }

  for( int i=0; i<(int)SSEs.size()+2; i++ ) { // intialize all values to False
    for( int j=0; j<=i; j++ ) {
      Strand_Pairings[i][j] = false;
    }
  }

  for( int i=1; i<(int)SSEs.size(); i++ ) { // step through SSEs and record pairings

    int ss1_beg = SSEs[i].beg_id;
    int ss1_end = SSEs[i].end_id;

    for( int j=0; j<i; j++ ) {

      int ss2_beg = SSEs[j].beg_id;
      int ss2_end = SSEs[j].end_id;

      int total_hbonds(0);

      for( int m=ss1_beg; m<=ss1_end; m++ ) {
	for( int n=ss2_beg; n<=ss2_end; n++ ) {

	  if( templ->get_backbone_HB_contact( m, n ) ) {
	    total_hbonds++;
	  }
	}
      }

      if( total_hbonds > 0 ) {
	Strand_Pairings[ SSEs[i].sse_id ][ SSEs[j].sse_id ] = true;
      }

    }

  }

}


template <class S1, class S2, class Etype>
vector<SSE_Data> SSSS<S1,S2,Etype>::
find_template_sses( int size ) 
{
  vector<SSE_Data> res; // to return

  //  int min_sse_len = 2;
  int min_sse_len = 3;

  SSE_Data tmp_sse;
  tmp_sse.beg_id = -1;
  tmp_sse.end_id = -1;

  for( int i=1; i<size-1; i++ ) {

    // check for first SSE starting at first residue
    if( (i==1) && ( templ->at(i)->rdata.sse_type != 331 ) ) {
      tmp_sse.beg_id = i;
    }

    // check for typical start of SSE
    if( (i>1) && ( templ->at(i)->rdata.sse_type != 331 ) &&
	!( templ->at(i-1)->rdata.sse_type != 331 ) ) {
      tmp_sse.beg_id = i;
    }

    // check for last SSE ending at second to last residue
    if( (i==size-2) && ( templ->at(i)->rdata.sse_type != 331 ) ) {
      tmp_sse.end_id = i;
    }

    // check for typical end of SSE
    if( (i<size-2) && ( templ->at(i)->rdata.sse_type != 331 ) && 
	!( templ->at(i+1)->rdata.sse_type != 331 ) ) {
      tmp_sse.end_id = i;
    }

    if( (tmp_sse.beg_id != -1) && (tmp_sse.end_id != -1) && 
	( tmp_sse.beg_id <= tmp_sse.end_id ) &&
	( tmp_sse.end_id - tmp_sse.beg_id + 1 >= min_sse_len ) ) {

      res.push_back( tmp_sse );
      tmp_sse.beg_id = -1;
      tmp_sse.end_id = -1;
    }

  }

  return res;

}

template <class S1, class S2, class Etype>
vector<Ali_Frag*> SSSS<S1,S2,Etype>::GetAlignedFrags()
{
  int i;

  if(Main_Frag_Selector==NULL) throw(string("Frag_Matrix undefined"));
  if(All_Frags==NULL) throw(string("Frag_Set undefined"));

  return All_Frags->export_all_frags();
  
}


template <class S1, class S2, class Etype>
int SSSS<S1,S2,Etype>::get_user_selection()
{
  int res;
  cerr << "Enter an integer:  ";
  cin >> res;
  return res;
}

template <class S1, class S2, class Etype>
void SSSS<S1,S2,Etype>::pause()
{
  string pause;
  cerr << "Hit a key and press Enter to continue." << endl;
  cin >> pause;
}


#endif  //_HMAP2_SSSS
