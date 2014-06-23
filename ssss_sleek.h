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

#include "frag_matrix.h"
#include "skel_set.h"
#include "skel_ali.h"
#include "ali_frag.h"

struct SSE_Data {

  // holds basic info on template SSEs

  int sse_id;
  int beg_id;
  int end_id;
};

using namespace std;

template <class S1, class S2, class Etype>
class SSSS : public Enumerator<S1,S2,Etype> {
  
  // SSSS = Sample Shifts in Secondary Structures
  // This class is a sub-optimal alignment enumerator.

  // It has two phases.  First, it scans through a DP matrix
  // looking for all possible alignments to each template secondary
  // structure element that would lead to decent (above a threshold)
  // global alignments.  The results is a list of alignment fragments
  // for each template SSE.  The second phase uses a different scoring
  // function to recombine these fragments into full alignments and
  // returns the top N alignments.


  // The structs below are a bit confusing.  An 'SSE_Frag' is made
  // for each alignment fragment found for a particular SSE.

 public:
  
  typedef AlignedPairList<S1,S2> SingleAlignment;
  typedef AlignedPair<S1,S2>     SinglePair;

  SSSS (const NOaliParams&, const Evaluator<S1,S2,Etype>&, 
	int, int, float, float, int, int, int, float ) ;

  ~SSSS();

  inline int estimateSize () const { return params->number_suboptimal; } ;
  void enumerate (DPMatrix<S1,S2,Etype>& dpm_fwd,
		  AlignmentSet<S1,S2,Etype>& as);

 private:
  
  // reads flags file to get locations of template sses
  //  vector<SSE_Data> find_sses_from_flags( int );
  vector<SSE_Data> find_template_sses( int );

  // takes in a skeleton_ali and returns a full alignment for final printing
  SingleAlignment convert_skel_to_APL( Skel_Ali* );
  void output_pir_ali( Skel_Ali*, int );

  const NOaliParams* params;

  const S1* query;
  const S2* templ;

  string query_seq;
  string templ_seq;
  int query_len, templ_len; // sequence lengths

  const Evaluator<S1,S2,Etype> evaluator;
  int mode;
  unsigned int max_skip; // max # of sses that can be skipped over when forming an alignment
  unsigned long long int max_alis_to_search;
  unsigned int max_subopt; // # of alignments that user wants returned from 'sample_shifts'
  unsigned int num_clusters; // # of clusters to group alignments into
  unsigned int ali_mode; // 0 = align fragments to entire template SSE
                           // 1 = align to part of template SSE

  float min_coverage;
  float min_SSE_CO;
  float max_avg_shift;  // when clustering, the maximum number of residues one alignment
                        // of a cluster can be from another member

  float shorter_seq_len;

  float min_ali_residues;
  float max_contact_dist;

  bool** contact;
  float templ_SSE_CO;

  float min_ali_separation;

  long int skels_created, skels_deleted;

  AlignedPairList<AASequence,AASequence> native_ali;
};

//constructor
template <class S1, class S2, class Etype>
SSSS<S1,S2,Etype>::
SSSS (const NOaliParams& p, const Evaluator<S1,S2,Etype>& eval,
      int num_alis_kept, int max_alis,
      float min_cov, float min_CO, int skip_sses,
      int num_c, int ali_how, float shift ) 
  : params(&p), 
  evaluator(eval),
  max_skip (skip_sses),
  max_alis_to_search(max_alis),
  max_subopt(num_alis_kept),
  num_clusters( num_c ),
  ali_mode( ali_how ),
  min_coverage( min_cov ),
  min_SSE_CO( min_CO ),
  max_avg_shift( shift )
{

    // setup initial values
    max_contact_dist = 8.f;
    //    min_SSE_CO = 0.60f;
    skels_created = 0;
    skels_deleted = 0;
}

//destructor
template <class S1, class S2, class Etype>
SSSS<S1,S2,Etype>::
~SSSS()
{}



template <class S1, class S2, class Etype>
void SSSS<S1,S2,Etype>::
enumerate (DPMatrix<S1,S2,Etype>& dpm_fwd,
	   AlignmentSet<S1,S2,Etype>& as)
{
  cerr << "ssss: enumerate: top" << endl;

  query_len = dpm_fwd.getQuerySize() - 1;
  templ_len = dpm_fwd.getTemplateSize() - 1;

  min_ali_separation = 0.1f * templ_len;

  shorter_seq_len = min( query_len-1, templ_len-1 );
  min_ali_residues = min_coverage * shorter_seq_len;

  // fill in 'query_seq' and 'templ_seq' public string variabls
  query_seq = *dpm_fwd.getQuerySequence()->getString();
  templ_seq = *dpm_fwd.getTemplateSequence()->getString();

  // fill 'query' and 'templ' Sequence variables
  query = dpm_fwd.getQuerySequence();
  templ = dpm_fwd.getTemplateSequence();

  // setup the similarity matrix
  float** sims = new float* [ query_seq.size() ];

  for( unsigned int i=0; i<query_seq.size(); i++ ) {
    sims[i] = new float[ templ_seq.size() ];
  }

  for( unsigned int i=0; i<query_seq.size(); i++ ) {
    for( unsigned int j=0; j<templ_seq.size(); j++ ) {
      sims[i][j] = dpm_fwd.getSim(i,j);
    }
  }

  // setup the CA_distance matrix
  float** ca_dists = new float* [ templ_seq.size() ];

  for( unsigned int i=0; i<templ_seq.size(); i++ ) {
    ca_dists[i] = new float [ templ_seq.size() ];
  }

  for( unsigned int i=0; i<templ_seq.size(); i++ ) {
    for( unsigned int j=0; j<templ_seq.size(); j++ ) {
      ca_dists[i][j] = ( templ->at(i)->rdata.ca - templ->at(j)->rdata.ca ).norm();
    }
  }

  // setup the contact matrix
  bool** templ_contacts = new bool* [ templ_seq.size() ];

  for( unsigned int i=0; i<templ_seq.size(); i++ ) {
    templ_contacts[i] = new bool [ templ_seq.size() ];
  }

  for( unsigned int i=0; i<templ_seq.size(); i++ ) {
    for( unsigned int j=0; j<templ_seq.size(); j++ ) {
      templ_contacts[i][j] = ( ca_dists[i][j] < max_contact_dist );
    }
  }

  // get the template sse boundaries
  vector<SSE_Data> SSEs = find_template_sses( templ_seq.size() );
  vector<vector<int> > sse_locs;

  for( unsigned int i=0; i<SSEs.size(); i++ ) {
    vector<int> tmp;
    tmp.push_back( SSEs[i].beg_id );
    tmp.push_back( SSEs[i].end_id );
    sse_locs.push_back( tmp );
  }


  cerr << "ssss: about to create Frag_Matrix" << endl;

  // setup the main Frag_Matrix
  Frag_Matrix Main( query_seq, templ_seq, sse_locs, (int)min_ali_residues, sims, ca_dists, ali_mode );
  //  Main.print_template_frags();

  cerr << "ssss:  created Frag Matrix" << endl;

  Main.fill_matrix_columns( 10 ); // 10 frags per sse to start with

  cerr << "ssss:  filled matrix columns" << endl;

  cerr << "ssss: Initial frag matrix:" << endl;
  Main.print_Active_Frags_matrix();

  cerr << "ssss: Initial INACTIVE frag matrix:" << endl;
  Main.print_Inactive_Frags_matrix();

  cerr << "ssss: about to find fragment connections" << endl;
  Main.find_fragment_connections();

  cerr << "ssss:  passed find_fragment_connection" << endl;

  //  string pause;

  Main.count_frag_children();

  cerr << "ssss:  past count_frag_children" << endl;

  while( true ) {
    unsigned long long int num_alis_to_search = Main.get_number_of_alis_to_search();

    cerr << "ssss: alis to search (before adding frag): " << num_alis_to_search << endl;

    if( num_alis_to_search > max_alis_to_search ) { break; }

    Main.fill_in_Frag_Matrix();

    cerr << "ssss: alis to search: " << Main.get_number_of_alis_to_search() << endl;
    Main.report_Matrix_statistics();
  }

  cerr << "ssss: Final frag matrix:" << endl;
  Main.print_Active_Frags_matrix();

  cerr << "ssss: about to find N-terminal connections" << endl;
  Main.find_N_terminal_connections( 0 );

  cerr << "ssss: Final number of alis to search: " << Main.get_number_of_alis_to_search() << endl;

  Skel_Set All_Skels( sse_locs, (int)min_ali_residues, min_SSE_CO,
		      (int)max_subopt, ca_dists, templ_contacts, 
		      max_avg_shift * templ_len,
		      query_seq, templ_seq, &Main );

  cerr << "ssss: about to find top skeletons" << endl;

  All_Skels.find_top_skeletons();

  cerr << "ssss: found top skeletons" << endl;

  if( max_avg_shift > 0.0f ) {
    All_Skels.cluster_alignments();
  }

  cerr << "ssss: clustered top skeletons" << endl;

  list<Skel_Ali*> Top_Skels_Found = All_Skels.get_Top_Skels();

  cerr << "ssss: Number of alignments found: " << Top_Skels_Found.size() << endl;

  as.clear();

  int ali_id(0);

  while( !Top_Skels_Found.empty() ) {

    Skel_Ali* tmp = Top_Skels_Found.front();

    //    cerr << "Skel to be converted:" << endl;
    //    tmp->print( query_seq, templ_seq, 25 );

    Top_Skels_Found.pop_front();

    output_pir_ali( tmp, ali_id );

    ali_id++;
  }

  //  cerr << "Number of alignments found: " << as.size() << endl;
  cerr << "Template size: " << templ_len - 1 << "/" << sse_locs.size() << endl;
}


template <class S1, class S2, class Etype>
void SSSS<S1,S2,Etype>::
output_pir_ali( Skel_Ali* sa, int ali_id )
{
  string t_seq, q_seq;
  int next_t_res(1), next_q_res(1);

  string gap_char = "-";
  string end_char = "*";

  cout << "#start" << endl;

  // step through all connections
  for( int i=1; i<(int)sa->members.size(); i++ ) {

    Frag_ID f_id = sa->members[i].prev_frag;
    Ali_Frag* af = sa->src->get_frag( f_id );

    int t_beg = sa->members[i-1].next_beg_res_idx;
    int t_end = sa->members[i].prev_end_res_idx;

    int q_beg = af->q( t_beg );
    int q_end = af->q( t_end );

    // insert zigzagged gaps representing loop (prior to fragment) in template and query


    int t_loop_beg = next_t_res;
    //    int q_loop_beg = next_q_res;

    int t_loop_end = t_beg - 1;
    int q_loop_end = q_beg - 1;


    if( t_loop_beg == t_loop_end ) { // template loop has only 1 residue

      // insert the one template loop residue
      t_seq.append( templ_seq.substr( next_t_res, 1 ) );
      q_seq.append( gap_char );
      next_t_res++;

      // insert all the remaining query loop residues
      while( next_q_res <= q_loop_end ) {
	t_seq.append( gap_char );
	q_seq.append( query_seq.substr( next_q_res, 1 ) );
	next_q_res++;
      }

    }
    else { // template loop is at least 2 residues long

      int t_loop_half_len;

      if( ( ( t_loop_end - t_loop_beg + 1 ) % 2 ) == 0 ) { // templ loop length is even
	t_loop_half_len = ( t_loop_end - t_loop_beg + 1 ) / 2;
      }
      else { // templ loop length is odd
	t_loop_half_len = ( t_loop_end - t_loop_beg ) / 2; // smaller half ( ie. (3)+4=7 )
      }

      // insert the first half of the template loop residues
      while( next_t_res<t_loop_beg+t_loop_half_len ) {
 	t_seq.append( templ_seq.substr( next_t_res, 1 ) );
	q_seq.append( gap_char );
	next_t_res++;
      }

      // insert all the remaining query loop residues
      while( next_q_res <= q_loop_end ) {
	t_seq.append( gap_char );
	q_seq.append( query_seq.substr( next_q_res, 1 ) );
	next_q_res++;
      }
      
      // insert the second half of the template loop residues
      while( next_t_res<=t_loop_end ) {
	t_seq.append( templ_seq.substr( next_t_res, 1 ) );
	q_seq.append( gap_char );
	next_t_res++;
      }

    }


    /*
    while( next_t_res < t_beg ) {
      t_seq.append( templ_seq.substr( next_t_res, 1 ) );
      q_seq.append( gap_char );
      next_t_res++;
    }

    while( next_q_res < q_beg ) {
      t_seq.append( gap_char );
      q_seq.append( query_seq.substr( next_q_res, 1 ) );
      next_q_res++;
    }
    */



    // insert matched pairs representing fragment
    for( int t=t_beg; t<=t_end; t++ ) {
      t_seq.append( templ_seq.substr( t, 1 ) );
      q_seq.append( query_seq.substr( af->q( t ), 1 ) );
    }

    // update next_t_res and next_t_res
    next_t_res = t_end + 1;
    next_q_res = q_end + 1;

  }

  // append last loop characters
  while( next_t_res < (int)( templ_seq.size() - 1 ) ) {
    t_seq.append( templ_seq.substr( next_t_res, 1 ) );
    q_seq.append( gap_char );
    next_t_res++;
  }

  while( next_q_res < (int)( query_seq.size() - 1 ) ) {
    t_seq.append( gap_char );
    q_seq.append( query_seq.substr( next_q_res, 1 ) );
    next_q_res++;
  }

  // append final '*' character for proper PIR files
  t_seq.append( end_char );
  q_seq.append( end_char );


  // print PIR-formatted files to cout
  int max_line_length(60);

  cout << ">P1;templ" << endl;
  cout << "structure:" << endl;

  int last_break( 0 );
  int next_break( max_line_length );

  // print all but the last line of the template
  while( (int)t_seq.size() > next_break ) {
    cout << t_seq.substr( last_break, max_line_length ) << endl;
    last_break = next_break;
    next_break += max_line_length;
  }

  // print the last line of the template
  cout << t_seq.substr( last_break ) << endl;
  last_break = next_break;
  next_break += max_line_length;

  cout << ">P1;query" << endl;
  cout << "sequence:mdl_" << ali_id << endl;

  last_break = 0;
  next_break = max_line_length;

  // print all but the last line of the query
  while( (int)q_seq.size() > next_break ) {
    cout << q_seq.substr( last_break, max_line_length ) << endl;
    last_break = next_break;
    next_break += max_line_length;
  }

  // print the last line of the query
  cout << q_seq.substr( last_break ) << endl;
  last_break = next_break;
  next_break += max_line_length;

  //  cout << endl;
  cout << "#end" << endl;

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

#endif  //_HMAP2_SSSS
