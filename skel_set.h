#ifndef _SKELSET
#define _SKELSET

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include "math.h"

#include "ssss_shared_defs.h"

#include "ali_dist.h"
#include "ali_frag.h"
#include "frag_set.h"
#include "ali_str_info.h"
#include "ali_strand_eval.h"
#include "skel_ali.h"
#include "UPGMA_Clusterer.h"

using namespace std;

class Skel_Set {

 private:

  list<Skel_Ali*> Top_Skels;

  int min_aligned_residues;
  float min_SSE_CO;
  float min_SSE_CO_fraction;
  int max_alis, max_bad_alis;

  bool** contacts;
  float** inter_ali_area;

  int* TSR_to_Nterm;
  int* TSR_to_Cterm;

  list<Skel_Ali*> Low_Coverage, Low_SSE_CO, Bad_Strands, Low_Score;

  ofstream ofs_cov, ofs_sse_co, ofs_strand_rules, ofs_score;

  string query_seq, templ_seq;

  bool tracking_mode;

  Ali_Dist* Measurer;
  UPGMA_Clusterer* clusterer;
  Frag_Set* Frags;
  Ali_Str_Info* Str_Data;

  Alignment_Strand_Evaluator* Strand_Eval;

  float max_cluster_size;

  unsigned long int skels_created, skels_deleted;

  unsigned long int num_culled_by_coverage, num_culled_by_contact_order;
  unsigned long int num_culled_by_strand_rules, num_culled_by_score;

  Skel_Ali* top_constrained_skel;

  // FUNCTIONS
  void pre_empt_low_coverage( Skel_Ali* );
  void handle_completed_skel( Skel_Ali* );
  void handle_completed_constrained_skel( Skel_Ali* );
  int find_next_post( Skel_Ali*, Skel_Ali*, int );
  bool passes_all_filters( Skel_Ali*, int& );
  void sort_top_skels( Skel_Ali* );
  void sort_culled_skels( Skel_Ali*, list<Skel_Ali*>&, int );

 public:

  // data
  vector<Skel_Ali*> Start_Skels;

  // FUNCTIONS
  Skel_Set() { }
  Skel_Set( int, float,
	    int, float,
	    Frag_Set*, Ali_Str_Info*, Alignment_Strand_Evaluator*, Ali_Dist* );
  ~Skel_Set();

  // general
  void find_top_skeletons();
  void find_top_constrained_skel( Skel_Ali* );

  void grow_skel( Skel_Ali* );
  void grow_constrained_skel( Skel_Ali*, Skel_Ali*, int );

  float find_template_SSE_CO();

  void handle_culled_skel_ali( Skel_Ali*, int );
  void send_culled_alis_to_files();



  void cluster_alignments();
  void get_exact_inter_ali_areas( const vector<Skel_Ali*>& );


  // access
  inline list<Skel_Ali*> get_Top_Skels() const { return Top_Skels; };
  Ali_Frag* get_frag( Frag_ID f ) const;
  Ali_Frag* get_frag( int s, int f ) const;

  // debugging
  void print_Top_Skels();
  void pause();
};

#endif  //_SKELSET
