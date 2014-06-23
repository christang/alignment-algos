#ifndef _ALIDIST
#define _ALIDIST

#include <iostream>
#include <fstream>
#include <vector>

#include "ssss_shared_defs.h"

using namespace std;

class Ali_Dist {

 private:

// DATA //
  // the alignments from which the area will be calculated
  vector<Res_Pair> main_ali;
  vector<Res_Pair> test_ali;

  // misc
  int templ_length; // the number of residues in the template sequence
  int query_length; // the number of residues in the query sequence
 public:
  vector<vector<float> > batch_dists;
 private:

// FUNCTIONS //
  // convert different alignment formats to vector<Res_Pair>
  void convert_strings_to_VRP( const string&, const string&, vector<Res_Pair>& );

  int get_sequence_length( const string& );
  int get_sequence_length( const vector<Res_Pair>&, char );

  bool extract_next_ali( ifstream&, string&, string& );
  void format_string_ends( string& );

  // area calculation
  int get_relative_position( float, float, const vector<Res_Pair>&);
  void insert_intersections( vector<Res_Pair>&, vector<Res_Pair>& );
  void insert_matching_points( vector<Res_Pair>&, vector<Res_Pair>& );
  float calculate_area_between_VRPs( vector<Res_Pair>&, vector<Res_Pair>& );

  // mutual coverage calculation
  float calc_templ_mutual_coverage();
  float calc_query_mutual_coverage();

 public:

// DATA //

// FUNCTIONS //
  // constructor/destructor
  Ali_Dist();
  ~Ali_Dist();

// GENERAL
  // load alignments from different source formats
  void load_main( string ); // load a fasta file
  void load_main( const vector<Res_Pair>& );
  void load_test( string );
  void load_test( const vector<Res_Pair>& );
  void load_PIR_file( string );
  void batch_compare_to_main_ali( string fn );

  float get_area_between_main_and_test();
  float get_dist_between_main_and_test();

  vector<Res_Pair> get_local_native_ali( int, int );
  float get_local_qt_shift( int, int );

  void print_batch_dists( ostream& = cout );

// OPERATORS



// ACCESS

// DEBUGGING
  void print_vrp( const vector<Res_Pair>& );
  void pause();
};

#endif  //_ALIDIST
