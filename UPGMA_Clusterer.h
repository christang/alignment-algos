#ifndef _UPGMA_CLUSTERER
#define _UPGMA_CLUSTERER

#include <iostream>
#include <string>
#include <vector>
#include <list>

#include "ssss_shared_defs.h"
#include "UPGMA_Tree.h"

using namespace std;

//class UPGMA_Tree;

class UPGMA_Clusterer {

 private:

  // DATA

  float** distance_orig; // original distance matrix
  float** distance_curr; // current matrix, reduced in size by one after each step

  vector<UPGMA_Tree> nodes; // set of nodes to cluster
  UPGMA_Tree* root; // base of tree, access point

  vector<vector<UPGMA_Tree*> > clusters;

  int num_points_orig, num_points_curr; // size of arrays

  int next_virtual_index;
  int next_node_index;

  int* w; // track # of real nodes each cluster represents
  int* node_index; // track indices of each node as they swap, etc.

  // FUNCTIONS

 public:

  UPGMA_Clusterer( float**, int );
  ~UPGMA_Clusterer();

  // DATA

  // FUNCTIONS
  void cluster();

  // general
  void find_clusters_under_threshold( float );
  void find_clusters_under_threshold( UPGMA_Tree*, float );
  void print_clusters_under_threshold();
  float find_closest_pair( int&, int& );
  void combine_nodes( int, int );
  void update_distance_matrix( int, int );
  void swap_cols( int, int );
  void delete_distance_matrix( float**, int );

  // access
  inline float get_dist_curr( int i, int j ) const { return( i>j ? distance_curr[i][j] :
							     distance_curr[j][i] ); }
  inline int get_num_clusters() const { return( clusters.size() ); }
  inline int get_num_members( int i ) const { return( clusters[i].size() ); }
  inline int get_member_index( int i, int j ) const { return clusters[i][j]->get_index(); }

  // debugging
  void print_distance_matrix( float**, int );

};

#endif  //_SKELSET
