#ifndef _UPGMA_TREE
#define _UPGMA_TREE

#include <iostream>
#include <string>
#include <vector>
#include <list>

#include "ssss_shared_defs.h"
#include "UPGMA_Tree.h"

using namespace std;

class UPGMA_Tree {

 private:

  // DATA
  UPGMA_Tree* l_node; // left child
  UPGMA_Tree* r_node; // right child

  float l_dist, r_dist;
  int index;
  int weight;

  float avg_leaf_dist;

  // FUNCTIONS

 public:

  UPGMA_Tree( UPGMA_Tree*, UPGMA_Tree*, float, float, int );
  ~UPGMA_Tree();

  // DATA


  // FUNCTIONS

  // general
  void set_avg_leaf_dist();
  bool it_is_terminal_node();
  void set_weight();
  //  void set_parent( UPGMA_Tree* );
  vector<UPGMA_Tree*> get_leaves();

  // access
  inline int get_index() const { return index; }
  inline float get_avg_leaf_dist() const { return avg_leaf_dist; }
  inline UPGMA_Tree* get_l_node() const { return l_node; }
  inline UPGMA_Tree* get_r_node() const { return r_node; }
  inline int get_weight() const { return weight; }

  // debugging
  void print_node();
  void print_tree();
};

#endif  //_SKELSET
