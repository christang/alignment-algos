#include "UPGMA_Tree.h"

using namespace std;

/*********************************************
* implementation of UPGMA_Tree class *
*********************************************/

// constructor 1
UPGMA_Tree::UPGMA_Tree( UPGMA_Tree* left, UPGMA_Tree* right, float ld, float rd, int id )
{
  //  cerr << "UPGMA_Tree constructor: top" << endl;

  index = id;
  l_node = left;
  r_node = right;

  l_dist = ld;
  r_dist = rd;

  /*
  if( l_node != (UPGMA_Tree*)NULL ) {
    cerr << "UPGMA_Tree constructor: l_node tree: " << endl;
    l_node->print_tree();
  }
  else {
    cerr << "UPGMA_Tree constructor: l_node is NULL" << endl;
  }

  if( r_node != (UPGMA_Tree*)NULL ) {
    cerr << "UPGMA_Tree constructor: r_node tree: " << endl;
    r_node->print_tree();
  }
  else {
    cerr << "UPGMA_Tree constructor: r_node is NULL" << endl;
  }
  */

  //  cerr << "UPGMA_Tree constructor: about to set avg leaf dist" << endl;
  set_avg_leaf_dist();
  //  cerr << "UPGMA_Tree constructor: avg leaf dist: " << get_avg_leaf_dist() << endl;

  //  cerr << "UPGMA_Tree constructor: about to set weight" << endl;
  set_weight();
  //  cerr << "UPGMA_Tree constructor: weight: " << get_weight() << endl;

  //  cerr << "UPGMA_Tree constructor: end" << endl;
}

UPGMA_Tree::~UPGMA_Tree()
{}


// GENERAL
void UPGMA_Tree::set_avg_leaf_dist()
{
  // Purpose: determine the average distance between the node and each leaf connected to it.

  //  cerr << "set_avg_leaf_dist: top" << endl;

  if( it_is_terminal_node() ) {
    //    cerr << "set_avg_leaf_dist: inside if" << endl;
    avg_leaf_dist = 0.f;
  }
  else {
    //    cerr << "set_avg_leaf_dist: inside else" << endl;
    avg_leaf_dist = ( ( l_node->get_weight() * ( l_dist + l_node->get_avg_leaf_dist() ) ) + 
		      ( r_node->get_weight() * ( r_dist + r_node->get_avg_leaf_dist() ) ) ) / 2.f;
  }

  //  cerr << "set_avg_leaf_dist: end" << endl;
}

bool UPGMA_Tree::it_is_terminal_node()
{
  return( ( l_node == (UPGMA_Tree*)NULL ) &&
	  ( r_node == (UPGMA_Tree*)NULL ) );
}

void UPGMA_Tree::set_weight()
{
  //  cerr << "set_weight: top" << endl;

  if( it_is_terminal_node() ) {
    //    cerr << "set_weight: about to return 1 for terminal node" << endl;
    weight = 1;
  }
  else {
    //    cerr << "set_weight: about to return sum for non-terminal node" << endl;
    weight = l_node->get_weight() + r_node->get_weight();
  }
}


vector<UPGMA_Tree*> UPGMA_Tree::get_leaves()
{
  vector<UPGMA_Tree*> res;

  if( it_is_terminal_node() ) {
    res.push_back( this );
    return res;
  }
  else if( l_node->it_is_terminal_node() && r_node->it_is_terminal_node() ) { // both are leaves
    res.push_back( l_node );
    res.push_back( r_node );

    return res;
  }
  else if( l_node->it_is_terminal_node() && !r_node->it_is_terminal_node() ) { // only left is leaf
    res.push_back( l_node );

    vector<UPGMA_Tree*> tmp_r;
    tmp_r = r_node->get_leaves();

    for( int i=0; i<(int)tmp_r.size(); i++ ) {
      res.push_back( tmp_r[i] );
    }

    return res;
  }
  else if( !l_node->it_is_terminal_node() && r_node->it_is_terminal_node() ) { // only right is leaf
    res.push_back( r_node );

    vector<UPGMA_Tree*> tmp_l;
    tmp_l = l_node->get_leaves();

    for( int i=0; i<(int)tmp_l.size(); i++ ) {
      res.push_back( tmp_l[i] );
    }

    return res;
  }
  else { // neither branch is a leaf
    vector<UPGMA_Tree*> tmp_l, tmp_r;
    tmp_l = l_node->get_leaves();
    tmp_r = r_node->get_leaves();

    for( int i=0; i<(int)tmp_l.size(); i++ ) {
      res.push_back( tmp_l[i] );
    }

    for( int i=0; i<(int)tmp_r.size(); i++ ) {
      res.push_back( tmp_r[i] );
    }

    return res;
  }

}


// DEBUGGING
void UPGMA_Tree::print_node()
{
  cerr << "---------------" << endl;
  cerr << "node index: " << index << endl;
  cerr << "avg_leaf_dist: " << avg_leaf_dist << endl;

  if( !it_is_terminal_node() ) {
    cerr << "left node index:  " << l_node->get_index() << endl;
    cerr << "left distance:  " << l_dist << endl;
    cerr << "right node index: " << r_node->get_index() << endl;
    cerr << "right distance: " << r_dist << endl;
    cerr << "---------------" << endl;
  }
  else {
    cerr << "terminal node!" << endl;
    cerr << "---------------" << endl;

  }
}

void UPGMA_Tree::print_tree()
{
  print_node();

  if( !it_is_terminal_node() ) {

    if( !l_node->it_is_terminal_node() ) {
      l_node->print_tree();
    }
    else {
      l_node->print_node();
    }

    if( !r_node->it_is_terminal_node() ) {
      r_node->print_tree();
    }
    else {
      r_node->print_node();
    }

  }

}
