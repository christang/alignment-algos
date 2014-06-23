#include "UPGMA_Clusterer.h"

using namespace std;

/*********************************************
* implementation of UPGMA_Clusterer class *
*********************************************/

// constructor 1
UPGMA_Clusterer::UPGMA_Clusterer( float** d, int size )
{

  // load distances
  distance_orig = d;
  num_points_orig = size;

  distance_curr = new float* [num_points_orig];

  for( int i=0; i<num_points_orig; i++ ) {
    distance_curr[i] = new float[i+1];
  }

  for( int i=0; i<num_points_orig; i++ ) {
    for( int j=0; j<=i; j++ ) {
      distance_curr[i][j] = distance_orig[i][j];
    }
  }

  num_points_curr = num_points_orig;

  // setup weight array
  w = new int [ num_points_orig ];
  for( int i=0; i<num_points_orig; i++ ) {
    w[i] = 1;
  }

  nodes.reserve( 2 * num_points_orig ); // enough space for all real and virtual nodes

  //  open_nodes.reserve( size );
  //  closed_nodes.reserve( size );

  // make real nodes
  for( int i=0; i<size; i++ ) {
    nodes.push_back( UPGMA_Tree( (UPGMA_Tree*)NULL, (UPGMA_Tree*)NULL, -1.f, -1.f, i ) );
  }

  // setup node_index array
  node_index = new int [ num_points_orig ];
  for( int i=0; i<num_points_orig; i++ ) {
    node_index[i] = i;
  }

  next_virtual_index = -1;
  next_node_index = num_points_orig;

}

// destructor
UPGMA_Clusterer::~UPGMA_Clusterer()
{
  // clean-up
  delete w;
  delete node_index;

  for( int i=0; i<num_points_curr; i++ ) {
  delete distance_curr[i];
  }
  delete distance_curr;

  delete root;
}

// MAIN

void UPGMA_Clusterer::cluster()
{

  cerr << "clustering..." << endl;

  while( num_points_curr > 2 ) {

    cerr << "num_points_curr: " << num_points_curr << "\t";

    int n1, n2;
    float curr_min_dist = find_closest_pair( n1, n2 );

    cerr << "curr_min_dist: " << curr_min_dist << endl;

    combine_nodes( n1, n2 );
  }

  // add one last node to form the root of the tree

  float min_dist = distance_curr[1][0];

  //  list<UPGMA_Tree>::iterator it = nodes.begin();

  UPGMA_Tree* node_0 = &nodes[ node_index[0] ];
  UPGMA_Tree* node_1 = &nodes[ node_index[1] ];

  root = new UPGMA_Tree( node_0,
			 node_1,
			 ( min_dist / 2.f ) - node_0->get_avg_leaf_dist(),
			 ( min_dist / 2.f ) - node_1->get_avg_leaf_dist(),
			 next_node_index );

  cerr << "root weight: " << root->get_weight() << endl;

}


// GENERAL


void UPGMA_Clusterer::find_clusters_under_threshold( float thresh )
{
  find_clusters_under_threshold( root, thresh );
}


void UPGMA_Clusterer::find_clusters_under_threshold( UPGMA_Tree* start, float thresh )
{

  // begin at start (assuming clustering has been done) and proceed downward recursively until
  // a node has an 'avg_leaf_dist' below threshold

  if( start->it_is_terminal_node() ) {
    vector<UPGMA_Tree*> tmp;
    tmp.push_back( start );
    clusters.push_back( tmp );
    return;
  }

  if( start->get_avg_leaf_dist() < thresh ) {

    // get the leaves under this node (members of the cluster) and save them
    vector<UPGMA_Tree*> tmp = start->get_leaves();
    clusters.push_back( tmp );
  }
  else {
    // recurse further down
    find_clusters_under_threshold( start->get_l_node(), thresh );
    find_clusters_under_threshold( start->get_r_node(), thresh );
  }

}


void UPGMA_Clusterer::print_clusters_under_threshold()
{
  for( int i=0; i<(int)clusters.size(); i++ ) {
    cerr << "cluster #" << i <<  endl;

    for( int j=0; j<(int)clusters[i].size(); j++ ) {
      cerr << "node: " << clusters[i][j]->get_index() << endl;
    }
    cerr << endl;

    vector<UPGMA_Tree*> clstr = clusters[i];

    for( int j=0; j<(int)clstr.size(); j++ ) {
      for( int k=0; k<(int)clstr.size(); k++ ) {
	
	if( k < j ) { cerr << "\t"; }
	else {

	  int min_idx, max_idx;

	  if( clstr[j]->get_index() < clstr[k]->get_index() ) {
	    min_idx = clstr[j]->get_index();
	    max_idx = clstr[k]->get_index();
	  }
	  else{
	    min_idx = clstr[k]->get_index();
	    max_idx = clstr[j]->get_index();
	  }

	  cerr << distance_orig[max_idx][min_idx] << "\t";
	}
	
      }
      cerr << endl;
      
    }

    cerr << endl << endl;
  }
  
}


float UPGMA_Clusterer::find_closest_pair( int& min1, int& min2 )
{
  // Purpose: return the indices of the two closest points in distance_curr

  float min_dist(999999.f);

  for( int i=1; i<num_points_curr; i++ ) {
    for( int j=0; j<i; j++ ) {
      if( distance_curr[i][j] < min_dist ) {
	min_dist = distance_curr[i][j];
	min1 = i;
	min2 = j;
      }
    }
  }

  return min_dist;
}

void UPGMA_Clusterer::combine_nodes( int n1, int n2 )
{
  // Purpose: Combine the nodes represented by the n1 and n2 entries in distance_curr
  //          into a single column in the 0th position.  Adjust 'open_nodes' to reflect
  //          this change.

  // n1 must be the smaller index
  if( n1 > n2 ) {
    int tmp = n1;
    n1 = n2;
    n2 = tmp;
  }

  // move rows n1 and n2 into positions 0 and 1 and adjust values in weight array
  if( n1 != 0 ) {
    swap_cols( n1, 0 ); // swap columns in distance matrix

    int tmp = w[n1]; // swap positions in weight array
    w[n1] = w[0];
    w[0] = tmp;

    tmp = node_index[n1]; // swap positions in node_index array
    node_index[n1] = node_index[0];
    node_index[0] = tmp;

  }

  if( n2 != 1 ) {
    swap_cols( n2, 1 );

    int tmp = w[n2];
    w[n2] = w[1];
    w[1] = tmp;

    tmp = node_index[n2];
    node_index[n2] = node_index[1];
    node_index[1] = tmp;

  }

  // replace the first two nodes in 'nodes' with a new one representing their average
  float min_dist = distance_curr[1][0];

  UPGMA_Tree* l_branch = &nodes[ node_index[0] ];
  UPGMA_Tree* r_branch = &nodes[ node_index[1] ];

  nodes.push_back( UPGMA_Tree( l_branch,
			       r_branch,
			       ( min_dist / 2.f ) - l_branch->get_avg_leaf_dist(),
			       ( min_dist / 2.f ) - r_branch->get_avg_leaf_dist(),
			       next_node_index
			       ) );
  

  // update node_index[] by, deleteing the first two elements, adding the
  // new (virtual) node to the front and moving every other element up by
  // one position (thus shortening the total length by one
  node_index[0] = next_node_index;
  for( int i=2; i<num_points_orig; i++ ) {
    node_index[i-1] = node_index[i];
  }

  next_node_index++; // increment the index for the next (virtual) node

  float dist_curr_new[ num_points_curr - 1 ][ num_points_curr - 1];

  int* w_new = new int[ num_points_curr - 1 ];

  // fill in dist_curr_new
  for( int i=2; i<num_points_curr; i++ ) {
    dist_curr_new[i-1][0] = ( ( w[0] * get_dist_curr(i,0) ) + 
			      ( w[1] * get_dist_curr(i,1) ) ) / ( w[0] + w[1] );
  }

  for( int i=2; i<num_points_curr; i++ ) {
    for( int j=2; j<num_points_curr; j++ ) {
      dist_curr_new[i-1][j-1] = get_dist_curr(i,j);
    }
  }

  for( int i=0; i<(num_points_curr-1); i++ ) {
    dist_curr_new[i][i] = 0;
  }

  // update weight array
  w_new[0] = w[0] + w[1];

  for( int i=2; i<num_points_curr; i++ ) {
    w_new[i-1] = w[i];
  }

  // delete old distance_curr and weight and copy new variables
  delete_distance_matrix( distance_curr, num_points_curr );
  delete w;
  num_points_curr--; // reduce the size of the matrix by one

  w = new int [ num_points_curr ];

  distance_curr = new float* [ num_points_curr ];
  for( int i=0; i<num_points_curr; i++ ) {
    distance_curr[i] = new float [ i+1 ];
  }

  for( int i=0; i<num_points_curr; i++ ) {
    for( int j=0; j<=i; j++ ) {
      distance_curr[i][j] = dist_curr_new[i][j];
    }
  }

  for( int i=0; i<num_points_curr; i++ ) {
    w[i] = w_new[i];
  }

}


void UPGMA_Clusterer::swap_cols( int n1, int n2 )
{
  // Purpose: make new version of distance_curr with n1 and n2 swapped

  float dist_curr_tmp[ num_points_curr ][ num_points_curr ];

  // initialize it to zero
  for( int i=0; i<num_points_curr; i++ ) {
    for( int j=0; j<num_points_curr; j++ ) {
      dist_curr_tmp[i][j] = 0;
    }
  }

  // fill in the rows with cols n1 and n2 at the front
  for( int i=1; i<num_points_curr; i++ ) {
    for( int j=0; j<i; j++ ) {
      if( i == n1 ) {
	if( j == n2 ) {
	  dist_curr_tmp[i][j] = get_dist_curr(j,i);
	}
	else {
	  dist_curr_tmp[i][j] = get_dist_curr(j,n2);
	}
      }
      else if( i == n2 ) {
	if(j == n1 ) {
	  dist_curr_tmp[i][j] = get_dist_curr(j,i);
	}
	else {
	  dist_curr_tmp[i][j] = get_dist_curr(j,n1);
	}
      }
      else {
	if( j == n1 ) {
	  dist_curr_tmp[i][j] = get_dist_curr(n2,i);
	}
	else if( j == n2 ) {
	  dist_curr_tmp[i][j] = get_dist_curr(n1,i);
	}
	else {
	  dist_curr_tmp[i][j] = get_dist_curr(j,i);
	}
      }
    }
  }

  // copy dist_curr_tmp back to distance_curr
  for( int i=0; i<num_points_curr; i++ ) {
    for( int j=0; j<=i; j++ ) {
      distance_curr[i][j] = dist_curr_tmp[i][j];
    }
  }

}


void UPGMA_Clusterer::update_distance_matrix( int a1, int a2 )
{


}


void UPGMA_Clusterer::delete_distance_matrix( float** d, int size )
{
  for( int i=0; i<size; i++ ) {
    delete d[i];
  }

  delete d;
}


// ACCESS

/*
UPGMA_Tree* UPGMA_Clusterer::get_node( int n )
{

  for( int i=0; i<(int)nodes.size(); i++ ) {
    if( nodes[i].get_index() == n ) {
      return &nodes[i];
    }
  }

  return (UPGMA_Tree*)NULL;
}
*/

// DEBUGGING
void UPGMA_Clusterer::print_distance_matrix( float** m, int size )
{
  /*
  for( int i=0; i<size; i++ ) {
    for( int j=0; j<i; j++ ) {
      cerr << "dist: d(" << i << "," << j << ") = " << m[i][j] << endl;
    }
  }
  */

  cerr << "matrix: " << endl;

  for( int j=0; j<size; j++ ) {
    for( int i=0; i<size; i++ ) {
      
      if( i < j ) {
	cerr << "\t";
      }
      else {
	cerr << m[i][j] << "\t";
      }
      
    }
    cerr << endl;
  }

}
