#include "stdio.h"
#include "clusterset.h"
#include <math.h>

using namespace std;

/***********************************************************
 * implementation of Cluster_Set class                       *
 **********************************************************/

// constructor 1
Cluster_Set::Cluster_Set( int size ) 

  : num_points( size ) // set 'num_points'
  
{

  // make space for 'distance' and 'distance_sq' matrices
  distance    = new float* [ num_points ];
  distance_sq = new float* [ num_points ];

  for ( int i=0; i<num_points; i++ ) {
    distance[i]    = new float [ i+1 ]; // avoid redundant distances
    distance_sq[i] = new float [ i+1 ];
  }

  for( int i=0; i<num_points; i++ ) {
    for( int j=0; j<=i; j++ ) {
      distance[i][j]    = -1.f; // initialize with -1
      distance_sq[i][j] = -1.f;
    }
  }

}


// destructor
Cluster_Set::~Cluster_Set() {
 
  // delete 2D distance matrix
  for( int i=0; i<num_points; i++ ) {
    delete distance[i];
    delete distance_sq[i];
  }

  delete distance;
  delete distance_sq;

}



// MAIN




// DEBUGGING

void Cluster_Set::print() {

  for( int i=0; i<num_points; i++ ) {
    for( int j=0; j<=i; j++ ) {
      cerr << "distance[" << i << "][" << j << "]: " << distance[i][j] << endl;
    }
  }


}
