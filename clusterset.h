#ifndef _CLUSTERSET
#define _CLUSTERSET

#include <iostream>
#include <string>
#include <vector>

#include "math.h"

//#include "alignment.h"

using namespace std;

class Cluster_Set {

 private:
  float** distance;
  float** distance_sq;
  int num_points;

 protected:

 public:

  // data


  // functions
  Cluster_Set( int size );
  ~Cluster_Set();

  // main

  // access
  inline float dist( int i, int j )    { return ( i>=j ? distance[i][j]    : distance[j][i] ); }
  inline float dist_sq( int i, int j ) { return ( i>=j ? distance_sq[i][j] : distance_sq[j][i] ); }
  inline void set_dist( int i, int j, float d )    { distance[i][j]    = d; }
  inline void set_dist_sq( int i, int j, float d ) { distance_sq[i][j] = pow( d, 2 ); }
  inline int get_num_points() { return num_points; }

  // debug
  void print();

};



#endif  //_CLUSTERSET
