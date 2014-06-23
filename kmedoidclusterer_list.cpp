#include "kmedoidclusterer.h"

#include <math.h>
#include <time.h>

using namespace std;

Cluster::Cluster() {
  centroid = -1;
  members.clear();
  variance = 0;
}

Cluster::Cluster( int cent ) {
  centroid = cent;
  members.clear();
  variance = 0;
}

/*********************************************
* implementation of K_Medoid_Clusterer class *
*********************************************/

// constructor 1
K_Medoid_Clusterer::K_Medoid_Clusterer( Cluster_Set* cs, int num )
  : 
  points( cs ),
  k_max( num )
{

  // setup 'curr' and 'next' states
  for( int i=0; i<k_max; i++ ) {
    curr.push_back( new Cluster );
    next.push_back( new Cluster );
  }

  num_points = points->get_num_points();   // initiate 'num_points'

  INVALID = new Cluster( -1 ); // setup the INVALID cluster (a flag)
  INVALID->centroid = -1;

}


// destructor
K_Medoid_Clusterer::~K_Medoid_Clusterer() {

  // delete states
  for( unsigned int i=0; i<curr.size(); i++ ) {
    delete curr[i];
  }

  for( unsigned int i=0; i<next.size(); i++ ) {
    delete next[i];
  }

  delete INVALID;
}


vector< vector<int> > K_Medoid_Clusterer::find_good_clustering( int n ) {

  // perform n clusterings, beginning with a different random assignment of clusters
  // each time.  return the best.

  srand( time(NULL) ); // seed the random number generator

  //  vector< vector<int> > best_clustering;

  clear_cluster_members( curr );
  clear_cluster_centroids( curr );

  randomly_choose_initial_clusters( curr );

  float min_variance = cluster( curr ); // set the min_distance variable with the first result
  //  best_clustering = save_clusters();
  save_state( curr );

  float curr_variance;

  for( int i=0; i<n; i++ ) {

    clear_cluster_members( curr );
    clear_cluster_centroids( curr );

    randomly_choose_initial_clusters( curr );

    curr_variance = cluster( curr );

    if( curr_variance < min_variance ) {
      min_variance = curr_variance;
      save_state( curr );
    }

  }

  //  return best_clustering;
  return saved_state;
}


vector< vector<int> > K_Medoid_Clusterer::simulated_annealing( float max_var ) {

  cerr << "max_var: " << max_var << endl;

  randomly_choose_initial_clusters( curr );
  update_centroids( curr );

  float e = cluster( curr ); // energy of 'curr' state
  float e_next; // energy of 'next' state
  float e_best; // lowest energy found so far

  long int max_iterations = curr.size() * curr.size();

  kT = e; // starting point for kT is the initial energy/variance
  cerr << "initial kT: " << kT << endl;

  copy_state( curr, best ); // initially, this is the lowest energy state
  e_best = e;

  while( kT > max_var ) {

    copy_state( best, curr ); // start each run with the best state found so far

    for( int i=0; i<max_iterations; i++ ) {

      copy_state( curr, next ); // setup neighboring state
      merge_two_clusters( next ); // get a neighboring state
      split_one_cluster ( next );
      e_next = cluster( next ); // km cluster and get energy of 'next'
    
      if( e_next < e_best ) { // energy is lowest so far
	copy_state( next, best );
	e_best = e_next;
      }
      
      if( exp( -( e_next - e ) / kT ) > random_p() ) { // then proceed to next state
	copy_state( next, curr ); // curr := next
	e = e_next;
      }
    }

    kT *= 0.9; // cool system
    cerr << "kT: " << kT << ", e_best: " << e_best << endl;
  }



  return saved_state;
}


float K_Medoid_Clusterer::cluster( vector<Cluster*>& vc ) {

  list<int> prev_centroids;
  list<int> curr_centroids = get_sorted_centroids( vc );

  while( prev_centroids != curr_centroids ) { // check if clusters changed

    prev_centroids = curr_centroids;

    update_centroids( vc );

    assign_all_points( vc );
  }

  return( get_total_variance( vc ) );
}


void K_Medoid_Clusterer::copy_state( const vector<Cluster*>& orig, vector<Cluster*>& copy ) {

  // make sure 'orig' and 'copy' have the same number of clusters (i.e., same size)

  while( copy.size() > orig.size() ) {
    delete copy.back(); // delete the last cluster pointed to
    copy.pop_back(); // drop the last element of the vector
  }

  while( copy.size() < orig.size() ) {
    copy.push_back( new Cluster ); // create a new cluster at the end of 'copy'
  }

  clear_state( copy );

  // make a fresh copy of the data pointed to in 'orig'
  for( unsigned int i=0; i<orig.size(); i++ ) {
    copy_cluster( orig[i], copy[i] );
  }
}

void K_Medoid_Clusterer::clear_state( vector<Cluster*>& vc ) {

  // reset all clusters
  for( unsigned int i=0; i<vc.size(); i++ ) {
    clear_cluster( vc[i] );
  }
}


void K_Medoid_Clusterer::save_state( const vector<Cluster*>& vc ) {

  // store state 'vc' in private variable 'saved_state'

  saved_state.clear();

  vector<int> tmp;

  for( unsigned int i=0; i<vc.size(); i++ ) {

    tmp.clear();
    tmp.push_back( vc[i]->centroid );

    // store all cluster members with the centroid at the front
    for( list<int>::iterator it = vc[i]->members.begin();
	 it != vc[i]->members.end(); it++ ) {

      if( *it != tmp.front() ) {
	tmp.push_back( *it );
      }

    }

    saved_state.push_back( tmp );

  }

}


void K_Medoid_Clusterer::load_state( vector<Cluster*>& vc ) {

  k_max = saved_state.size();

  float variance;

  // delete current 'vc'
  clear_state( vc );

  // create new 'vc' from 'saved_state'
  for( unsigned int i=0; i<saved_state.size(); i++ ) {

    vc.push_back( new Cluster( saved_state[i].front() ) );

    variance = 0;

    for( unsigned int j=0; j<saved_state[i].size(); j++ ) {
      vc.back()->members.push_back( saved_state[i][j] );
      variance += points->dist_sq( saved_state[i].front(), saved_state[i][j] );
    }

    vc.back()->variance = variance / vc.back()->members.size() ;
  }
}


void K_Medoid_Clusterer::choose_initial_clusters( vector<Cluster*>& vc) {

  clear_state( vc );

  int step = num_points / k_max;

  cerr << "step: " << step << endl;

  for( int i=0; i<(k_max-1); i++ ) {
    for( int j=0; j<step; j++ ) {
      vc[i]->members.push_back( i*step + j );
    }
  }

  for( int j=(k_max-1)*step; j<num_points; j++ ) {
    vc[k_max-1]->members.push_back( j );
  }


}


void K_Medoid_Clusterer::randomly_choose_initial_clusters( vector<Cluster*>& vc ) {

  clear_cluster_members( vc );
  clear_cluster_centroids( vc );

  // make sure each cluster has at least one member
  for( unsigned int i=0; i<vc.size(); i++ ) {
    vc[i]->members.push_back( i );
  }

  for( int i=vc.size(); i<num_points; i++ ) {
    get_random_cluster( vc )->members.push_back( i );
  }

  // set variance
  for( unsigned int i=0; i<vc.size(); i++ ) {
    vc[i]->variance = get_cluster_variance( vc[i] );
  }
}

void K_Medoid_Clusterer::assign_all_points( vector<Cluster*>& vc ) {

  // clear all previous members from cluster lists
  clear_cluster_members( vc );

  // put each point in the cluster with the centroid it is nearest
  for( int i=0; i<num_points; i++ ) {
    put_with_nearest_centroid( i, vc );
  }

}


void K_Medoid_Clusterer::put_with_nearest_centroid( int p, vector<Cluster*>& vc ) {

  float min_centroid_dist_sq = points->dist_sq( p, vc[0]->centroid );
  int min_cluster = 0;

  float dist_sq;

  for( unsigned int i=1; i<vc.size(); i++ ) {

    dist_sq = points->dist_sq( p, vc[i]->centroid );

    if( dist_sq < min_centroid_dist_sq ) {
      min_centroid_dist_sq = dist_sq;
      min_cluster = i;
    }
  }

  // update cluster variance
  vc[ min_cluster ]->variance += ( min_centroid_dist_sq - vc[ min_cluster ]->variance ) / 
    ( vc[ min_cluster ]->members.size() + 1 );

  vc[ min_cluster ]->members.push_back( p ); // put new point in cluster
}


void K_Medoid_Clusterer::update_centroids( vector<Cluster*>& vc ) {

  // look at the members of each cluster and choose the point that is the best centroid
  // i.e., has the smallest total squared distance to each other point in the cluster

  //  float min_cluster_distance;
  //  float min_cluster_variance, curr_cluster_variance;

  for( unsigned int i=0; i<vc.size(); i++ ) {
    update_cluster_centroid( vc[i] );
  }

}


void K_Medoid_Clusterer::update_cluster_centroid( Cluster* clstr ) {

    // initially, min_cluster_variance is the sum of squared distances between the 
    // first member of the cluster and all other members

    float min_cluster_variance = 0;

    list<int>::iterator begin = clstr->members.begin();
    list<int>::iterator end   = clstr->members.end();
    list<int>::iterator min   = begin;
    list<int>::iterator it;

    for( it=begin; it != end; it++ ) {
      min_cluster_variance += points->dist_sq( *it, *min );
    }

    // now that min_cluster_variance has a value, do the same for all other members of
    // the cluster and keep track of the one with the smallest variance

    list<int>::iterator curr;
    curr++; // start testing from the 1st position
    float curr_cluster_variance;

    for( curr=begin; curr != end; curr++ ) {

      curr_cluster_variance = 0; // reset

      for( it = begin; it != end; it++ ) {
      	curr_cluster_variance += points->dist_sq( *it, *curr );
      }

      if( curr_cluster_variance < min_cluster_variance ) {
	min_cluster_variance = curr_cluster_variance;
	min = curr;
      }

    }

    clstr->centroid = *min; // update the centroid with the new minimum
    clstr->variance = min_cluster_variance / clstr->members.size();
}


float K_Medoid_Clusterer::get_cluster_variance( Cluster* c ) {

  if( c->members.size() == 0 ) { return -1.f; }

  float variance = 0;

  for( list<int>::iterator it = c->members.begin();
       it != c->members.end(); it++ ) {
    variance += points->dist_sq( c->centroid, *it );
  }

  return variance / c->members.size();
}


void K_Medoid_Clusterer::remove_cluster( vector<Cluster*>& vc, Cluster* c ) {

  vector<Cluster*>::iterator it;

  for( it=vc.begin(); it!=vc.end(); it++ ) {
    if( *it == c ) {
      vc.erase( it );
      return;
    }
  }
}

void K_Medoid_Clusterer::copy_cluster( Cluster* orig, Cluster* copy ) {

  //  clear_cluster( copy );

  // copy cluster members
  list<int>::iterator it;
  for( it=orig->members.begin(); it!=orig->members.end(); it++ ) {
    copy->members.push_back( *it );
  }

  // copy centroid and variance
  copy->centroid = orig->centroid;
  copy->variance = orig->variance;

}


void K_Medoid_Clusterer::clear_cluster( Cluster* c ) {
  c->members.clear();
  c->centroid = -1;
  c->variance = 0;
}

void K_Medoid_Clusterer::merge_two_clusters( vector<Cluster*>& vc ) {

  // merge clusters
  vector<Cluster*> to_merge = choose_clusters_to_merge( vc );
  if( to_merge[0] == (Cluster*)NULL ) { 
    to_merge = get_nearest_clusters( vc );
  }
  Cluster* merged = merge_clusters( to_merge[0], to_merge[1] );
    
  // remove original clusters that merged from state and add the new one
  remove_cluster( vc, to_merge[0] );
  remove_cluster( vc, to_merge[1] );
  vc.push_back( merged );
}


void K_Medoid_Clusterer::split_one_cluster ( vector<Cluster*>& vc ) {

  //  cerr << "top of split_one_cluster" << endl;

  // split cluster
  Cluster* to_split = choose_cluster_to_split( vc );
  //  cerr << "chose cluster" << endl;
  if( to_split->centroid == -1 ) { // warning flag for bad cluster choice
    //    cerr << "to_split is INVALID" << endl;
    to_split = get_broadest_cluster( vc );
    //    cerr << "after get_broadest, to_split is: " << endl;
    //    print_cluster( to_split );
    //    string pause; cin >> pause; 
  }

  //  cerr << "past get_broadest_cluster (if necessary)" << endl;
  //  cerr << "to_split is now: " << endl;
  //  print_cluster( to_split );

  vector<Cluster*> split = split_cluster( to_split );
  
  //  cerr << "split cluster" << endl;

  // remove original cluster that was split and add the two new ones
  remove_cluster( vc, to_split );
  vc.push_back( split[0] );
  vc.push_back( split[1] );

  //  cerr << "bottom of split_one_cluster" << endl;
}


Cluster* K_Medoid_Clusterer::merge_clusters( Cluster* clstr1, Cluster* clstr2 ) {

  Cluster* res = new Cluster;

  // copy 'clstr1' into 'res'
  copy_cluster( clstr1, res );

  // add 'clstr2' members to 'res'
  for( list<int>::iterator it=clstr2->members.begin(); it!=clstr2->members.end(); it++ ) {
    res->members.push_back( *it );
  }

  // find the proper centroid for res
  update_cluster_centroid( res );

  return res;
}


vector<Cluster*> K_Medoid_Clusterer::split_cluster( Cluster* clstr ) {

  // split 'clstr' into two new clusters.  then delete the original 'clstr'
  // and put the two new clusters at the end of the 'clusters' vector

  // split the members of 'clstr' by finding the two points which are furthest apart

  //  cerr << "size of cluster to split: " << clstr->members.size() << endl;

  if( clstr->members.size() <= 1 ) {
    cerr << "Cannot split a cluster with only one member." << endl;
    cerr << "Cluster centroid is: " << clstr->centroid << endl;    
    cerr << "Exiting." << endl;
    exit(-1);
  }

  list<int>::iterator begin  = clstr->members.begin();
  list<int>::iterator second = clstr->members.begin(); second++;
  list<int>::iterator end    = clstr->members.end();
  list<int>::iterator penult = clstr->members.end(); penult--;

  list<int>::iterator far1 = begin;
  list<int>::iterator far2 = second;
  float max_distance = points->dist( *far1, *far2 );
  float cur_distance;
  
  for( list<int>::iterator it = begin; it != penult; it++ ) {

    for( list<int>::iterator jt = it; jt != end; jt++ ) { // note: slightly redundant

      cur_distance = points->dist( *it, *jt );

      if( cur_distance > max_distance ) {
	max_distance = cur_distance;
	far1 = it; far2 = jt;
      }

    }
  }

  //  cerr << "found far1 and far2: " << *far1 << ", " << *far2 << endl;

  // create a vector of two new clusters with the centroids the two furthest apart points found above

  vector<Cluster*> res;
  res.push_back( new Cluster( *far1 ) );
  res.push_back( new Cluster( *far2 ) );

  //  cerr << "created res clusters" << endl;

  // move each point from the original 'clstr' to either 'res[0]' or 'res[1]' depending
  // on which has the nearest centroid

  for( list<int>::iterator it = begin; it!= end; it++ ) {
    put_with_nearest_centroid( *it, res );
  }

  //  cerr << "about to return res" << endl;

  return res;
}


float K_Medoid_Clusterer::get_total_variance( const vector<Cluster*>& vc) {

  float variance = 0;

  for( unsigned int i=0; i<vc.size(); i++ ) {
    variance += vc[i]->variance * (float)vc[i]->members.size();
  }

  return variance / (float)num_points;
}


float K_Medoid_Clusterer::get_inter_centroid_distance_sq( Cluster* c1, Cluster* c2 ) {
  return( points->dist_sq( c1->centroid, c2->centroid ) );
}


Cluster* K_Medoid_Clusterer::choose_cluster_to_split( const vector<Cluster*>& vc ) {

  float total_variance = get_total_variance( vc );

  Cluster* cand = get_random_cluster( vc );

  int i=0;
  int max_attempts = 2*vc.size();

  // delta is???

  // keep getting new random clusters until the condition is met
  while( i<max_attempts ) {

    if( ( ( cand->variance / total_variance ) > random_p() ) &&
	( cand->members.size() > 1 ) ) {
      break;
    }

    cand = get_random_cluster( vc );
    i++;
  }

  if( i >= max_attempts ) { 
    //    cerr << "Could not find valid cluster to split." << endl;
    cand = INVALID;
  }
  
  //  cerr << "about to return candidate, which has centroid: " << cand->centroid << endl;

  return cand;
}

vector<Cluster*> K_Medoid_Clusterer::choose_clusters_to_merge( const vector<Cluster*>& vc ) {

  // get two random clusters (that must be different)
  Cluster* cand1 = get_random_cluster( vc );
  Cluster* cand2 = cand1;
  while( cand2 == cand1 ) { cand2 = get_random_cluster( vc ); }

  // use square of distance between candidate's centroids as delta

  int i=0;
  int max_attempts = 10*vc.size();

  // keep getting new pairs of random clusters until the condition is met
  while( ( i<max_attempts ) &&
	 ( exp( get_inter_centroid_distance_sq( cand1, cand2 ) / kT ) < random_p() ) ) {

    // get new candidates
    cand1 = get_random_cluster( vc );
    cand2 = cand1;
    while( cand2 == cand1 ) { cand2 = get_random_cluster( vc ); }
    i++;
  }

  vector<Cluster*> res;

  if( i >= max_attempts ) { // nothing found, nothing chosen
    res.push_back( (Cluster*)NULL );     res.push_back( (Cluster*)NULL );
  }
  else {
    res.push_back( cand1 ); res.push_back( cand2 );

  }

  return res;
}


vector<Cluster*> K_Medoid_Clusterer::get_nearest_clusters( const vector<Cluster*>& vc ) {

  float min_distance = points->dist( vc[0]->centroid, vc[1]->centroid );
  int min_idx1 = 0;
  int min_idx2 = 1;
  float dist;

  for( unsigned int i=0; i<vc.size()-1; i++ ) {
    for( unsigned int j=i+1; j<vc.size(); j++ ) {

      dist = points->dist( vc[i]->centroid, vc[j]->centroid );

      if( dist < min_distance ) {
	min_distance = dist;
	min_idx1 = i;
	min_idx2 = j;
      }
    }
  }

  vector<Cluster*> res;
  res.push_back( vc[min_idx1] );
  res.push_back( vc[min_idx2] );

  return res;
}


Cluster* K_Medoid_Clusterer::get_broadest_cluster( const vector<Cluster*>& vc ) {

  float max_variance = -1.f;
  int max_idx = -1;
  float var;

  for( unsigned int i=1; i<vc.size(); i++ ) {

    var = vc[i]->variance;

    if( ( var > max_variance ) && ( vc[i]->members.size() > 1 ) ) {
      max_variance = var;
      max_idx = i;
    }
  }

  if( max_idx == -1 ) {
    cerr << "no cluster found with more than one member! exiting." << endl;
    exit(-1);
  }

  return vc[max_idx];
}


Cluster* K_Medoid_Clusterer::get_random_cluster( const vector<Cluster*>& vc ) {
  return( vc[ rand() % vc.size() ] );
}


list<int> K_Medoid_Clusterer::get_sorted_centroids( const vector<Cluster*>& vc ) {

  list<int> res;

  for( unsigned int i=0; i<vc.size(); i++ ) {
    res.push_back( vc[i]->centroid );
  }

  res.sort();

  return res;
}


void K_Medoid_Clusterer::clear_cluster_members( vector<Cluster*>& vc ) {

  for( unsigned int i=0; i<vc.size(); i++ ) {
    vc[i]->members.clear();
  }

}


void K_Medoid_Clusterer::clear_cluster_centroids( vector<Cluster*>& vc ) {

  for( unsigned int i=0; i<vc.size(); i++ ) {
    vc[i]->centroid = -1;
    vc[i]->variance = 0;
  }

}


// DEBUGGING

void K_Medoid_Clusterer::print_cluster( Cluster* c ) {

  cerr << "Centroid : " << c->centroid << ", Variance: " << c->variance;
  cerr << ", #: " << c->members.size() << endl;

  cerr << "{ ";  
  for( list<int>::iterator it = c->members.begin(); it != c->members.end(); it++ ) {
      cerr << *it << " ";
  }
  cerr << "}" << endl;
}


void K_Medoid_Clusterer::print_state( const vector<Cluster*>& vc ) {

  for( unsigned int i=0; i<vc.size(); i++ ) {
    print_cluster( vc[i] );
    cerr << endl;
  }
}
