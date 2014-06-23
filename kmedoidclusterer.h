#ifndef _KMEDCLUST
#define _KMEDCLUST

#include <string>
#include <vector>
#include <list>

#include "clusterset.h"

using namespace std;

struct Cluster {
  vector<int> members;
  int centroid;
  float variance;

  Cluster();
  Cluster( int );
  
  inline bool operator< (const Cluster& c) 
  { return members.size() < c.members.size(); }
};

class K_Medoid_Clusterer {

 private:
  Cluster_Set* points; // original source of points
  int num_points; // number of points to cluster
  int k_max; // maximum number of allowed clusters
  float kT; // temperature
  vector<Cluster*> curr;
  vector<Cluster*> next; // proposed next state in simulated annealing algorithm
  vector<Cluster*> temp; // temporary state
  vector<Cluster*> best; // lowest energy state
  vector<Cluster*> final; // final state, lowest energy for the smallest number 
                          // of clusters passing clusters passing 'below_max_var'
  list<int> prev_centroids;
  Cluster* INVALID;

 protected:

 public:

  // data


  // FUNCTIONS
  K_Medoid_Clusterer( Cluster_Set*, int );
  ~K_Medoid_Clusterer();

  // overall clustering
  vector< vector<int> > find_good_clustering( int ); // do in # of clusterings and return the best
  vector< vector<int> > simulated_annealing( float );
  float cluster( vector<Cluster*>& ); // perform k-medoids clustering
  void copy_state( const vector<Cluster*>&, vector<Cluster*>& );
  void clear_state( vector<Cluster*>& );

  // helper functions
  vector<vector<int> > output_state( const vector<Cluster*>& );
  void input_state( const vector<vector<int> >&, vector<Cluster*>& );
  void choose_initial_clusters( vector<Cluster*>& ); // arbitrarily cluster points at first
  void randomly_choose_initial_clusters( vector<Cluster*>& ); // randomly cluster to start
  Cluster* get_random_cluster( const vector<Cluster*>& );
  inline float random_p() { return (float)( rand() % 100 ) / 100.f; }
  list<int> get_sorted_centroids( const vector<Cluster*>& );
  void clear_cluster_members( vector<Cluster*>& );
  void clear_cluster_centroids( vector<Cluster*>& );

  // cluster_functions
  void assign_all_points( vector<Cluster*>& ); 
  void put_with_nearest_centroid( int, vector<Cluster*>& ); // add a point to the right cluster
  void update_centroids( vector<Cluster*>& );
  void update_cluster_centroid( Cluster* );
  float get_cluster_variance( Cluster* );
  void remove_cluster( vector<Cluster*>&, Cluster* );
  void copy_cluster( Cluster*, Cluster* );
  void clear_cluster( Cluster* );
  bool below_max_var( const vector<Cluster*>&, float );

  // simulated annealing functions
  Cluster* merge_clusters( Cluster*, Cluster* );
  vector<Cluster*> split_cluster( Cluster* );
  Cluster* choose_cluster_to_split( const vector<Cluster*>& );
  vector<Cluster*> choose_clusters_to_merge( const vector<Cluster*>& );
  vector<Cluster*> get_nearest_clusters( const vector<Cluster*>& );
  Cluster* get_broadest_cluster( const vector<Cluster*>& );
  void merge_two_clusters( vector<Cluster*>& );
  void split_one_cluster ( vector<Cluster*>& );

  float get_total_variance( const vector<Cluster*>& );
  float get_inter_centroid_distance_sq( Cluster*, Cluster* );

  // debugging
  void print_cluster( Cluster* );
  void print_state( const vector<Cluster*>& );

};



#endif  //_KMEDCLUST
