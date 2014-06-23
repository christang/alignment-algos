#include "ali_str_info.h"

using namespace std;

/*********************************************
* implementation of Ali_Str_Info class *
*********************************************/

// constructor 1
Ali_Str_Info::Ali_Str_Info()
{}

Ali_Str_Info::~Ali_Str_Info()
{

  // delete similarity matrix, sims
  for( int i=0; i<query_len; i++ ) {
    delete sims[i];
  }
  delete sims;

  // delete the inter-residue distance matrix
  for( int i=0; i<templ_len; i++ ) {
    delete cb_dists[i];
  }
  delete cb_dists;

  // delete the contact matrix
  for( int i=0; i<templ_len; i++ ) {
    delete templ_contacts[i];
  }
  delete templ_contacts;

  // delete the query_predicted_loops array
  delete query_predicted_loops;

  // delete the TSR arrays
  delete TSR_to_N;
  delete TSR_to_C;

}

// OPERATORS

// SETUP
void Ali_Str_Info::load_seq_lengths( int t_len, int q_len )
{
  // load the template and query sequence lengths
  // necessary to do this before loading any other arrays as these ...
  // lengths determine the boundaries of the array

  templ_len = t_len;
  query_len = q_len;
}


void Ali_Str_Info::load_seq_strings( string templ, string query )
{
  // load the sequence strings
  templ_seq = templ;
  query_seq = query;
}


void Ali_Str_Info::load_sims( float** data )
{
  // setup the similarity matrix, sims**
  sims = new float* [ query_len ];

  for( int i=0; i<query_len; i++ ) {
    sims[i] = new float[ templ_len ];
  }

  for( int i=0; i<query_len; i++ ) {
    for( int j=0; j<templ_len; j++ ) {
      sims[i][j] = data[i][j];
    }
  }

}


void Ali_Str_Info::load_cb_dists( float** data )
{
  // setup the CB_distance matrix, cb_dists**
  cb_dists = new float* [ templ_len ];

  for( int i=0; i<templ_len; i++ ) {
    cb_dists[i] = new float [ templ_len ];
  }

  for( int i=0; i<templ_len; i++ ) {
    for( int j=0; j<templ_len; j++ ) {
      cb_dists[i][j] = data[i][j];
    }
  }

}


void Ali_Str_Info::load_contacts( bool** data )
{
  // setup the contact matrix, templ_contacts**
  templ_contacts = new bool* [ templ_len ];

  for( int i=0; i<templ_len; i++ ) {
    templ_contacts[i] = new bool [ templ_len ];
  }

  for( int i=0; i<(int)templ_len; i++ ) {
    for( int j=0; j<(int)templ_len; j++ ) {
      templ_contacts[i][j] = data[i][j];
    }
  }

}


void Ali_Str_Info::load_query_predicted_loops( bool* data )
{
  // setup the query predicted loops array
  query_predicted_loops = new bool [ query_len ];

  for( int i=0; i<query_len; i++ ) {
    query_predicted_loops[i] = data[i];
  }

}


void Ali_Str_Info::load_SSE_Data( vector<SSE_Data> data )
{
  templ_SSEs = data;
  num_templ_SSEs = templ_SSEs.size();
}

void Ali_Str_Info::load_TSRs( int* data_N, int* data_C )
{
  // setup the N- and C-terminal TSR arrays

  TSR_to_N = new int [ templ_len ];
  TSR_to_C = new int [ templ_len ];

  for( int i=0; i<templ_len; i++ ) {
    TSR_to_N[i] = data_N[i];
    TSR_to_C[i] = data_C[i];
  }

}


// ACCESS

// DEBUGGING
