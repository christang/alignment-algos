#include "ali_strand_eval.h"

using namespace std;

/*********************************************
* implementation of Alignment_Strand_Evaluator class *
*********************************************/

// constructor 1
Alignment_Strand_Evaluator::Alignment_Strand_Evaluator()
{}

// destructor
Alignment_Strand_Evaluator::~Alignment_Strand_Evaluator(  )
{}

// GENERAL
void Alignment_Strand_Evaluator::load_SSE_contacts( int size, bool** contacts )
{
  num_sses = size;
  SSE_contacts = contacts;
}

void Alignment_Strand_Evaluator::load_All_Strands( vector<SSE_Data> SSEs )
{
  All_Strands.clear();

  for( int i=0; i<(int)SSEs.size(); i++ )
    {
      if( SSEs[i].ss_type == 330 ) // SSE is a strand
	{ 
	  All_Strands.push_back( SSEs[i].sse_id );
	}

    }

}

void Alignment_Strand_Evaluator::determine_rules()
{

  for( int i=0; i<(int)All_Strands.size(); i++ ) {

    int num_partners(0);

    for( int j=0; j<i; j++ ) {
      if( SSE_contacts[ All_Strands[i] ][ All_Strands[j] ] ) {
	num_partners++;
      }
    }

    for( int k=i+1; k<(int)All_Strands.size(); k++ ) {
      if( SSE_contacts[ All_Strands[k] ][ All_Strands[i] ] ) {
	num_partners++;
      }
    }

    if( num_partners == 1 ) {
      Edge_Strands.push_back( All_Strands[i] );
    }
    else if( num_partners > 1 ) {
      Core_Strands.push_back( All_Strands[i] );
    }

  }

  // setup All_Strands_Paired rule to ensure no strands are unpaired in an alignment

  for( int i=0; i<(int)All_Strands.size(); i++ ) {

    list<int> tmp;

    tmp.push_back( All_Strands[i] );

    for( int j=0; j<i; j++ ) {
      if( SSE_contacts[ All_Strands[i] ][ All_Strands[j] ] ) {
	tmp.push_back( All_Strands[j] );
      }
    }

    for( int k=i+1; k<(int)All_Strands.size(); k++ ) {
      if( SSE_contacts[ All_Strands[k] ][ All_Strands[i] ] ) {
	tmp.push_back( All_Strands[k] );
      }
    }

    All_Strands_Paired.push_back( tmp );
  }

  // setup No_Missing_Cores rule to ensure that no core strand is left out while
  // its two flanking strands are present

  for( int i=0; i<(int)Core_Strands.size(); i++ ) {

    vector<int> All_Partners;

    for( int j=0; j<(int)All_Strands.size(); j++ ) {
      if( Core_Strands[i] > All_Strands[j] ) {
	if( SSE_contacts[ Core_Strands[i] ][ All_Strands[j] ] ) {
	  All_Partners.push_back( All_Strands[j] );
	}
      }
      else if( All_Strands[j] > Core_Strands[i] ) {
	if( SSE_contacts[ All_Strands[j] ][ Core_Strands[i] ] ) {
	  All_Partners.push_back( All_Strands[j] );
	}
      }
    }

    for( int j=1; j<(int)All_Partners.size(); j++ ) {
      for( int k=0; k<j; k++ ) {
	list<int> tmp;
	tmp.push_back( All_Partners[k] );
	tmp.push_back( All_Partners[j] );
	tmp.push_back( Core_Strands[i] );
	No_Missing_Cores.push_back( tmp );
      }
    }

  }

}


bool Alignment_Strand_Evaluator::ali_passes_rules( const list<int>& sse_id_list )
{

  // check All_Strands_Paired

  for( int i=0; i<(int)All_Strands_Paired.size(); i++ ) {

    if( list_contains( sse_id_list, All_Strands_Paired[i].front() ) ) { // this strand is in ali

      bool sse_is_in_list = false;
      list<int>::const_iterator it = All_Strands_Paired[i].begin(); it++;

      while( !sse_is_in_list && it!=All_Strands_Paired[i].end() ) {

	sse_is_in_list = list_contains( sse_id_list, *it );

	it++;
      }

      if( !sse_is_in_list ) { return false; } // failed rule All_Strands_Paired[i]

    }

  }

  // passed every rule in All_Strands_Paired

  // next, check No_Missing_Cores

  for( int i=0; i<(int)No_Missing_Cores.size(); i++ ) {

    list<int>::const_iterator it = No_Missing_Cores[i].begin();
    int strand1, strand2, core_strand;

    strand1 = *it; it++;
    strand2 = *it; it++;
    core_strand = *it;

    if( list_contains( sse_id_list, strand1 ) && list_contains( sse_id_list, strand2 ) ) {

      if( !list_contains( sse_id_list, core_strand ) ) {
	return false;
      }
    }
  }

  // passed every rule in No_Missing_Cores

  // ... more rules?

  // at this point, you've passed every rule, so return true
  return true; 


}


bool Alignment_Strand_Evaluator::list_contains( const list<int>& lst, int n )
{

  for( list<int>::const_iterator it = lst.begin(); it!=lst.end(); it++ ) {

    if( *it == n ) {
      return true;
    }

  }

  return false;
}
