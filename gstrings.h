/**
 *  Package HMAP2.1
 *  File: gstrings.h
 *  Desc: Generate gapped strings for displaying alignments
 *
 *  Created on 11/26/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_GSTRINGS
#define _HMAP2_GSTRINGS

#include <algorithm>
#include <assert.h>
#include <cctype>
#include <string>

#include "alignment.h"

using namespace std;

class SequenceGaps {

  static const char gapChar;

public:

  template <class S1, class S2, class Etype>
  SequenceGaps (const AlignmentSet<S1,S2,Etype>& as,
		valarray<bool> mask = valarray<bool>(0));

  template <class S1, class S2, class Etype>
  SequenceGaps (const AlignmentSet<S1,S2,Etype>& as,
		int qs, int ts,
		valarray<bool> mask = valarray<bool>(0));

  void build (const string& sequence,
	      string& result,
	      char gc = gapChar);

  template <class S1, class S2>
  void build (const string& sequence,
	      const AlignedPairList<S1,S2>& a, 
	      string& result,
	      char gc = gapChar);

private:
  
  template <class S1, class S2, class Etype>
  void buildAnchors (const AlignmentSet<S1,S2,Etype>& as,
		     valarray<bool> mask);

  int  query_len;
  int  template_len;
  vector<int>   anchors;
  int         gap_total;

};

template <class S1, class S2, class Etype>
SequenceGaps::SequenceGaps (const AlignmentSet<S1,S2,Etype>& as,
			    valarray<bool> mask)
  : query_len    ((int)as.getQuerySequence()->size()),
    template_len ((int)as.getTemplateSequence()->size()),
    anchors      (template_len-1,0),
    gap_total    (0)
{ buildAnchors (as,mask);}


template <class S1, class S2, class Etype>
SequenceGaps::SequenceGaps (const AlignmentSet<S1,S2,Etype>& as,
			    int qs, int ts, valarray<bool> mask)
  : query_len    (qs),
    template_len (ts),
    anchors      (template_len-1,0),
    gap_total    (0)
{ buildAnchors (as,mask);}


template <class S1, class S2, class Etype>
void SequenceGaps::buildAnchors (const AlignmentSet<S1,S2,Etype>& as,
				 valarray<bool> mask)
{ 
  typedef AlignedPair<S1,S2>     Pair;
  typedef AlignedPairList<S1,S2> Alignment;
  
  bool do_all=false;
  if ((int)mask.size()==0) do_all = true;
  assert (do_all || mask.size()==as.size());

  for (int i=0; i<(int)as.size(); ++i) {
    if (do_all || mask[i]) {
      typename list<Pair>::const_iterator pair_it  = as[i].begin();
      typename list<Pair>::const_iterator pair_end = as[i].end();
      const Pair* prev = &*pair_it; ++pair_it;
      for (; pair_it!=pair_end; ++pair_it) {
	// Revised to deal with zigzag alignments
	if (pair_it->query_idx() != prev->query_idx()+1) {
	  int x = pair_it->query_idx();
	  int y = x - prev->query_idx();
	  int gap = y-1;
	  anchors[prev->template_idx()]=
	    max(anchors[prev->template_idx()],gap);
	}
	prev = &*pair_it;
      }
    }
  }

  for (int i=0; i<template_len-1; ++i) gap_total+=anchors[i];

}


template <class S1, class S2>
void SequenceGaps::build (const string& seq,
			  const AlignedPairList<S1,S2>& ali,
			  string& value,
			  char gc)
{
  string result = "";

  assert (query_len==(int)seq.size());

  typedef AlignedPair<S1,S2> Pair;

  typename list<Pair>::const_iterator pair_it  = ali.begin();

  for (int j=0; j<template_len-1; ++j) {
    int a_gap = anchors[j]+1;
    const Pair* curr = &*pair_it;
    if (curr->template_idx()==j) {
      ++pair_it;
      int a = curr->template_idx();
      int x = curr->query_idx();

      int b = pair_it->template_idx();
      int y = pair_it->query_idx();
      if (b-a == 1 || y-x == 1) {
	string subseq = seq.substr(x,y-x);
	result.append(subseq);
	a_gap -= y-x;
      } else { // Deal with zigzag alignments
	string subseq = seq.substr(x,y-x);
	int (*pf)(int)=tolower; 
	transform(subseq.begin()+1,subseq.end(),
		  subseq.begin()+1,pf);
	result.append(subseq);
	a_gap -= y-x;
      } 
    }
    
    result.append(a_gap,gc);
  }
  
  int z = template_len+gap_total-(int)result.size();
  if (z>1) result.append (z-1,gc);

  result.append(1,seq[(int)seq.size()-1]);
  value = result;
}

#endif  //_HMAP2_GSTRINGS
