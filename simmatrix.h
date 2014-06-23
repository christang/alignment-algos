/**
 *  Package HMAP2.1
 *  File: simmatrix.h
 *  Desc: Create similarity matrix for dynamic programming
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_SIMMATRIX
#define _HMAP2_SIMMATRIX

#include "matrix.h"

class SimilarityMatrix : public matrix<float>
{

  typedef matrix<float> Base;

public:
  
  template <class S1, class S2, class Etype>
  SimilarityMatrix (const S1& query_seq,
		    const S2& templ_seq,
		    const Evaluator<S1,S2,Etype>& eval);
  
private:

  template <class S1, class S2, class Etype>
  void build (const S1& query_seq,
	      const S2& templ_seq,
	      const Evaluator<S1,S2,Etype>& eval);
  
};


template <class S1, class S2, class Etype>
SimilarityMatrix::SimilarityMatrix (const S1& qs,
				    const S2& ts,
				    const Evaluator<S1,S2,Etype>& eval  )
  : Base      (qs.size(),ts.size())
{
  build (qs,ts,eval);
  eval.post_process (*this);
}

template <class S1, class S2, class Etype>
void SimilarityMatrix::build (const S1& qs,
			      const S2& ts,
			      const Evaluator<S1,S2,Etype>& eval  )
{
  int q_last = rows()-1;
  int t_last = cols()-1;

  for (int i=0; i<=q_last; ++i) {
    operator() (i,0) = 0.f;
    operator() (i,t_last) = 0.f;
  }

  for (int j=0; j<=t_last; ++j) {
    operator() (0,j) = 0.f;
    operator() (q_last,j) = 0.f;
  }

  for (int i=1; i<q_last; ++i) {
    for (int j=1; j<t_last; ++j) {
      operator() (i,j) = eval.similarity(qs,ts,i,j);
    }
  }
}


#endif  // _HMAP2_SIMMATRIX
