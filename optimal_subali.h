/**
 *  Package HMAP2.1
 *  File: optimal_subali.h
 *  Desc: An enumerator for generating the optimal alignment. 
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_OPTIMAL_SUBALI
#define _HMAP2_OPTIMAL_SUBALI

#include "alib.h"
#include "alignment.h"
#include "enumerator.h"

#include <iostream>

template <class S1, class S2, class Etype>
class Optimal_Subali : public Enumerator<S1,S2,Etype> {
public:

  Optimal_Subali (int, int, int, int);  

  inline int estimateSize () const { return 1; } ;

  void enumerate( DPMatrix<S1,S2,Etype>& dpm,
		  AlignmentSet<S1,S2,Etype>& as );

private:
  bool islocal;

  int q1_end; // endpoints of sub-alignment region
  int t1_end;
  int q2_beg;
  int t2_beg;


};


template <class S1, class S2, class Etype>
Optimal_Subali<S1,S2,Etype>::Optimal_Subali (int q1, int t1, int q2, int t2)

  : q1_end (q1),
     t1_end (t1),
     q2_beg (q2),
     t2_beg (t2)
     
{ }



// Perform standard traceback to extract optimal alignment

template <class S1, class S2, class Etype>
void Optimal_Subali<S1,S2,Etype>::enumerate( DPMatrix<S1,S2,Etype>& dpm,
					     AlignmentSet<S1,S2,Etype>& as )
{
  int k = as.size();
  as.resize(k+1);

  int q_last = q2_beg;
  int t_last = t2_beg;
  
  as[k].score = dpm.getCell(q_last,t_last)->score;
  as[k].append(q_last,t_last);

  // Iterate traceback
  while (q_last>q1_end) {
    const DPCell* cell = dpm.getCell(q_last,t_last);
    q_last = cell->prev_query_idx;
    t_last = cell->prev_template_idx;
    as[k].prepend(q_last,t_last);
  }

  // The first aligned pair should be (q1_end, t1_end)
  if (q_last!=q1_end || t_last!=t1_end ) throw string ("Illegal alignment start pair");

}


#endif  // _HMAP2_OPTIMAL_SUBALI
