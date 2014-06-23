/**
 *  Package HMAP2.1
 *  File: optimal_rev.h
 *  Desc: An enumerator for generating the optimal reverse alignment. 
 *
 *  Created on 12/15/06
 *  Author: kuziemko @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_OPTIMAL_REV
#define _HMAP2_OPTIMAL_REV

#include "alib.h"
#include "alignment.h"
#include "enumerator.h"

#include <iostream>

template <class S1, class S2, class Etype>
class Optimal_Rev : public Enumerator<S1,S2,Etype> {
public:

  Optimal_Rev (align_t type = global);  

  inline int estimateSize () const { return 1; } ;
  void enumerate (const DPMatrix<S1,S2,Etype>& dpm,
		  AlignmentSet<S1,S2,Etype>& as) const;
  void enumerate_local (const DPMatrix<S1,S2,Etype>& dpm,
			AlignmentSet<S1,S2,Etype>& as) const;
  void find_max (const DPMatrix<S1,S2,Etype>& dpm,
		 int* q_last, int* t_last, float* score) const;

private:
  bool islocal;
};


template <class S1, class S2, class Etype>
Optimal_Rev<S1,S2,Etype>::Optimal_Rev (align_t type)
  : islocal (type==local) { }

// Perform standard traceback to extract optimal alignment

template <class S1, class S2, class Etype>
void Optimal_Rev<S1,S2,Etype>::enumerate (const DPMatrix<S1,S2,Etype>& dpm,
					  AlignmentSet<S1,S2,Etype>& as) const
{
  if (islocal) {
    enumerate_local (dpm, as);
    return;
  }

  int k = as.size();
  as.resize(k+1);

  int q_last = dpm.getQuerySize() - 1;
  int t_last = dpm.getTemplateSize() - 1;
  
  int q_first = 0;
  int t_first = 0;

  as[k].score = dpm.getCell(q_first,t_first)->score;
  as[k].append(q_first,t_first);

  // Iterate traceback
  while (q_first < q_last) {
    const DPCell* cell = dpm.getCell(q_first,t_first);
    q_first = cell->prev_query_idx;
    t_first = cell->prev_template_idx;
    as[k].append(q_first,t_first);
  }

  // The first aligned pair should be (0,0)
  if (q_first!=q_last || t_first!=t_last) throw string ("Illegal alignment start pair");
}

template <class S1, class S2, class Etype>
void Optimal_Rev<S1,S2,Etype>::
enumerate_local (const DPMatrix<S1,S2,Etype>& dpm,
		 AlignmentSet<S1,S2,Etype>& as) const
{
  int k = as.size();
  as.resize(k+1);

  int q_last = dpm.getQuerySize() - 1;
  int t_last = dpm.getTemplateSize() - 1;

  int q_first = 0;
  int t_first = 0;

  float score = 0.f;

  as[k].append(q_first,t_first);
  find_max (dpm, &q_first, &t_first, &score);

  as[k].score = score;
  as[k].append(q_first,t_first);

  // Iterate traceback
  while (q_first < q_last) {
    const DPCell* cell = dpm.getCell(q_first,t_first);
    q_first = cell->prev_query_idx;
    t_first = cell->prev_template_idx;
    if (dpm.getCell(q_first,t_first)->score <= 0.f) break;
    as[k].append(q_first,t_first);
  }

  if (q_first!=q_last && t_first!=t_last) as[k].append (q_last,t_last);
} 

template <class S1, class S2, class Etype>
void Optimal_Rev<S1,S2,Etype>::find_max (const DPMatrix<S1,S2,Etype>& dpm,
				     int* q, int* t, float *s) const
{
  *q = 0;
  *t = 0;
  *s = dpm.getCell(*q,*t)->score;

  for (int i=dpm.getQuerySize()-1; i>0; --i) {
    for (int j=dpm.getTemplateSize()-1; j>0; --j) {
      if (*s < dpm.getCell(i,j)->score) {
	*q = i;
	*t = j;
	*s = dpm.getCell(i,j)->score;
      }
    }
  }
}

#endif  // _HMAP2_OPTIMAL_REV
