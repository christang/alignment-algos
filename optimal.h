/**
 *  Package HMAP2.1
 *  File: optimal.h
 *  Desc: An enumerator for generating the optimal alignment. 
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_OPTIMAL
#define _HMAP2_OPTIMAL

#include "alib.h"
#include "alignment.h"
#include "enumerator.h"

#include <iostream>

template <class S1, class S2, class Etype>
class Optimal : public Enumerator<S1,S2,Etype> {
public:

  Optimal (align_t type = global);  

  inline int estimateSize () const { return 1; } ;
  void enumerate (DPMatrix<S1,S2,Etype>& dpm,
		  AlignmentSet<S1,S2,Etype>& as);
  void enumerate_local (DPMatrix<S1,S2,Etype>& dpm,
			AlignmentSet<S1,S2,Etype>& as);
  void find_max (const DPMatrix<S1,S2,Etype>& dpm,
		 int* q_last, int* t_last, float* score) const;

private:
  bool islocal;
};


template <class S1, class S2, class Etype>
Optimal<S1,S2,Etype>::Optimal (align_t type)
  : islocal (type==local) { }

// Perform standard traceback to extract optimal alignment

template <class S1, class S2, class Etype>
void Optimal<S1,S2,Etype>::enumerate (DPMatrix<S1,S2,Etype>& dpm,
				      AlignmentSet<S1,S2,Etype>& as) 
{
  if (islocal) {
    enumerate_local (dpm, as);
    return;
  }

  int k = as.size();
  as.resize(k+1);

  int q_last = dpm.getQuerySize() - 1;
  int t_last = dpm.getTemplateSize() - 1;
  
  as[k].score = dpm.getCell(q_last,t_last)->score;
  as[k].append(q_last,t_last);

  // Iterate traceback
  while (q_last>0) {
    const DPCell* cell = dpm.getCell(q_last,t_last);
    q_last = cell->prev_query_idx;
    t_last = cell->prev_template_idx;
    as[k].prepend(q_last,t_last);
  }

  // The first aligned pair should be (0,0)
  if (q_last!=0 || t_last!=0) throw string ("Illegal alignment start pair");
}

template <class S1, class S2, class Etype>
void Optimal<S1,S2,Etype>::
enumerate_local (DPMatrix<S1,S2,Etype>& dpm,
		 AlignmentSet<S1,S2,Etype>& as)
{
  int k = as.size();
  as.resize(k+1);

  int q_last = dpm.getQuerySize() - 1;
  int t_last = dpm.getTemplateSize() - 1;
  float score = 0.f;

  as[k].append(q_last,t_last);
  find_max (dpm, &q_last, &t_last, &score);

  as[k].score = score;
  as[k].prepend(q_last,t_last);

  // Iterate traceback
  while (q_last>0) {
    const DPCell* cell = dpm.getCell(q_last,t_last);
    q_last = cell->prev_query_idx;
    t_last = cell->prev_template_idx;
    if (dpm.getCell(q_last,t_last)->score <= 0.f) break;
    as[k].prepend(q_last,t_last);
  }

  if (q_last!=0 && t_last!=0) as[k].prepend (0,0);
} 

template <class S1, class S2, class Etype>
void Optimal<S1,S2,Etype>::find_max (const DPMatrix<S1,S2,Etype>& dpm,
				     int* q, int* t, float *s) const
{
  *q = dpm.getQuerySize()-2;
  *t = dpm.getTemplateSize()-2;
  *s = dpm.getCell(*q,*t)->score;

  for (int i=0; i<dpm.getQuerySize()-1; ++i) {
    for (int j=0; j<dpm.getTemplateSize()-1; ++j) {
      if (*s < dpm.getCell(i,j)->score) {
	*q = i;
	*t = j;
	*s = dpm.getCell(i,j)->score;
      }
    }
  }
}

#endif  // _HMAP2_OPTIMAL
