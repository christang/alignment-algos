/**
 *  Package HMAP2.1
 *  File: dpmatrix.h
 *  Desc: Template class for dynamic programming matrix.
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_DPMATRIX
#define _HMAP2_DPMATRIX

#include <iostream>

#include "alib.h"
#include "evaluator.h"
#include "matrix.h"
#include "simmatrix.h"

enum direction_t {
  fwd = 1,
  rev = 2
};

struct DPCell {
  int   prev_query_idx;
  int   prev_template_idx;
  int   query_idx;
  int   template_idx;
  float score;
  static const int null;
  DPCell();
  void setTB (int, int, float);
};

template <class S1, class S2, class Etype>
class DPMatrix {

public:
  
  DPMatrix (const S1& query_seq,
	    const S2& templ_seq  );
  DPMatrix (const S1& query_seq,
	    const S2& templ_seq,
	    const Evaluator<S1,S2,Etype>& eval,
	    direction_t direction = fwd,
	    align_t type = global);

  DPMatrix (const S1& query_seq,
	    const S2& templ_seq,
	    const Evaluator<S1,S2,Etype>& eval,
	    int q1_end, int q2_beg,
	    int t1_end, int t2_beg,
	    direction_t dir = fwd,
	    align_t type = global);

  ~DPMatrix ();

  void setEvaluator (const Evaluator<S1,S2,Etype>& eval,
		     direction_t direction  );
  const DPCell* getCell (int query_pos,
			 int templ_pos  ) const;
  direction_t getDirection () const;
  int getQuerySize () const;
  int getTemplateSize () const;
  const Evaluator<S1,S2,Etype>* getEvaluator () const;
  inline const S1*    getQuerySequence () const { return query_seq; } ;
  inline const S2* getTemplateSequence () const { return templ_seq; } ;
  inline float getSim (int i, int j) const 
  { return simmatrix->operator()(i,j); } ;
  void reevaluate ();

protected:

  template <class S> void saveSeq (const S& source, const S** target  );
  void saveEvaluator (const Evaluator<S1,S2,Etype>& source,
		      const Evaluator<S1,S2,Etype>** target  );

  void initMtxMem ();
  void initMtxVal ();
  void resetMtxVal ();

  void build ();

  //AK
  void build_subdpm(int, int, int, int);

  void build_forw_dpm_nonlinear_gaps (int q0, int q1,
				      int t0, int t1  );
  void build_forw_local_dpm_nonlinear_gaps (int q0, int q1,
					    int t0, int t1  );
  
  void build_rev_dpm_nonlinear_gaps (int q0, int q1,
				     int t0, int t1  );
  void build_rev_local_dpm_nonlinear_gaps (int q0, int q1,
					   int t0, int t1  );

  void build_forw_dpm_linear_gaps ();
  void build_rev_dpm_linear_gaps ();
  
  const S1* query_seq;
  const S2* templ_seq;
  const Evaluator<S1,S2,Etype>* evaluator;
  direction_t direction;
  bool islocal;

  matrix<DPCell>* dpmatrix;
  SimilarityMatrix* simmatrix;

};


template <class S1, class S2, class Etype>
ostream& operator<< (ostream& o, const DPMatrix<S1,S2,Etype>& dpm) 
{
  int q_len = dpm.getQuerySequence()->size();
  int t_len = dpm.getTemplateSequence()->size();

  for (int i=0; i<q_len; ++i) {
    for (int j=0; j<t_len; ++j) {
      o << dpm.getCell(i,j)->score << "\t";
    }
    o << endl;
  }
  return o;
}


template <class S1, class S2, class Etype>
DPMatrix<S1,S2,Etype>::DPMatrix (const S1& qs,
				 const S2& ts  )
  : evaluator (0),
    direction (fwd),
    islocal   (false),
    dpmatrix  (0),
    simmatrix (0)
{
  saveSeq<S1> (qs,query_seq);
  saveSeq<S2> (ts,templ_seq);
  initMtxMem ();
  initMtxVal ();
}

template <class S1, class S2, class Etype>
DPMatrix<S1,S2,Etype>::DPMatrix (const S1& qs,
				 const S2& ts,
				 const Evaluator<S1,S2,Etype>& eval,
				 direction_t dir,
				 align_t type  )
  : evaluator (0),
    direction (dir),
    islocal   (type==local),
    dpmatrix  (0),
    simmatrix (0)
{
  saveSeq<S1> (qs,&query_seq);
  saveSeq<S2> (ts,&templ_seq);
  saveEvaluator (eval,&evaluator);
  initMtxMem ();
  initMtxVal ();
  build ();
}


template <class S1, class S2, class Etype>
DPMatrix<S1,S2,Etype>::DPMatrix (const S1& qs,
				 const S2& ts,
				 const Evaluator<S1,S2,Etype>& eval,
				 int q1_end, int t1_end,
				 int q2_beg, int t2_beg,
				 direction_t dir,
				 align_t type )

  : evaluator (0),
    direction (dir),
    islocal   (type==local),
    dpmatrix  (0),
    simmatrix (0)
{
  saveSeq<S1> (qs,&query_seq);
  saveSeq<S2> (ts,&templ_seq);
  saveEvaluator (eval,&evaluator);
  initMtxMem ();
  initMtxVal ();
  build_subdpm ( q1_end, t1_end, q2_beg, t2_beg );
}



template <class S1, class S2, class Etype>
DPMatrix<S1,S2,Etype>::~DPMatrix ()
{
  //delete query_seq;
  //delete templ_seq;
  //if (evaluator) delete evaluator;
  if (dpmatrix)  delete dpmatrix;
  if (simmatrix) delete simmatrix;
}

template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::setEvaluator (const Evaluator<S1,S2,Etype>& eval,
					  direction_t dir)
{
  direction = dir;
  //if (evaluator) delete evaluator;
  saveEvaluator (eval,&evaluator);
  reevaluate ();
}

template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::reevaluate ()
{
  resetMtxVal ();
  build ();
}

template <class S1, class S2, class Etype> template <class S>
void DPMatrix<S1,S2,Etype>::saveSeq (const S& source, const S** target  )
{ *target = &source; }

template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::saveEvaluator(const Evaluator<S1,S2,Etype>& source,
					  const Evaluator<S1,S2,Etype>** target)
{ *target = &source; }

template <class S1, class S2, class Etype>
inline const DPCell* DPMatrix<S1,S2,Etype>::getCell(int query_pos,
						    int templ_pos  ) const
{ return &((*dpmatrix)[query_pos][templ_pos]); }

template <class S1, class S2, class Etype>
direction_t DPMatrix<S1,S2,Etype>::getDirection () const
{ return direction; }

template <class S1, class S2, class Etype>
const Evaluator<S1,S2,Etype>* DPMatrix<S1,S2,Etype>::getEvaluator () const
{ return evaluator; }

template <class S1, class S2, class Etype>
int DPMatrix<S1,S2,Etype>::getQuerySize () const
{ return (int)query_seq->size(); }

template <class S1, class S2, class Etype>
int DPMatrix<S1,S2,Etype>::getTemplateSize () const
{ return (int)templ_seq->size(); }

template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::initMtxMem () 
{
  if (dpmatrix) delete dpmatrix;
  
  int sz1 = (int)query_seq->size();  // row size
  int sz2 = (int)templ_seq->size();  // col size

  dpmatrix = new matrix<DPCell> (sz1,sz2);
}

template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::initMtxVal ()
{
  int sz1 = (int)query_seq->size();
  int sz2 = (int)templ_seq->size();

  for (int i=0; i<sz1; ++i) {
    for (int j=0; j<sz2; ++j) {
      (*dpmatrix)[i][j].query_idx = i;
      (*dpmatrix)[i][j].template_idx = j;
    }
  }
}

template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::resetMtxVal ()
{
  int sz1 = query_seq->size();
  int sz2 = templ_seq->size();
  
  for (int i=0; i<sz1; ++i) {
    for (int j=0; j<sz2; ++j) {
      (*dpmatrix)[i][j].prev_query_idx = DPCell::null;
      (*dpmatrix)[i][j].prev_template_idx = DPCell::null;
      (*dpmatrix)[i][j].score = 0;
    }
  }
}


template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::build ()
{
  if( simmatrix ) {
    delete simmatrix;
  }

  evaluator->pre_calculate (*query_seq, *templ_seq);
  simmatrix = new SimilarityMatrix (*query_seq,
				    *templ_seq,
				    *evaluator);

  int sz1 = (int)query_seq->size();
  int sz2 = (int)templ_seq->size();

  (*dpmatrix)[0][0].score = 0;
  (*dpmatrix)[sz1-1][sz2-1].score = 0;
    
  if (direction==fwd) {
    if (islocal) build_forw_local_dpm_nonlinear_gaps (0,sz1-1,0,sz2-1);
    else build_forw_dpm_nonlinear_gaps (0,sz1-1,0,sz2-1);
  } else {
    if (islocal) build_rev_local_dpm_nonlinear_gaps (0,sz1-1,0,sz2-1);
    else build_rev_dpm_nonlinear_gaps (0,sz1-1,0,sz2-1);
  }
  
}

template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::build_subdpm( int q1_end, int t1_end,
					   int q2_beg, int t2_beg )
{

  if( simmatrix ) {
    delete simmatrix;
  }

  evaluator->pre_calculate (*query_seq, *templ_seq);
  simmatrix = new SimilarityMatrix (*query_seq,
				    *templ_seq,
				    *evaluator);

  (*dpmatrix)[q1_end][t1_end].score = 0;
  (*dpmatrix)[q2_beg][t2_beg].score = 0;

  if (direction==fwd) {
    if (islocal) {
      build_forw_local_dpm_nonlinear_gaps (q1_end,q2_beg,t1_end,t2_beg);
    }
    else {
      build_forw_dpm_nonlinear_gaps (q1_end,q2_beg,t1_end,t2_beg);
    }
  }
  else {
    if (islocal) {
      build_rev_local_dpm_nonlinear_gaps (q1_end,q2_beg,t1_end,t2_beg);
    }
    else {
      build_rev_dpm_nonlinear_gaps (q1_end,q2_beg,t1_end,t2_beg);
    }
  }

}


template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::
build_forw_dpm_nonlinear_gaps (int q0, int q1, int t0, int t1)
{
  if (q1<=q0 || t1<=t0) 
    throw string ("Illegal bounds building DPM");

  float s_initial = (*dpmatrix)[q0][t0].score;

  int q0_p1 = q0 + 1;
  int t0_p1 = t0 + 1;
  int q0_p2 = q0 + 2;
  int t0_p2 = t0 + 2;
  int q1_m1 = q1 - 1;
  int t1_m1 = t1 - 1;
  
  float s;
  
  // Special case #1: boundary conditions force deletion 
  if (q1==q0_p1) {
    s = s_initial;
    s -= evaluator->deletion(*query_seq,*templ_seq,q0,q1,t0,t1);
    s += (*simmatrix)[q1][t1];
    (*dpmatrix)[q1][t1].setTB(q0,t0,s);
    return;
  }

  // Special case #2: boundary conditions force insertion
  if (t1==t0_p1) {
    s = s_initial;
    s -= evaluator->insertion(*query_seq,*templ_seq,q0,q1,t0,t1);
    s += (*simmatrix)[q1][t1];
    (*dpmatrix)[q1][t1].setTB(q0,t0,s);
    return;
  }

  // Dynamic programming in the standard case:
  // Start by evaluating the boundary cells (*s), i.e.:

  //        |t0|       |t1|  
  //    ----+--+-------+--+----
  //      q0|M |       |  |
  //    ----+--+-------+--+----
  //        |  |*******|  |  
  //        |  |*      |  |  
  //        |  |*      |  |  
  //        |  |*      |  |  
  //    ----+--+-------+--+----
  //      q1|  |       |  |
  //    ----+--+-------+--+----
  //        |  |       |  |  

  // Match state @ (+1,+1)
  s = s_initial + (*simmatrix)[q0_p1][t0_p1];
  (*dpmatrix)[q0_p1][t0_p1].setTB(q0,t0,s);

  // Top row boundary
  for (int j=t0_p2; j<t1; ++j) {
    s = s_initial;
    s -= evaluator->deletion(*query_seq,*templ_seq,q0,q0_p1,t0,j);
    s += (*simmatrix)[q0_p1][j];
    (*dpmatrix)[q0_p1][j].setTB(q0,t0,s);
  }

  // Left col boundary
  for (int i=q0_p2; i<q1; ++i) {
    s = s_initial;
    s -= evaluator->insertion(*query_seq,*templ_seq,q0,i,t0,t0_p1);
    s += (*simmatrix)[i][t0_p1];
    (*dpmatrix)[i][t0_p1].setTB(q0,t0,s);
   }

  int opt_i;
  int opt_j;
  float opt_s;

  // Evaluate non-boundary cells

  //        |t0|       |t1|  
  //    ----+--+-------+--+----
  //      q0|M |       |  |
  //    ----+--+-------+--+----
  //        |  |\_     |  |  
  //        |  |  \_   |  |  
  //        |  |    *  |  |  
  //        |  |       |  |  
  //    ----+--+-------+--+----
  //      q1|  |       |  |
  //    ----+--+-------+--+----
  //        |  |       |  |  

  for (int i=q0_p2; i<q1; ++i) {
    for (int j=t0_p2; j<t1; ++j) {

      int i_m1 = i-1;
      int j_m1 = j-1;
      
      // Match case
      opt_i = i_m1;
      opt_j = j_m1;
      opt_s = (*dpmatrix)[opt_i][opt_j].score+(*simmatrix)[i][j];
      
      // Deletion cases
      for (int k=t0_p1; k<j_m1; ++k) {
	s = (*dpmatrix)[i_m1][k].score;
	s -= evaluator->deletion(*query_seq,*templ_seq,i_m1,i,k,j);
	s += (*simmatrix)[i][j];
	if (s > opt_s) {
	  opt_i = i_m1;
	  opt_j = k;
	  opt_s = s;
	}
      }

      // Insertion cases
      for (int k=q0_p1; k<i_m1; ++k) {
	s = (*dpmatrix)[k][j_m1].score;
	s -= evaluator->insertion(*query_seq,*templ_seq,k,i,j_m1,j);
	s += (*simmatrix)[i][j];
	if (s > opt_s) {
	  opt_i = k;
	  opt_j = j_m1;
	  opt_s = s;
	}
      }

      // Set optimal TB
      (*dpmatrix)[i][j].setTB (opt_i, opt_j, opt_s);

    }
  }
  
  // Finally, find the optimal path, assuming that
  // the pair (q1,t1) is a match

  //        |t0|       |t1|  
  //    ----+--+-------+--+----
  //      q0|M |       |  |
  //    ----+--+-------+--+----
  //        |  |\_____?|  |  
  //        |  | |\_  \|  |  
  //        |  | |  \ ?|  |  
  //        |  |?\???\?|  |  
  //    ----+--+-------+--+----
  //      q1|  |       |M*|
  //    ----+--+-------+--+----
  //        |  |       |  |  

  // Match state @ (q1,t1)
  opt_i = q1_m1;
  opt_j = t1_m1;
  opt_s = (*dpmatrix)[opt_i][opt_j].score+(*simmatrix)[q1][t1];
      
  // Bottom row boundary
  for (int k=t0_p1; k<t1; ++k) {
    s = (*dpmatrix)[q1_m1][k].score;
    s -= evaluator->deletion(*query_seq,*templ_seq,q1_m1,q1,k,t1);
    s += (*simmatrix)[q1][t1];
    if (s > opt_s) {
      opt_i = q1_m1;
      opt_j = k;
      opt_s = s;
    }
  }
  
  // Right col boundary
  for (int k=q0_p1; k<q1; ++k) {
    s = (*dpmatrix)[k][t1_m1].score;
    s -= evaluator->insertion(*query_seq,*templ_seq,k,q1,t1_m1,t1);
    s += (*simmatrix)[q1][t1];
    if (s > opt_s) {
      opt_i = k;
      opt_j = t1_m1;
      opt_s = s;
    }
  }
  
  // Set optimal TB
  (*dpmatrix)[q1][t1].setTB (opt_i, opt_j, opt_s);

}

template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::
build_forw_local_dpm_nonlinear_gaps (int q0, int q1,int t0, int t1)
{

  if (q1<=q0 || t1<=t0) 
    throw string ("Illegal bounds building DPM");

  float s_initial = (*dpmatrix)[q0][t0].score;

  int q0_p1 = q0 + 1;
  int t0_p1 = t0 + 1;
  int q0_p2 = q0 + 2;
  int t0_p2 = t0 + 2;
  int q1_m1 = q1 - 1;
  int t1_m1 = t1 - 1;
  
  float s;
  
  // Special case #1: boundary conditions force deletion 
  if (q1==q0_p1) {
    s = s_initial;
    s -= evaluator->deletion(*query_seq,*templ_seq,q0,q1,t0,t1);
    s += (*simmatrix)[q1][t1];
    (*dpmatrix)[q1][t1].setTB(q0,t0,s);
    return;
  }

  // Special case #2: boundary conditions force insertion
  if (t1==t0_p1) {
    s = s_initial;
    s -= evaluator->insertion(*query_seq,*templ_seq,q0,q1,t0,t1);
    s += (*simmatrix)[q1][t1];
    (*dpmatrix)[q1][t1].setTB(q0,t0,s);
    return;
  }

  // Dynamic programming in the standard case:
  // Start by evaluating the boundary cells (*s), i.e.:

  // Match state @ (+1,+1)
  s = s_initial + (*simmatrix)[q0_p1][t0_p1];
  s = max (0.f,s);
  (*dpmatrix)[q0_p1][t0_p1].setTB(q0,t0,s);

  // Top row boundary
  for (int j=t0_p2; j<t1; ++j) {
    s = s_initial;
    s -= evaluator->deletion(*query_seq,*templ_seq,q0,q0_p1,t0,j);
    s += (*simmatrix)[q0_p1][j];
    s = max (0.f,s);
    (*dpmatrix)[q0_p1][j].setTB(q0,t0,s);
  }

  // Left col boundary
  for (int i=q0_p2; i<q1; ++i) {
    s = s_initial;
    s -= evaluator->insertion(*query_seq,*templ_seq,q0,i,t0,t0_p1);
    s += (*simmatrix)[i][t0_p1];
    s = max (0.f,s);
    (*dpmatrix)[i][t0_p1].setTB(q0,t0,s);
   }

  int opt_i;
  int opt_j;
  float opt_s;

  // Evaluate non-boundary cells

  for (int i=q0_p2; i<q1; ++i) {
    for (int j=t0_p2; j<t1; ++j) {

      int i_m1 = i-1;
      int j_m1 = j-1;
      
      // Match case
      opt_i = i_m1;
      opt_j = j_m1;
      opt_s = (*dpmatrix)[opt_i][opt_j].score+(*simmatrix)[i][j];
      opt_s = max (0.f,opt_s);
   
      // Deletion cases
      for (int k=t0_p1; k<j_m1; ++k) {
	s = (*dpmatrix)[i_m1][k].score;
	s -= evaluator->deletion(*query_seq,*templ_seq,i_m1,i,k,j);
	s += (*simmatrix)[i][j];
	s = max (0.f,s);
	if (s > opt_s) {
	  opt_i = i_m1;
	  opt_j = k;
	  opt_s = s;
	}
      }

      // Insertion cases
      for (int k=q0_p1; k<i_m1; ++k) {
	s = (*dpmatrix)[k][j_m1].score;
	s -= evaluator->insertion(*query_seq,*templ_seq,k,i,j_m1,j);
	s += (*simmatrix)[i][j];
	s = max (0.f,s);
	if (s > opt_s) {
	  opt_i = k;
	  opt_j = j_m1;
	  opt_s = s;
	}
      }

      // Set optimal TB
      (*dpmatrix)[i][j].setTB (opt_i, opt_j, opt_s);

    }
  }
  
  // Finally, find the optimal path, assuming that
  // the pair (q1,t1) is a match

  // Match state @ (q1,t1)
  opt_i = q1_m1;
  opt_j = t1_m1;
  opt_s = (*dpmatrix)[opt_i][opt_j].score+(*simmatrix)[q1][t1];
  opt_s = max (0.f,opt_s);
      
  // Bottom row boundary
  for (int k=t0_p1; k<t1; ++k) {
    s = (*dpmatrix)[q1_m1][k].score;
    s -= evaluator->deletion(*query_seq,*templ_seq,q1_m1,q1,k,t1);
    s += (*simmatrix)[q1][t1];
    s = max (0.f,s);
    if (s > opt_s) {
      opt_i = q1_m1;
      opt_j = k;
      opt_s = s;
    }
  }
  
  // Right col boundary
  for (int k=q0_p1; k<q1; ++k) {
    s = (*dpmatrix)[k][t1_m1].score;
    s -= evaluator->insertion(*query_seq,*templ_seq,k,q1,t1_m1,t1);
    s += (*simmatrix)[q1][t1];
    s = max (0.f,s);
    if (s > opt_s) {
      opt_i = k;
      opt_j = t1_m1;
      opt_s = s;
    }
  }
  
  // Set optimal TB
  (*dpmatrix)[q1][t1].setTB (opt_i, opt_j, opt_s);

}

template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::build_rev_dpm_nonlinear_gaps (int q0, int q1,
							  int t0, int t1)
{

  cerr << "starting to build rev non-local" << endl;

  if (q1<=q0 || t1<=t0) 
    throw string ("Illegal bounds building DPM");

  float s_initial = (*dpmatrix)[q1][t1].score;

  int q0_p1 = q0 + 1;
  int t0_p1 = t0 + 1;
  int q1_m2 = q1 - 2;
  int t1_m2 = t1 - 2;
  int q1_m1 = q1 - 1;
  int t1_m1 = t1 - 1;
  
  float s;

  // Special case #1: boundary conditions force deletion 
  if (q1==q0_p1) {
    s = s_initial;
    s -= evaluator->deletion(*query_seq,*templ_seq,q0,q1,t0,t1);
    s += (*simmatrix)[q0][t0];
    (*dpmatrix)[q0][t0].setTB(q1,t1,s);
    return;
  }

  // Special case #2: boundary conditions force insertion
  if (t1==t0_p1) {
    s = s_initial;
    s -= evaluator->insertion(*query_seq,*templ_seq,q0,q1,t0,t1);
    s += (*simmatrix)[q0][t0];
    (*dpmatrix)[q0][t0].setTB(q1,t1,s);
    return;
  }

  // Dynamic programming in the standard case:
  // Start by evaluating the boundary cells (*'s), i.e.:

  //        |t0|       |t1|  
  //    ----+--+-------+--+----
  //      q0|  |       |  |
  //    ----+--+-------+--+----
  //        |  |      *|  |  
  //        |  |      *|  |  
  //        |  |      *|  |  
  //        |  |*******|  |  
  //    ----+--+-------+--+----
  //      q1|  |       |M |
  //    ----+--+-------+--+----
  //        |  |       |  |  

  // Match state @ (+1,+1)
  s = s_initial + (*simmatrix)[q1_m1][t1_m1];
  (*dpmatrix)[q1_m1][t1_m1].setTB(q1,t1,s);

  // Bottom row boundary
  for (int j=t1_m2; j>t0; --j) {
    s = s_initial;
    s -= evaluator->deletion(*query_seq,*templ_seq,q1_m1,q1,j,t1);
    s += (*simmatrix)[q1_m1][j];
    (*dpmatrix)[q1_m1][j].setTB(q1,t1,s);
  }

  // Right col boundary
  for (int i=q1_m2; i>q0; --i) {
    s = s_initial;
    s -= evaluator->insertion(*query_seq,*templ_seq,i,q1,t1_m1,t1);
    s += (*simmatrix)[i][t1_m1];
    (*dpmatrix)[i][t1_m1].setTB(q1,t1,s);
  }


  int opt_i;
  int opt_j;
  float opt_s;

  // Evaluate non-boundary cells

  //        |t0|       |t1|  
  //    ----+--+-------+--+----
  //      q0|  |       |  |
  //    ----+--+-------+--+----
  //        |  |       |  |  
  //        |  |   *   |  |  
  //        |  |    \_ |  |  
  //        |  |      \|  |  
  //    ----+--+-------+--+----
  //      q1|  |       |M |
  //    ----+--+-------+--+----
  //        |  |       |  |  

  for (int i=q1_m2; i>q0; --i) {
    for (int j=t1_m2; j>t0; --j) {

      int i_p1 = i+1;
      int j_p1 = j+1;
      
      // Match case
      opt_i = i_p1;
      opt_j = j_p1;
      opt_s = (*dpmatrix)[opt_i][opt_j].score+(*simmatrix)[i][j];
      
      // Deletion cases
      for (int k=t1_m1; k>j_p1; --k) {
	s = (*dpmatrix)[i_p1][k].score;
	s -= evaluator->deletion(*query_seq,*templ_seq,i,i_p1,j,k);
	s += (*simmatrix)[i][j];
	if (s > opt_s) {
	  opt_i = i_p1;
	  opt_j = k;
	  opt_s = s;
	}
      }

      // Insertion cases
      for (int k=q1_m1; k>i_p1; --k) {
	s = (*dpmatrix)[k][j_p1].score;
	s -= evaluator->insertion(*query_seq,*templ_seq,i,k,j,j_p1);
	s += (*simmatrix)[i][j];
	if (s > opt_s) {
	  opt_i = k;
	  opt_j = j_p1;
	  opt_s = s;
	}
      }

      // Set optimal TB
      (*dpmatrix)[i][j].setTB (opt_i, opt_j, opt_s);

    }
  }

  // Finally, find the optimal last step, assuming that
  // the pair (q0,t0) is a match

  //        |t0|       |t1|  
  //    ----+--+-------+--+----
  //      q0|M*|       |  |
  //    ----+--+-------+--+----
  //        |  | _     |  |  
  //        |  | |\_   |  |  
  //        |  | |  \_ |  |  
  //        |  | \    \|  |  
  //    ----+--+-------+--+----
  //      q1|  |       |M |
  //    ----+--+-------+--+----
  //        |  |       |  |  


  // Match state @ (q0,t0)
  opt_i = q0_p1;
  opt_j = t0_p1;
  opt_s = (*dpmatrix)[opt_i][opt_j].score+(*simmatrix)[q0][t0];
      
  // Top row boundary
  for (int k=t1_m1; k>t0; --k) {
    s = (*dpmatrix)[q0_p1][k].score;
    s -= evaluator->deletion(*query_seq,*templ_seq,q0,q0_p1,t0,k);
    s += (*simmatrix)[q0][t0];
    if (s > opt_s) {
      opt_i = q0_p1;
      opt_j = k;
      opt_s = s;
    }
  }
  
  // Left col boundary
  for (int k=q1_m1; k>q0; --k) {
    s = (*dpmatrix)[k][t0_p1].score;
    s -= evaluator->insertion(*query_seq,*templ_seq,q0,k,t0,t0_p1);
    s += (*simmatrix)[q0][t0];
    if (s > opt_s) {
      opt_i = k;
      opt_j = t1_m1;
      opt_s = s;
    }
  }
  
  // Set optimal TB
  (*dpmatrix)[q0][t0].setTB (opt_i, opt_j, opt_s);

  
}

template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::
build_rev_local_dpm_nonlinear_gaps (int q0, int q1,int t0, int t1)
{

  if (q1<=q0 || t1<=t0) 
    throw string ("Illegal bounds building DPM");

  float s_initial = (*dpmatrix)[q1][t1].score;

  int q0_p1 = q0 + 1;
  int t0_p1 = t0 + 1;
  int q1_m2 = q1 - 2;
  int t1_m2 = t1 - 2;
  int q1_m1 = q1 - 1;
  int t1_m1 = t1 - 1;
  
  float s;
  
  // Special case #1: boundary conditions force deletion 
  if (q1==q0_p1) {
    s = s_initial;
    s -= evaluator->deletion(*query_seq,*templ_seq,q0,q1,t0,t1);
    s += (*simmatrix)[q0][t0];
    (*dpmatrix)[q0][t0].setTB(q1,t1,s);
    return;
  }

  // Special case #2: boundary conditions force insertion
  if (t1==t0_p1) {
    s = s_initial;
    s -= evaluator->insertion(*query_seq,*templ_seq,q0,q1,t0,t1);
    s += (*simmatrix)[q0][t0];
    (*dpmatrix)[q0][t0].setTB(q1,t1,s);
    return;
  }

  // Dynamic programming in the standard case:
  // Start by evaluating the boundary cells (*s), i.e.:

  // Match state @ (+1,+1)
  s = s_initial + (*simmatrix)[q1_m1][t1_m1];
  s = max (0.f,s);
  (*dpmatrix)[q1_m1][t1_m1].setTB(q1,t1,s);

  // Top row boundary
  for (int j=t1_m2; j>t0; --j) {
    s = s_initial;
    s -= evaluator->deletion(*query_seq,*templ_seq,q1_m1,q1,j,t1);
    s += (*simmatrix)[q1_m1][j];
    s = max (0.f,s);
    (*dpmatrix)[q1_m1][j].setTB(q1,t1,s);
  }

  // Left col boundary
  for (int i=q1_m2; i>q0; --i) {
    s = s_initial;
    s -= evaluator->insertion(*query_seq,*templ_seq,i,q1,t1_m1,t1);
    s += (*simmatrix)[i][t1_m1];
    s = max (0.f,s);
    (*dpmatrix)[i][t1_m1].setTB(q1,t1,s);
   }

  int opt_i;
  int opt_j;
  float opt_s;

  // Evaluate non-boundary cells

  for (int i=q1_m2; i>q0; --i) {
    for (int j=t1_m2; j>t0; --j) {

      int i_p1 = i+1;
      int j_p1 = j+1;
      
      // Match case
      opt_i = i_p1;
      opt_j = j_p1;
      opt_s = (*dpmatrix)[opt_i][opt_j].score+(*simmatrix)[i][j];
      opt_s = max (0.f,opt_s);
   
      // Deletion cases
      for (int k=t1_m1; k>j_p1; --k) {
	s = (*dpmatrix)[i_p1][k].score;
	s -= evaluator->deletion(*query_seq,*templ_seq,i,i_p1,j,k);
	s += (*simmatrix)[i][j];
	s = max (0.f,s);
	if (s > opt_s) {
	  opt_i = i_p1;
	  opt_j = k;
	  opt_s = s;
	}
      }

      // Insertion cases
      for (int k=q1_m1; k>i_p1; --k) {
	s = (*dpmatrix)[k][j_p1].score;
	s -= evaluator->insertion(*query_seq,*templ_seq,i,k,j,j_p1);
	s += (*simmatrix)[i][j];
	s = max (0.f,s);
	if (s > opt_s) {
	  opt_i = k;
	  opt_j = j_p1;
	  opt_s = s;
	}
      }

      // Set optimal TB
      (*dpmatrix)[i][j].setTB (opt_i, opt_j, opt_s);

    }
  }
  
  // Finally, find the optimal path, assuming that
  // the pair (q0,t0) is a match

  // Match state @ (q0,t0)
  opt_i = q0_p1;
  opt_j = t0_p1;
  opt_s = (*dpmatrix)[opt_i][opt_j].score+(*simmatrix)[q0][t0];
  opt_s = max (0.f,opt_s);
      
  // Bottom row boundary
  for (int k=t1_m1; k>t0; --k) {
    s = (*dpmatrix)[q0_p1][k].score;
    s -= evaluator->deletion(*query_seq,*templ_seq,q0,q0_p1,t0,k);
    s += (*simmatrix)[q0][t0];
    s = max (0.f,s);
    if (s > opt_s) {
      opt_i = q0_p1;
      opt_j = k;
      opt_s = s;
    }
  }
  
  // Right col boundary
  for (int k=q1_m1; k>q0; --k) {
    s = (*dpmatrix)[k][t0_p1].score;
    s -= evaluator->insertion(*query_seq,*templ_seq,q0,k,t0,t0_p1);
    s += (*simmatrix)[q0][t0];
    s = max (0.f,s);
    if (s > opt_s) {
      opt_i = k;
      opt_j = t0_p1;
      opt_s = s;
    }
  }
  
  // Set optimal TB
  (*dpmatrix)[q0][t0].setTB (opt_i, opt_j, opt_s);
  
}

template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::build_forw_dpm_linear_gaps ()
{
  throw string ("Under construction!~Please use nonlinear gap algorithm");
}

template <class S1, class S2, class Etype>
void DPMatrix<S1,S2,Etype>::build_rev_dpm_linear_gaps ()
{
  throw string ("Under construction!~Please use nonlinear gap algorithm");
}

#endif  // _HMAP2_DPMATRIX
