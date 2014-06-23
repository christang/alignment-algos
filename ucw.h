/**
 *  Package HMAP2.1
 *  File: ucw.h
 *  Desc: Unconstrained near-optimal alignment enumeration
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_UCW
#define _HMAP2_UCW

#include <string>

#include "alignment.h"
#include "enumerator.h"
#include "noalib.h"
#include "sflags.h"

using namespace std;

template <class S1, class S2, class Etype>
class UnconstrainedNearOptimal : public Enumerator<S1,S2,Etype> {
public:
  
  typedef AlignedPairList<S1,S2> SingleAlignment;
  typedef AlignedPair<S1,S2>     SinglePair;

  UnconstrainedNearOptimal (const NOaliParams& p) ;

  unsigned int user_limit;

  inline int estimateSize () const { return params->number_suboptimal; } ;

  void enumerate (DPMatrix<S1,S2,Etype>& dpm,
		  AlignmentSet<S1,S2,Etype>& as);

private:

  void branch (const DPMatrix<S1,S2,Etype>& dpm,
	       AlignmentSet<S1,S2,Etype>& as,
	       int q0, int t0, int k0, float thresh);

  void opt_path (const DPMatrix<S1,S2,Etype>& dpm,
		 AlignmentSet<S1,S2,Etype>& as,
		 int q0, int t0, int k0, float thresh);

  const NOaliParams* params;
  bool warn_user;

};

template <class S1, class S2, class Etype>
UnconstrainedNearOptimal<S1,S2,Etype>::
UnconstrainedNearOptimal (const NOaliParams& p) 
  : params(&p), warn_user(true) {}

// Perform constrained Waterman-style branching tracebacks

template <class S1, class S2, class Etype>
void UnconstrainedNearOptimal<S1,S2,Etype>::
enumerate (DPMatrix<S1,S2,Etype>& dpm,
	   AlignmentSet<S1,S2,Etype>& as)
{
  int q_last = dpm.getQuerySize() - 1;
  int t_last = dpm.getTemplateSize() - 1;

  warn_user = true;
  user_limit = 100000;  
                         // This should be a parameter, but I'll 
                         // shortcut this for now.

  as.reserve(min(estimateSize()*20,int(user_limit)));

  as.push_back(SingleAlignment());
  int k_last = as.size() - 1;

  float threshold = 
    (1.f - params->delta_ratio) * dpm.getCell(q_last,t_last)->score;
  threshold = min(threshold,dpm.getCell(q_last,t_last)->score-0.1f);
  branch (dpm, as, q_last, t_last, k_last, threshold);
  as.sortSet (params->number_suboptimal);
}

template <class S1, class S2, class Etype>
void UnconstrainedNearOptimal<S1,S2,Etype>::
branch (const DPMatrix<S1,S2,Etype>& dpm,
	AlignmentSet<S1,S2,Etype>& as,
	int q0, int t0, int k0, float threshold)
{
  if (q0==1 || t0==1) {

    // This is the base case.  If either q0 or t0 reaches 1, we
    // terminate the recursion and complete the k0'th alignment.

    as[k0].prepend (q0,t0);
    as[k0].prepend (0,0);
    as[k0].score += dpm.getCell(q0,t0)->score;
    
  } else {

    int k = k0;
    float f,r,g;
    SingleAlignment curr (as[k0]);
    
    // First, check if user limits on number of alignments was exceeded.
    // Force the optimal alignment to the end if true.
    
    if (as.size() > user_limit) {
      if (warn_user) {
	warn_user = false;
	cerr << "Number of alignments exceeding user limits (";
	cerr << user_limit << ")" << endl;
	cerr << "User is advised to lower current threshold (";
	cerr << params->delta_ratio  << ") or further constrain alignments";
	cerr << endl << endl;
      }
      opt_path (dpm,as,q0,t0,k0,threshold);
      return;
    }
    
    // Check if more space is needed.
    
    if (as.size() > 0.9 * as.capacity() && as.capacity() < user_limit)
      as.reserve (min<unsigned int>(as.capacity()*2,user_limit));
    
    // Begin the recursive step.

    const S1* qseq = dpm.getQuerySequence();
    const S2* tseq = dpm.getTemplateSequence();
    const Evaluator<S1,S2,Etype>* e = dpm.getEvaluator();
    
    // Match case.  f and r are used to store the forward and reverse
    // traceback scores.  Their sum must exceed threshold to traverse
    // the branch.  This is Waterman's condition for near-optimal 
    // alignments.
    
    r = curr.score + dpm.getSim (q0,t0);
    f = dpm.getCell(q0-1,t0-1)->score;
    
    if (f + r > threshold) {       // check the Waterman condition
      if ((int)as.size() == k) as.push_back (curr);
      as[k].prepend (q0,t0);
      as[k].score = r;
      branch (dpm,as,q0-1,t0-1,k,threshold);
      k=(int)as.size();
    }
    
    // Deletions in the query w.r.t. the template
    
    for (int i = t0-2; i>0; --i) {
      f = dpm.getCell(q0-1,i)->score;
      g = e->deletion (*qseq,*tseq,q0-1,q0,i,t0);
      
      if (f + r - g > threshold) {   // check the Waterman condition
	if ((int)as.size() == k) as.push_back (curr);
	as[k].prepend (q0,t0);
	as[k].score = r - g;
	branch (dpm,as,q0-1,i,k,threshold);
	k=(int)as.size();
      }
    }

    // Insertions in the query w.r.t. the template
    
    for (int j = q0-2; j>0; --j) {
      f = dpm.getCell(j,t0-1)->score;
      g = e->insertion (*qseq,*tseq,j,q0,t0-1,t0);
      
      if (f + r - g > threshold) {   // check the Waterman condition
	if ((int)as.size() == k) as.push_back (curr);
	as[k].prepend (q0,t0);
	as[k].score = r - g;
	branch (dpm,as,j,t0-1,k,threshold);
	k=(int)as.size();
      }
    }
    
    if (k == k0) {
      
      // The score fell below threshold after extending current branch
      // It was good once, so just take optimal path to beginning.
      
      opt_path (dpm,as,q0,t0,k0,threshold);
      return;
    }
  }
}

template <class S1, class S2, class Etype>
void UnconstrainedNearOptimal<S1,S2,Etype>::
opt_path (const DPMatrix<S1,S2,Etype>& dpm,
	  AlignmentSet<S1,S2,Etype>& as,
	  int q0, int t0, int k0, float threshold)
{

  // Now we are not allowed to branch.  Force the optimal path until we 
  // reach the base case.
  
  const S1* qseq = dpm.getQuerySequence();
  const S2* tseq = dpm.getTemplateSequence();
  const Evaluator<S1,S2,Etype>* e = dpm.getEvaluator();
  
  int pq=-1,pt=-1;
  
  while ( t0>1 && q0>1 ) {

    // cerr << t0 << ":" << q0 << endl;
    
    as[k0].prepend (q0,t0);                  // Extend the alignment
    as[k0].score += dpm.getSim (q0,t0);      // Add the next match
    
    pq = dpm.getCell(q0,t0)->prev_query_idx;
    pt = dpm.getCell(q0,t0)->prev_template_idx;
      
    float g;                                 // Subtract g.p., if any
    if (q0-pq==1) g = e->deletion  (*qseq,*tseq,pq,q0,pt,t0);
    else          g = e->insertion (*qseq,*tseq,pq,q0,pt,t0);
    
    as[k0].score -= g;
    
    t0=pt;                                   // Update indices
    q0=pq;
    
  } 

  // Do base case
  
  as[k0].prepend (q0,t0);
  as[k0].prepend (0,0);
  as[k0].score += dpm.getCell(q0,t0)->score;

}

#endif  //_HMAP2_CW
