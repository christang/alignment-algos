/**
 *  Package HMAP2.1
 *  File: cw.h
 *  Desc: Constrained near-optimal alignment enumeration
 *
 *  Created on 11/9/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_CW
#define _HMAP2_CW

#include <algorithm>
#include <string>

#include "alignment.h"
#include "enumerator.h"
#include "noalib.h"
#include "sflags.h"

using namespace std;

template <class S1, class S2, class Etype>
class ConstrainedNearOptimal : public Enumerator<S1,S2,Etype> {
public:
  
  typedef AlignedPairList<S1,S2> SingleAlignment;
  typedef AlignedPair<S1,S2>     SinglePair;

  ConstrainedNearOptimal (const NOaliParams& p,
			  const SuboptFlags& f) ;

  unsigned int user_limit;

  inline int estimateSize () const { return params->number_suboptimal; } ;

  void enumerate (DPMatrix<S1,S2,Etype>& dpm,
		  AlignmentSet<S1,S2,Etype>& as);

private:

  void branch (const DPMatrix<S1,S2,Etype>& dpm,
	       AlignmentSet<S1,S2,Etype>& as,
	       int q0, int t0, int k0, float thresh, bool force_opt=false);

  void opt_path (const DPMatrix<S1,S2,Etype>& dpm,
		 AlignmentSet<S1,S2,Etype>& as,
		 int q0, int t0, int k0, float thresh, bool force_opt=false);

  const NOaliParams* params;
  const SuboptFlags* subopt;
  bool warn_user;

};

template <class S1, class S2, class Etype>
ConstrainedNearOptimal<S1,S2,Etype>::
ConstrainedNearOptimal (const NOaliParams& p,
			const SuboptFlags& f) 
  : params(&p), subopt(&f), warn_user(true) {}

// Perform constrained Waterman-style branching tracebacks

template <class S1, class S2, class Etype>
void ConstrainedNearOptimal<S1,S2,Etype>::
enumerate (DPMatrix<S1,S2,Etype>& dpm,
	   AlignmentSet<S1,S2,Etype>& as)
{
  int q_last = dpm.getQuerySize() - 1;
  int t_last = dpm.getTemplateSize() - 1;

  warn_user = true;
  user_limit = 1000000;  
                         // This should be a parameter, but I'll 
                         // shortcut this for now.

  as.reserve(min(estimateSize()*20,int(user_limit)));

  as.push_back(SingleAlignment());
  as.back().uid = 0;
  int k_last = as.size() - 1;

  float threshold = 
    (1.f - params->delta_ratio) * dpm.getCell(q_last,t_last)->score;
  threshold = min(threshold,dpm.getCell(q_last,t_last)->score-0.1f);
  branch (dpm, as, q_last, t_last, k_last, threshold);
  cerr<<"Ali#="<<as.size()<<endl;
  as.sortSet (params->number_suboptimal);
}

template <class S1, class S2, class Etype>
void ConstrainedNearOptimal<S1,S2,Etype>::
branch (const DPMatrix<S1,S2,Etype>& dpm,
	AlignmentSet<S1,S2,Etype>& as,
	int q0, int t0, int k0, float threshold, bool force_opt)
{
  if (q0==1 || t0==1) {

    // This is the base case.  If either q0 or t0 reaches 1, we
    // terminate the recursion and complete the k0'th alignment.

    as[k0].prepend (q0,t0);
    as[k0].prepend (0,0);
    as[k0].score += dpm.getCell(q0,t0)->score;
    
  } else {

    // There are two options here.  We are allowed to branch or 
    // follow the optimal path.  Branching is controlled by the
    // subopt array.

    if (!force_opt) {
      
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
	force_opt = true;
 	opt_path (dpm,as,q0,t0,k0,threshold,true);
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
	opt_path (dpm,as,q0-1,t0-1,k,threshold,force_opt);
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
	  opt_path (dpm,as,q0-1,i,k,threshold,force_opt);
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
	  opt_path (dpm,as,j,t0-1,k,threshold,force_opt);
	  k=(int)as.size();
	}
      }

      if (k == k0) {

      	// The score fell below threshold after extending current branch
	// It was good once, so just take optimal path to beginning.

	opt_path (dpm,as,q0,t0,k0,threshold,true);
	return;
      }

    } else {  
      
	opt_path (dpm,as,q0,t0,k0,threshold,force_opt);
	return;
    }
  
  }

}

template <class S1, class S2, class Etype>
void ConstrainedNearOptimal<S1,S2,Etype>::
opt_path (const DPMatrix<S1,S2,Etype>& dpm,
	  AlignmentSet<S1,S2,Etype>& as,
	  int q0, int t0, int k0, float threshold, bool force_opt)
{
  if (q0==1 || t0==1) {

    // This is the base case.  If either q0 or t0 reaches 1, we
    // terminate the recursion and complete the k0'th alignment.

    as[k0].prepend (q0,t0);
    as[k0].prepend (0,0);
    as[k0].score += dpm.getCell(q0,t0)->score;
    
  } else {

    // Now we are not allowed to branch.  Force the optimal path until we 
    // reach a branch point, or check to see if we have reached the base 
    // case.
    
    const S1* qseq = dpm.getQuerySequence();
    const S2* tseq = dpm.getTemplateSequence();
    const Evaluator<S1,S2,Etype>* e = dpm.getEvaluator();
    
    int pq=-1,pt=-1;

    // Below is a choice of two branch point control rules
    // Uncomment related lines to specify rules

    // #1 Branch point allowed by a change in subopt flag state
    bool flag=!(*subopt)[t0]; 

    // #2 Branch point allowed by  entering the next subopt=true region
    // bool flag=false;

    while ( t0>1 && q0>1 ) {

      // as[k0].print_pairs();
      // cerr << t0 << ":" << q0 << " - " << flag << (*subopt)[t0] << endl;
    
      // #1 Check for change in subopt flag state
      if ( !force_opt && (*subopt)[t0]==flag ) break; 

      // #2 Toggle flag with subopt=false; when flag=true, allow branch
      // if ( !(*subopt)[t0] ) flag=true;
      // if ( !force_opt && (*subopt)[t0] && flag ) break; 
      
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

    // Do base case or evaluate branching.
    
    branch (dpm,as,pq,pt,k0,threshold,force_opt);
      
  }
}

#endif  //_HMAP2_CW
