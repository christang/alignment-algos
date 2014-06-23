/**
 *  Package HMAP2.1
 *  File: kscw.h
 *  Desc: K-sorted constrained near-optimal alignment enumeration
 *
 *  Created on 11/24/07
 *  Author: cltang @ honig lab
 *
 *  copyright 2007.  all rights reserved.
 *
 */

#ifndef _HMAP2_KSCW
#define _HMAP2_KSCW

#include <string>

#include "alignment.h"
#include "enumerator.h"
#include "noalib.h"
#include "sflags.h"

using namespace std;

template <class S1, class S2, class Etype>
class KSConstrainedNearOptimal : public Enumerator<S1,S2,Etype> {
public:
  
  typedef AlignedPairList<S1,S2> SingleAlignment;
  typedef AlignedPair<S1,S2>     SinglePair;

  KSConstrainedNearOptimal (const NOaliParams& p,
			    const SuboptFlags& f) ;

  int estimateSize () const;

  void enumerate (DPMatrix<S1,S2,Etype>& dpm,
		  AlignmentSet<S1,S2,Etype>& as);

  struct op_data {
    unsigned int limit;
    int q0, t0, k0;
    float score, thresh, new_r;
    op_data (unsigned int l,int q,int t,int k,float th,float s=0.f,float n=0.f) :
      limit(l), q0(q), t0(t), k0(k), score(s), thresh(th), new_r(n) {};
    inline bool operator< (const op_data& a) const 
      { return score > a.score; } ;
  };

private:  

  void branch (const DPMatrix<S1,S2,Etype>& dpm,
	       AlignmentSet<S1,S2,Etype>& as,
	       op_data& data);

  void opt_path (const DPMatrix<S1,S2,Etype>& dpm,
		 AlignmentSet<S1,S2,Etype>& as,
		 op_data& data, bool force_opt=false);

  int numBranch () const;

  const NOaliParams* params;
  const SuboptFlags* subopt;
  bool warn_user;

};

template <class S1, class S2, class Etype>
KSConstrainedNearOptimal<S1,S2,Etype>::
KSConstrainedNearOptimal (const NOaliParams& p,
			  const SuboptFlags& f) 
  : params(&p),subopt(&f),warn_user(true) {}

template <class S1, class S2, class Etype>
int KSConstrainedNearOptimal<S1,S2,Etype>::estimateSize() const
{
  return params->number_suboptimal;
}

template <class S1, class S2, class Etype>
int KSConstrainedNearOptimal<S1,S2,Etype>::numBranch() const
{
  int c=1;
  bool f=(*subopt)[0];
  for (unsigned int i=1; i<subopt->size(); ++i) {
    if ((*subopt)[i]!=f) { ++c; f=!f; }
  }
  return c;
}

// Factorial

int factorial (int f) 
{ int r=1; for (int i=2;i<=f;++i) r*=i; return r; }

// For debugging stuff

ostream& operator<< (ostream& os, KSConstrainedNearOptimal<HMAPSequence,SMAPSequence,Gn2Eval>::op_data& op)
{ os<<"limit="<<op.limit<<",q0="<<op.q0<<",t0="<<op.t0<<",k0="<<op.k0
    <<",s="<<op.score<<",ns="<<op.new_r<<",t="<<op.thresh; return os; } ;

ostream& operator<< (ostream& os, KSConstrainedNearOptimal<HMAPSequence,SMAPSequence,Hmap2Eval>::op_data& op)
{ os<<"limit="<<op.limit<<",q0="<<op.q0<<",t0="<<op.t0<<",k0="<<op.k0
    <<",s="<<op.score<<",ns="<<op.new_r<<",t="<<op.thresh; return os; } ;

// Perform constrained Waterman-style branching tracebacks

template <class S1, class S2, class Etype>
void KSConstrainedNearOptimal<S1,S2,Etype>::
enumerate (DPMatrix<S1,S2,Etype>& dpm,
	   AlignmentSet<S1,S2,Etype>& as)
{
  int q_last = dpm.getQuerySize() - 1;
  int t_last = dpm.getTemplateSize() - 1;
  
  warn_user = true;

  as.reserve(min(estimateSize()*2,int(params->user_limit)));

  as.push_back(SingleAlignment());
  as.back().uid = 1;
  int k_last = as.size() - 1;

  float threshold = 
    (1.f - params->delta_ratio) * dpm.getCell(q_last,t_last)->score;
  threshold = min(threshold,dpm.getCell(q_last,t_last)->score-0.1f);
  op_data op (params->k_limit,q_last,t_last,k_last,threshold);
  branch (dpm, as, op);

  cerr<<"Ali#="<<as.size()<<endl;
  as.sortSet (params->number_suboptimal);
}

template <class S1, class S2, class Etype>
void KSConstrainedNearOptimal<S1,S2,Etype>::
branch (const DPMatrix<S1,S2,Etype>& dpm,
	AlignmentSet<S1,S2,Etype>& as,
	op_data& op)
{
  unsigned int k_limit(op.limit);
  int q0(op.q0), t0(op.t0), k0(op.k0);
  float threshold(op.thresh);

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

    float f,r,g,sum;
    SingleAlignment curr (as[k0]);
    
    // Operations are stored in k_sort and sorted when it is
    // filled.  Only the highest scoring branches up to k_limit
    // in number will be allowed to extend their alignments.

    vector<op_data> k_sort;
    k_sort.reserve (q0+t0);
    
    // First, check if user limits on number of alignments was exceeded.
    // Force the optimal alignment to the end if yes.

    if (as.size() > params->user_limit) {
      if (warn_user) {
	warn_user = false;
	cerr << "Number of alignments exceeding user limits (";
	cerr << params->user_limit << ")" << endl;
	cerr << "User is advised to lower current threshold (";
	cerr << params->delta_ratio  << ") or further constrain alignments";
	cerr << endl << endl;
      }
      opt_path (dpm,as,op,true);
      return;
    }
    
    // Check if more space is needed.

    if (as.size() > 0.9 * as.capacity() && as.capacity() < params->user_limit)
     	as.reserve (min(as.capacity()*2,params->user_limit));
    
    // Begin the recursive step.

    const S1* qseq = dpm.getQuerySequence();
    const S2* tseq = dpm.getTemplateSequence();
    const Evaluator<S1,S2,Etype>* e = dpm.getEvaluator();
    
    // Match case.  f and r are used to store the forward and reverse
    // traceback scores.  Their sum must exceed threshold to traverse
    // the branch.  This is Waterman's condition for near-optimal 
    // alignments.

    // Save each operation in k_sort, with limit=k_limit/2.
    // (Only the best alignment branch will keep the original k_limit.
    // See below how we make this work.)

    r = curr.score + dpm.getSim (q0,t0);
    f = dpm.getCell(q0-1,t0-1)->score;
    sum = f + r;

    if (sum > threshold)        // check the Waterman condition
      k_sort.push_back (op_data (k_limit/2,q0-1,t0-1,k0,threshold,sum,r));
    
    // Deletions in the query w.r.t. the template
    
    for (int i = t0-2; i>0; --i) {
      f = dpm.getCell(q0-1,i)->score;
      g = e->deletion (*qseq,*tseq,q0-1,q0,i,t0);
      sum = f + r - g;
      
      if (sum > threshold)     // check the Waterman condition
	k_sort.push_back (op_data (k_limit/2,q0-1,i,k0,threshold,sum,r-g));
    }
    
    // Insertions in the query w.r.t. the template
    
    for (int j = q0-2; j>0; --j) {
      f = dpm.getCell(j,t0-1)->score;
      g = e->insertion (*qseq,*tseq,j,q0,t0-1,t0);
      sum = f + r - g;
      
      if (sum > threshold)     // check the Waterman condition
	k_sort.push_back (op_data (k_limit/2,j,t0-1,k0,threshold,sum,r-g));
    }
    
    // Check for operations

    if (k_sort.size() == 0) {
      
      // The score fell below threshold after extending current branch
      // It was good once, so just take optimal path to beginning.
      
      op_data new_op (1,q0,t0,k0,threshold);
      opt_path (dpm,as,new_op,true);
      return;

    } else {

      // Sort and keep only the top k_limit operations.

      if (k_sort.size()>k_limit) {
	partial_sort (k_sort.begin(),k_sort.begin()+k_limit,k_sort.end());
	k_sort.erase (k_sort.begin()+k_limit,k_sort.end());
	// cerr << "sz=" << k_sort.size() << endl;
      } else {
	sort (k_sort.begin(),k_sort.end());
      }

      // Process the sorted operations.
      // (Extend good alignments, increment alignment indices, send
      // operations on to next step.)
      

      typename vector<op_data>::iterator it=k_sort.begin();
      it->limit*=2;  // Let limit=k_limit only for first op.

      int count=0;
      for (int k=k0;it!=k_sort.end();++it) {
	it->k0=k;
	if ((int)as.size() == k) { as.push_back (curr); as[k].uid=k; }
	as[k].prepend (q0,t0);
	as[k].score = it->new_r;
	cerr << "(" << count << ") uid = " << as[k0].uid << " -> " << as[k].uid << endl;
	cerr << *it << endl;
	opt_path (dpm,as,*it);
	k=(int)as.size();
	count ++;
      }

    }
  }
}

template <class S1, class S2, class Etype>
void KSConstrainedNearOptimal<S1,S2,Etype>::
opt_path (const DPMatrix<S1,S2,Etype>& dpm,
	  AlignmentSet<S1,S2,Etype>& as,
	  op_data& op, bool force_opt)
{
  unsigned int k_limit(op.limit);
  int q0(op.q0), t0(op.t0), k0(op.k0);
  float threshold(op.thresh);

  // Check if k_limit = 1.  If yes, force optimal alignment.
  if (k_limit<=1) force_opt = true;

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

    // #1 Branch point allowed by a change in subopt flag state
    bool flag=!(*subopt)[t0]; 

    while ( t0>1 && q0>1 ) {

      // as[k0].print_pairs(); cout << endl;
      // cerr << t0 << ":" << q0 << " - " << flag << (*subopt)[t0] << endl;
    
      // #1 Check for change in subopt flag state
      if ( !force_opt && (*subopt)[t0]==flag ) break; 

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

    op_data op(k_limit,pq,pt,k0,threshold);
    //cerr << "-o- " << op << endl;
    branch (dpm,as,op);
      
  }
}

#endif  //_HMAP2_KSCW
