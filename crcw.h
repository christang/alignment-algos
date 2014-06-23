/**
 *  Package HMAP2.1
 *  File: crcw.h
 *  Desc: Controlled redundancy near-optimal alignment enumeration
 *
 *  Created on 11/24/07
 *  Author: cltang @ honig lab
 *
 *  copyright 2007.  all rights reserved.
 *
 */

#ifndef _HMAP2_CRCW
#define _HMAP2_CRCW

#include <string>

#include <iomanip>

#include "alignment.h"
#include "enumerator.h"
#include "noalib.h"
#include "sflags.h"

using namespace std;

template <class S1, class S2, class Etype>
class CRConstrainedNearOptimal : public Enumerator<S1,S2,Etype> {
public:
  
  typedef AlignedPairList<S1,S2> SingleAlignment;
  typedef AlignedPair<S1,S2>     SinglePair;

  CRConstrainedNearOptimal (const NOaliParams& p,
			    const SuboptFlags& f) ;

  // width of tree = p->k_limit
  // limit on maximum operations sorted = p->sort_limit
  // limit on maximum alignments generated = p->user_limit
  float threshold;

  int estimateSize () const;

  void enumerate (DPMatrix<S1,S2,Etype>& dpm,
		  AlignmentSet<S1,S2,Etype>& as);

  struct op_data {  // operation
    unsigned int limit, index;
    int q0, t0, k0;
    float score, new_r;
    op_data (unsigned int l,int q,int t,int k,float s=0.f,float n=0.f) :
      limit(l), index(0), q0(q), t0(t), k0(k), score(s), new_r(n) {};
    inline bool operator< (const op_data& a) const 
      { return score > a.score; } ;
  };

private:

  void branch (const DPMatrix<S1,S2,Etype>& dpm,
	       AlignmentSet<S1,S2,Etype>& as,
	       op_data& data);

  void filter_and_extend (const DPMatrix<S1,S2,Etype>& dpm,
			  AlignmentSet<S1,S2,Etype>& as,
			  int q0, int t0,
			  vector<op_data>& v_op);

  void force_opt_path (const DPMatrix<S1,S2,Etype>& dpm,
		       AlignmentSet<S1,S2,Etype>& as,
		       op_data& data);

  void   init_mem (unsigned int t_last);
  void reinit_mem (unsigned int t_last, unsigned int o_last);
  void   free_mem ();

  int numBranch () const;

  const NOaliParams* params;
  const SuboptFlags* subopt;
  bool warn_user;

  int **alignments;
  int *regions;

  unsigned int count_redundant;
  unsigned int count_subpaths;

};

template <class S1, class S2, class Etype>
CRConstrainedNearOptimal<S1,S2,Etype>::
CRConstrainedNearOptimal (const NOaliParams& p,
			  const SuboptFlags& f) 
  : params(&p),subopt(&f),warn_user(true) {}

template <class S1, class S2, class Etype>
int CRConstrainedNearOptimal<S1,S2,Etype>::estimateSize() const
{
  return params->number_suboptimal;
}

template <class S1, class S2, class Etype>
int CRConstrainedNearOptimal<S1,S2,Etype>::numBranch() const
{
  int c=1;
  bool f=(*subopt)[0];
  for (unsigned int i=1; i<subopt->size(); ++i) {
    if ((*subopt)[i]!=f) { ++c; f=!f; }
  }
  return c;
}

// For debugging stuff

class HMAPSequence;
class SMAPSequence;
class Gn2Eval;
class Hmap2Eval;

ostream& operator<< (ostream& os, CRConstrainedNearOptimal<HMAPSequence,SMAPSequence,Gn2Eval>::op_data& op)
{ os<<"limit="<<op.limit<<",q0="<<op.q0<<",t0="<<op.t0<<",k0="<<op.k0
    <<",s="<<op.score<<",ns="<<op.new_r; return os; } ;

ostream& operator<< (ostream& os, CRConstrainedNearOptimal<HMAPSequence,SMAPSequence,Hmap2Eval>::op_data& op)
{ os<<"limit="<<op.limit<<",q0="<<op.q0<<",t0="<<op.t0<<",k0="<<op.k0
    <<",s="<<op.score<<",ns="<<op.new_r; return os; } ;

//void pause () 
//{ string x; cout<< "pause"; cin >> x; cout<< endl; }

// Perform constrained Waterman-style branching tracebacks

template <class S1, class S2, class Etype>
void CRConstrainedNearOptimal<S1,S2,Etype>::
enumerate (DPMatrix<S1,S2,Etype>& dpm,
	   AlignmentSet<S1,S2,Etype>& as)
{
  int q_last = dpm.getQuerySize() - 1;
  int t_last = dpm.getTemplateSize() - 1;

  init_mem (t_last);
 
  warn_user = true;
  as.reserve(min(estimateSize()*2,int(params->user_limit)));

  as.push_back(SingleAlignment());
  as.back().uid=1;

  int init = as.size() - 1;

  threshold =(1.f - params->delta_ratio) * dpm.getCell(q_last,t_last)->score;
  threshold = min(threshold,dpm.getCell(q_last,t_last)->score-0.1f);

  op_data op (params->k_limit,q_last,t_last,init);
  count_redundant = 0;
  count_subpaths = 0;
  branch (dpm, as, op);

  cerr << "Removed " << count_redundant << " subpaths with more than "
       << params->max_overlap*100 << "% overlap. Started with " << count_subpaths 
       << "." << endl;

  cerr<<"Number of alignments before sorting: "<<as.size()<<"."<<endl;
  as.sortSet (params->number_suboptimal);

  free_mem ();
}

template <class S1, class S2, class Etype>
void CRConstrainedNearOptimal<S1,S2,Etype>::init_mem (unsigned int t_last)
{
  try {
    alignments = new int* [params->sort_limit];
    for (unsigned int i=0; i<params->sort_limit; ++i)
      alignments[i] = new int [t_last];
    
    int state = 0;
    regions = new int [t_last];
    for (unsigned int i=0; i<subopt->size()-1; ++i) {
      if ((*subopt)[i+1]!=(*subopt)[i]) ++state;
      regions[i] = state;
    }
  } catch (bad_alloc) {
    cerr << "Memory fault!" << endl;
    exit(-1);
  }
}

template <class S1, class S2, class Etype>
void CRConstrainedNearOptimal<S1,S2,Etype>::reinit_mem (unsigned int t_last,
							unsigned int o_last)
{
  for (unsigned int i=0; i<o_last; ++i)
    for (unsigned int j=0; j<t_last; ++j)
      alignments[i][j] = -1;
}

template <class S1, class S2, class Etype>
void CRConstrainedNearOptimal<S1,S2,Etype>::free_mem ()
{
  delete[] regions;
  for (unsigned int i=0; i<params->sort_limit; ++i) delete[] alignments[i];
  delete[] alignments;
}

template <class S1, class S2, class Etype>
void CRConstrainedNearOptimal<S1,S2,Etype>::
branch (const DPMatrix<S1,S2,Etype>& dpm,
	AlignmentSet<S1,S2,Etype>& as,
	op_data& op)
{
  unsigned int k_limit(op.limit);
  int q0(op.q0), t0(op.t0), k0(op.k0);

  // There are two options here.  We are allowed to branch or 
  // follow the optimal path.  Branching is controlled by the
  // subopt array and k_limit.
  
  // cerr << "Beg: " << op << endl;
      
  if ( k_limit < 2 ) { force_opt_path (dpm,as,op); return; }

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
    cerr << "user_limit: " << op << endl;
    force_opt_path (dpm,as,op); return;
  }
  
  // Check if more space is needed.
  
  if (as.size() > 0.9 * as.capacity() && as.capacity() < params->user_limit)
    as.reserve (min(as.capacity()*2,params->user_limit));
  
  // Operations are stored in all_op and sorted when it is
  // filled.  Only the highest scoring branches up to k_limit
  // in number will be allowed to extend their alignments.
  
  vector<op_data> all_op;
  all_op.reserve (q0+t0);
  
  // Begin the recursive step.
  
  const S1* qseq = dpm.getQuerySequence();
  const S2* tseq = dpm.getTemplateSequence();
  const Evaluator<S1,S2,Etype>* e = dpm.getEvaluator();
  
  // Match case.  f and r are used to store the forward and reverse
  // traceback scores.  Their sum must exceed threshold to traverse
  // the branch.  This is Waterman's condition for near-optimal 
  // alignments.
  
  float f,r,g,sum;
  SingleAlignment curr (as[k0]);
  
  // Save each operation in all_op, with limit=k_limit-1.
  // (Only the best alignment branch will keep the original k_limit.
  // See below how we make this work.)
  
  r = curr.score + dpm.getSim (q0,t0);
  f = dpm.getCell(q0-1,t0-1)->score;
  sum = f + r;
  
  if (sum > threshold)        // check the Waterman condition
    all_op.push_back (op_data (k_limit,q0-1,t0-1,k0,sum,r));
  
  // Deletions in the query w.r.t. the template
  
  for (int i = t0-2; i>0; --i) {
    f = dpm.getCell(q0-1,i)->score;
    g = e->deletion (*qseq,*tseq,q0-1,q0,i,t0);
    sum = f + r - g;
    
    if (sum > threshold)     // check the Waterman condition
      all_op.push_back (op_data (k_limit,q0-1,i,k0,sum,r-g));
  }
  
  // Insertions in the query w.r.t. the template
  
  for (int j = q0-2; j>0; --j) {

    f = dpm.getCell(j,t0-1)->score;
    g = e->insertion (*qseq,*tseq,j,q0,t0-1,t0);
    sum = f + r - g;
    
    if (sum > threshold)     // check the Waterman condition
      all_op.push_back (op_data (k_limit,j,t0-1,k0,sum,r-g));
    }
  
  // Check for operations
  
  if (all_op.size() == 0) {
    
    // The score fell below threshold after extending current branch
    // It was good once, so just take optimal path to beginning.
    
    force_opt_path (dpm,as,op);
    return;
    
  } else {
    
    // Sort and process the operations.

    if (all_op.size()>params->sort_limit) {
      partial_sort (all_op.begin(),all_op.begin()+params->sort_limit,all_op.end());
      all_op.erase (all_op.begin()+params->sort_limit,all_op.end());
    } else {
      sort (all_op.begin(),all_op.end());
    }
    
    filter_and_extend (dpm,as,q0,t0,all_op);
    
    // Continue recurrance for alignments that pass the 
    // filter criteria (where all_op.k0 != -1).
    
    typename vector<op_data>::iterator it=all_op.begin();

    int count = 0;
    for (;it!=all_op.end();++it) {
      // if ( it->k0 > -1 ) cerr << "(" << count << ") uid = " << as[k0].uid << " -> " << as[it->k0].uid << endl;
      // cerr << *it << endl;
      if ( it->k0 > -1 ) branch (dpm,as,*it);
      // cerr << "Fin: " << *it << endl;
      ++count ;
    }
    
  }
}

// Filters incoming v_op operations for alignments that 
// are non-redundant and maximal in score, then fills the 
// AlignmentSet with these alignments.  k0 is set to -1 for
// alignment operations that have finished.

template <class S1, class S2, class Etype>
void CRConstrainedNearOptimal<S1,S2,Etype>::
filter_and_extend (const DPMatrix<S1,S2,Etype>& dpm,
		   AlignmentSet<S1,S2,Etype>& as,
		   int q0, int t0,
		   vector<op_data>& v_op)
{

  // Define the completion length for finished alignments.

  int end_alignment=2;

  // Allocate memory.

  bool *filter = new bool [v_op.size()];  // redundancy filter, true=keep
  int  *p_rq   = new  int [v_op.size()];  // subpath start query index
  int  *p_rt   = new  int [v_op.size()];  // subpath start templ index
  int  *l_sp   = new  int [v_op.size()];  // length of subpath
  int  *state  = new  int [v_op.size()];  // end state(region) of subpath
  float *rs   = new float [v_op.size()];  // accumulating reverse score

  // Fill up alignment matrix

  reinit_mem (t0,v_op.size());

  count_subpaths += v_op.size();

  const S1* qseq = dpm.getQuerySequence();
  const S2* tseq = dpm.getTemplateSequence();
  const Evaluator<S1,S2,Etype>* e = dpm.getEvaluator();

  // cerr << "F&E: " << qseq->olc(q0) << q0 << tseq->olc(t0) << t0 << endl;

  for (unsigned int i=0; i<v_op.size(); ++i) {

    float g;

    int pq,pt;
    
    v_op[i].index = i;
    int q = v_op[i].q0;
    int t = v_op[i].t0;
    l_sp[i] = 1;
    state[i] = regions[t-1];   // Starting region
    rs[i] = v_op[i].new_r;

    while (q>0 && t>0 && regions[t-1] == state[i]) {

      // Store alignment and score

      alignments[i][t-1]=q;
      ++l_sp[i];

      pq = dpm.getCell(q,t)->prev_query_idx;
      pt = dpm.getCell(q,t)->prev_template_idx;
      if (q-pq==1) g = e->deletion  (*qseq,*tseq,pq,q,pt,t);
      else         g = e->insertion (*qseq,*tseq,pq,q,pt,t);
      rs[i] += dpm.getSim (q,t);
      rs[i] -= g;

      // Increment alignment

      q = pq;
      t = pt;
      
    }

    p_rq [i] = q;   // Base case when pq and pt are 0.
    p_rt [i] = t;
    state[i] = regions[t-1];   // Set ending region

  }

//   for (unsigned int i=0; i<v_op.size(); ++i)
//     { cerr << "[" << i << "] (" << p_rq [i] << "," << p_rt [i] << ") ";
//     for (int j=0; j<t0-1; ++j) cerr << setw(3) << alignments[i][j] << " ";
//     cerr << " (" << q0 << "," << t0 << ")" << endl; }

  // pause();

  // Set filter flags

  for (unsigned int i=1; i<v_op.size(); ++i) filter[i]=false;
  
  filter [0] = true;
  unsigned int count = 0;
  unsigned int accepted = 1;
  unsigned int lim = v_op.back().limit;

  for (unsigned int i=1; i<v_op.size() && accepted < lim; ++i) {
    filter [i] = true;
    for (unsigned int j=0; j<i; ++j) {
      if (filter [i] && filter [j] && state [i] == state [j]) {
	float overlap = 0.f;
	float overlap_max = params->max_overlap*(float)l_sp[j];
	//cerr << i << " with "<< j << " (" << p_rq[i] << "==" << p_rq[j] << "&&" <<  p_rt[i]<<"=="<< p_rt[j] << ") "; 
	if (p_rq[i]==p_rq[j] && p_rt[i]==p_rt[j]) { 
	  ++overlap; 
	  //cerr << "* "; 
	}
	for (int k=t0-1; k>=p_rt[i]; --k) {
	  if ((alignments[i][k]>-1 && alignments[j][k]>-1 &&
	       alignments[i][k]==alignments[j][k])) {
	    ++overlap; //cerr << alignments[i][k] << "|" << k << " ";
	    if (overlap > overlap_max) {
	      //cerr << "Rej: " <<" ov=" << overlap << "|" << overlap_max << "/" << l_sp[j] << endl; 
	      filter [i] = false;
	      ++count;
	      break;
	    }
	  }
	}
	//cerr << endl;
      }
    }
    if (filter [i]) ++accepted;
    //cerr << "no. accepted " << accepted << endl;
  }

  count_redundant += count;

//    for (unsigned int i=0; i<v_op.size(); ++i) { 
//      if (filter[i]) {
//        cerr << "Accept [" << i << "] (" << p_rq[i] << "," << p_rt[i]<<") ";
//        for (int j=0; j<t0-1; ++j) cerr << setw(3) << alignments[i][j]<<" ";
//        cerr << " (" << q0 << "," << t0 << ")" << endl; } }
   
  // Filter operations

  vector<op_data> tmp_op;
  tmp_op.reserve(count);

  accepted = 0;
  for (unsigned int i=0; i<v_op.size() && accepted<lim; ++i)
    if (filter[i]) {
      tmp_op.push_back(v_op[i]);
      ++accepted;
    }
  tmp_op.swap(v_op);

  //cerr << "size v_op=" << v_op.size() << endl;

  for (unsigned int i=1; i<v_op.size(); ++i) 
    v_op[i].limit=max(2u,lim/2);  // Let limit=k_limit only for optimal (e.g.
		                  // where i>=1). All other limit=k_limit/2.

  // cerr << "v_op=" << v_op.size() << endl;

  // Fill alignment set using filters

  int k = v_op[0].k0;
  SingleAlignment curr (as[k]);

  for (unsigned int i=0; i<v_op.size(); ++i) {

    int q0_i = v_op[i].index;

    // Add this alignment to the set.
    
    if (k == (int) as.size()) { as.push_back (curr); as[k].uid=k; }
    
    as[k].prepend (q0,t0);
    for (int j=t0-1; j>p_rt[q0_i]; --j) {
      int ali_q0 = alignments[q0_i][j-1];
      if (ali_q0 > -1) as[k].prepend (ali_q0,j);
    }
    as[k].score = rs[q0_i];
    
    // Set new q0,t0,k0 for next round
    
    v_op[i].q0=p_rq[q0_i];
    v_op[i].t0=p_rt[q0_i];
    v_op[i].k0=k;
    
    // Check for base condition and finalize if it is.

    // cerr << "Fill1: " << v_op[i] << endl;
    // as[k].print_pairs();
    
    if (p_rq[q0_i] <= end_alignment || p_rt[q0_i] <= end_alignment) {
      force_opt_path (dpm,as,v_op[i]);
      v_op[i].k0 = -1;
    }

    // Check for branch width = 1 and finalize if it is

    //else if (v_op.size() == 1) {
    //  force_opt_path (dpm,as,v_op[i]);
    //  v_op[i].k0 = -1;
    //}
    
    // cerr << "Fill2: " << v_op[i] << endl;
    // as[k].print_pairs();
    k = (int) as.size();
      
  }

  // Clean up

  delete[] rs;
  delete[] state;
  delete[] l_sp;
  delete[] p_rt;
  delete[] p_rq;
  delete[] filter;

}

template <class S1, class S2, class Etype>
void CRConstrainedNearOptimal<S1,S2,Etype>::
force_opt_path (const DPMatrix<S1,S2,Etype>& dpm,
		AlignmentSet<S1,S2,Etype>& as,
		op_data& op)
{

  // Force the optimal path until we have reached the base case.
  
  const S1* qseq = dpm.getQuerySequence();
  const S2* tseq = dpm.getTemplateSequence();
  const Evaluator<S1,S2,Etype>* e = dpm.getEvaluator();
  
  int pq=-1,pt=-1;
  int q0(op.q0), t0(op.t0), k0(op.k0);
  
  while ( t0>0 && q0>0 ) {
    
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

  as[k0].prepend (0,0);

  //cerr << "Opt: " << op << endl;
  //as[k0].print_pairs();

}

#endif  //_HMAP2_CRCW
