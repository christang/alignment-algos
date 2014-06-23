/**
 *  Package HMAP2.1
 *  File: hmap2eval.h
 *  Desc: Define the hmap scoring function
 *        (btwn HMAP query and SMAP template)
 *
 *  Deposited to cvs on 11/19/07
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_HMAP2EVAL
#define _HMAP2_HMAP2EVAL

#include "gn2_eval.h"

class Hmap2Eval :
  public Evaluator<HMAPSequence,SMAPSequence,Hmap2Eval> 
{

public:

  Hmap2Eval (Gn2Params& p);

  inline float similarity (const HMAPSequence& q,
			   const SMAPSequence& t,
			   int q_pos, int t_pos) const 
  {
    float ip = dot_product  (q[q_pos]->aa_profile,t[t_pos]->aa_profile);
    float pc = pearson_corr (q[q_pos]->sse_values,t[t_pos]->sse_values);
    
    float sim =
      ip * exp (params->alpha * pc *
		q[q_pos]->sse_confid *
		t[t_pos]->sse_confid);
    return sim;
  }
  
  inline float deletion   (const HMAPSequence& q,
			   const SMAPSequence& t,
			   int q_pos1, int q_pos2, 
			   int t_pos1, int t_pos2) const
  {
    int dist = t_pos2 - t_pos1;
    if (dist < 2) return 0;

    float gi, ge;
    gi = min(t[t_pos1]->gap_init(),t[t_pos2]->gap_init());
    ge = min(t[t_pos1]->gap_extn(),t[t_pos2]->gap_extn());

    switch (params->align_type) {
    case global:                   // overhangs penalized
    case global_local:             // overhangs penalized in template
      return gi + ge*(dist-2);
    case local:                    // overhangs not penalized
    case semi_local:               // overhangs not penalized
    case local_global:             // overhangs penalized in query
      if (t[t_pos1]->isHead() || t[t_pos2]->isTail())
	return 0;
      else
	return gi + ge*(dist-2);
    default:
      throw string ("Illegal gap style");
    }
  }
  
  inline float insertion  (const HMAPSequence& q,
			   const SMAPSequence& t,
			   int q_pos1, int q_pos2, 
			   int t_pos1, int t_pos2) const
  {
    int dist = q_pos2 - q_pos1;
    if (dist < 2) return 0;

    float gi, ge;
    gi = min(t[t_pos1]->gap_init(),t[t_pos2]->gap_init());
    ge = min(t[t_pos1]->gap_extn(),t[t_pos2]->gap_extn());

    switch (params->align_type) {
    case global:                   // overhangs penalized
    case local_global:             // overhangs penalized in query
      return gi + ge*(dist-2);
    case local:                    // overhangs not penalized
    case semi_local:               // overhangs not penalized
    case global_local:             // overhangs penalized in template
      if (q[q_pos1]->isHead() || q[q_pos2]->isTail())
	return 0;
      else
	return gi + ge*(dist-2);
    default:
      throw string ("Illegal gap style");
    }
  }

  void pre_calculate (const HMAPSequence& s1, const SMAPSequence& s2) const;
  inline void post_process (SimilarityMatrix& s) const {  
    norm_elements (s,s,1,s.rows()-1,1,s.cols()-1);
    shift_elements (s,s,1,s.rows()-1,1,s.cols()-1,-params->zero_shift);
  } ;

private:
  
  Gn2Params* params;
  mutable vector <float> v_gi,v_ge,v_cn;
  mutable vector <vector <float> > vv_gi,vv_ge,vv_cd;

};

#endif  // _HMAP2_HMAP2EVAL
