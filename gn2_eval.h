/**
 *  Package HMAP2.1
 *  File: gn2eval.h
 *  Desc: Define the gn2 evaluator scoring function
 *
 *  Deposited to cvs on 11/19/07
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_GN2EVAL
#define _HMAP2_GN2EVAL

#include "hmap_eval.h"
#include "gn2lib_seq.h"
#include "struct.h"

class Gn2Params : 
  public HMAPaliParams 
{
  
public:
  
  Gn2Params ();

  void read (ParamStore* p);

  vector <float> ss_lods;

  static const float default_gn2_gi_coil;
  static const float default_gn2_ge_coil;
  static const float default_gn2_gi_ss;
  static const float default_gn2_ge_ss;
  static const float default_aa_weight;
  static const float default_ss_weight;
  static const float default_cn_weight;
  static const float default_hp_weight;
  static const float default_hb_weight;
  static const float default_ic_weight;
  static const float default_dd_constr;
  static const float default_gn2_shift;
  static const bool  default_ss_dependent_gp;  

  float gap_init_coil;
  float gap_extn_coil;
  float gap_init_ss;
  float gap_extn_ss;
  float aa_weight;
  float ss_weight;
  float cn_weight;
  float hp_weight;
  float hb_weight;
  float ic_weight;
  float dd_constr;
  float gn2_shift;
  bool  ss_dependent_gp;

};


class Gn2Eval :
  public Evaluator<HMAPSequence,SMAPSequence,Gn2Eval> 
{

public:

  Gn2Eval (Gn2Params& p);

  inline float similarity (const HMAPSequence& q,
			   const SMAPSequence& t,
			   int q_pos, int t_pos) const 
  {
    float ip = norm_dot_product (q[q_pos]->aa_profile,t[t_pos]->aa_profile);
    unsigned int lods_idx = t[t_pos]->lods_type * 12 + q[q_pos]->lods_type;
    
    float log_aa = 0.543f / ( 2.85f - exp (ip) ) - 0.738f; 
    //float log_aa = -0.416f + exp ( -7.262f + exp ( 1.1f + ip ) ); 

    float log_ss = params->ss_lods [ lods_idx ];
    float log_cn = 2.f * t.weighted_contact_number [ t_pos ] - 0.9f ;
    float log_hp = exp(exp(-abs(q[q_pos]->hydropathy-t[t_pos]->hydropathy))
		       *(0.75f+0.3f*abs(t[t_pos]->hydropathy-0.22f))) - 1.8f;

    //float downweight_loops = 0.2f;
    //if (t[t_pos]->p_coil() > 0.5) log_aa *= downweight_loops;

    float sim = 
      params->gn2_shift +
      params->aa_weight * log_aa +
      params->ss_weight * log_ss +
      params->cn_weight * log_cn +
      params->hp_weight * log_hp;
    
    return sim;
  }
  
  inline float deletion   (const HMAPSequence& q,
			   const SMAPSequence& t,
			   int q_pos1, int q_pos2, 
			   int t_pos1, int t_pos2) const
  {
    int di = t_pos2 - t_pos1;
    if (di < 2) return 0;

    int p1 = t_pos1;     // i1,i2 = indices to nodes
    int p2 = t_pos2 - 2; // p1,p2 = indices to matrices
    
    // GAP PENALTY:
    float GP = 8100.f;
    if (t.distance [p2][p1] < 18.f)
        // Affine score                          // Const. hb and dist
      GP = vv_gi[p2][p1] + vv_ge[p2][p1] * (di-2) + vv_cd[p2][p1];

    switch (params->align_type) {
    case global:                   // overhangs penalized
    case global_local:             // overhangs penalized in template
      return GP;
    case local:                    // overhangs not penalized
    case semi_local:               // overhangs not penalized
    case local_global:             // overhangs penalized in query
      if (t[t_pos1]->isHead() || t[t_pos2]->isTail())
	return 0;
      else
	return GP;
    default:
      throw string ("Invalid align_type");
    }
  }
  
  inline float insertion  (const HMAPSequence& q,
			   const SMAPSequence& t,
			   int q_pos1, int q_pos2, 
			   int t_pos1, int t_pos2) const
  {
    int di = q_pos2-q_pos1;
    if (di < 2) return 0;

    // GAP PENALTY:
            // Affine score                        // Penalty within high CN
    float GP = v_gi[t_pos1] + v_ge[t_pos1] * (di-2) + v_cn[t_pos1];

    switch (params->align_type) {
    case global:                   // overhangs penalized
    case local_global:             // overhangs penalized in query
      return GP;
    case local:                    // overhangs not penalized
    case semi_local:               // overhangs not penalized
    case global_local:             // overhangs penalized in template
      if (q[q_pos1]->isHead() || q[q_pos2]->isTail())
	return 0;
      else
	return GP;
    default:
      throw string ("Invalid align_type");
    }
  }

  void pre_calculate (const HMAPSequence& s1, const SMAPSequence& s2) const;
  inline void post_process (SimilarityMatrix& s) const {} ;

private:
  
  Gn2Params* params;
  mutable vector <float> v_gi,v_ge,v_cn;
  mutable vector <vector <float> > vv_gi,vv_ge,vv_cd;

};

#endif  // _HMAP2_GN2EVAL
