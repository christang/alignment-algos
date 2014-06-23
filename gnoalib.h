/**
 *  Package HMAP2.1
 *  File: gnoalib.h
 *  Desc: Define the gnoali evaluator scoring function
 *
 *  Created on 11/10/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#ifndef _HMAP2_GNOALIB
#define _HMAP2_GNOALIB

#include "hmapalib.h"
#include "gn2lib_seq.h"
#include "struct.h"

class GnoaliParams : 
  public HMAPaliParams 
{
  
public:
  
  GnoaliParams ();

  void read (ParamStore* p);

  static const float default_di_par1;
  static const float default_di_par2;
  static const float default_di_par3;
  static const float default_hb_par1;
  static const float default_hb_par2;
  static const float default_hb_par3;
  static const float default_ac_par1;
  static const float default_ac_par2;
  static const float default_ac_par3;
  static const float default_igo_alpha;
  static const float default_igo_beta;

  float di_par1;
  float di_par2;
  float di_par3;
  float hb_par1;
  float hb_par2;
  float hb_par3;
  float ac_par1;
  float ac_par2;
  float ac_par3;
  float igo_alpha;
  float igo_beta;

};


#define  dist_0       di_par1
#define  dist_off     di_par2
#define  dist_scale   di_par3
#define  hb_0         hb_par1
#define  hb_off       hb_par2
#define  hb_scale     hb_par3
#define  acc_0        ac_par1
#define  acc_off      ac_par2
#define  acc_scale    ac_par3
#define  qgo_alpha    igo_alpha
#define  qgo_beta     igo_beta

class GnoaliEval :
  public Evaluator<HMAPSequence,SMAPSequence,GnoaliEval> 
{

public:

  GnoaliEval (GnoaliParams& p);

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
    int i1 = t_pos1;
    int i2 = t_pos2;
    int di = i2 - i1;
    if (di < 2) return 0;

    int p1 = i1;     // i1,i2 = indices to nodes
    int p2 = i2 - 2; // p1,p2 = indices to matrices
    
    // broken-hbonds score:
    float br = float(t.brokenhb[p2][p1])/float(di-1);
    float b0 = br+params->hb_0; float bp = b0*b0/params->hb_scale;
    
    // distance-based score:
    float rd1 = t.distance[p2][p1];
    float rd2 = max(t.distance2[p2][p1]-7.0f,0.0f) - max(rd1-7.0f,0.0f);
    float ra = 0.735759f;
    float rd = 0;
    int sd = abs(t[t_pos2]->rdata.isse-t[t_pos1]->rdata.isse);
    if (sd>1) ra = exp(t.angle[p2][p1]) * 2;
    else      rd = exp(2*rd2/params->dist_scale);
    
    float gp = exp ((rd1+params->dist_0)/params->dist_scale) * ra + rd;
    
    // sse-based offset:
    float ro = (t[t_pos1]->rdata.isse>=0 && t[t_pos1]->rdata.isse==t[t_pos2]->rdata.isse)?params->dist_off:0;
    
    //cout << "(p2=" << p2 << ",p1=" << p1 << ",broken=" << t.brokenhb[p2][p1] << ",rd1=" << rd1 << ",rd2=" << rd2 << "," << bp << ",ra=" << ra << ",gp=" << gp << ",ro=" << ro << ") = " << (params->hb_off + bp)+(ro + gp) << endl;
    //       if (t1->data->residue->isInSameSSE(*t2->data->residue)) {
    //      cout << ".";
    //       } else {
    //      cout << "/";
    //       }
    
    switch (params->align_type) {
    case global:                   // overhangs penalized
    case global_local:             // overhangs penalized in template
      return (params->hb_off + bp)+(ro + gp);
    case local:                    // overhangs not penalized
    case semi_local:               // overhangs not penalized
    case local_global:             // overhangs penalized in query
      if (t[t_pos1]->isHead() || t[t_pos2]->isTail())
	return 0;
      else
	return (params->hb_off + bp)+(ro + gp);
    default:
      throw string ("Invalid align_type");
    }
  }
  
  inline float insertion  (const HMAPSequence& q,
			   const SMAPSequence& t,
			   int q_pos1, int q_pos2, 
			   int t_pos1, int t_pos2) const
  {
    int dist = q_pos2-q_pos1;
    if (dist < 2) return 0;
    
    // accessibility-based insertion penalty
    double a1 = t[t_pos1]->rdata.accessibility;
    double a2 = t[t_pos2]->rdata.accessibility;
    float ga = exp ((params->acc_0+(a1+a2)/2.0f)/params->acc_scale);
    
    // sse-based offset
    float ao = 0;
    if (t[t_pos1]->rdata.isse>=0 && t[t_pos1]->rdata.isse==t[t_pos2]->rdata.isse) {
      if (t[t_pos1]->rdata.sse_type==TC_Helix) ao = params->qgo_alpha;
      else ao = params->qgo_beta;
    } else {
    }
    
    switch (params->align_type) {
    case global:                   // overhangs penalized
    case local_global:             // overhangs penalized in query
      return ao + ga*(dist-1);
    case local:                    // overhangs not penalized
    case semi_local:               // overhangs not penalized
    case global_local:             // overhangs penalized in template
      if (q[q_pos1]->isHead() || q[q_pos2]->isTail())
	return 0;
      else
	return ao + ga*(dist-1);
    default:
      throw string ("Invalid align_type");
    }
  }

  void pre_calculate (const HMAPSequence& s1, const SMAPSequence& s2) const;
  void post_process (SimilarityMatrix& s) const;

private:
  
  GnoaliParams* params;

};

#endif  // _HMAP2_GNOALIB
