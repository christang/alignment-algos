/**
 *  Package HMAP2.1
 *  File: gn2lib.cpp
 *  Desc: Define the gn2 evaluator scoring function
 *
 *  Deposited to cvs on 11/19/07
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include "gn2_eval.h"

const float Gn2Params::default_gn2_gi_coil =  1.2;
const float Gn2Params::default_gn2_ge_coil =  0.08;
const float Gn2Params::default_gn2_gi_ss   =  100;
const float Gn2Params::default_gn2_ge_ss   =  1.;
const float Gn2Params::default_aa_weight   =  1.00;
const float Gn2Params::default_ss_weight   =  2.2;
const float Gn2Params::default_cn_weight   =  3.4;
const float Gn2Params::default_hp_weight   =  1.2;
const float Gn2Params::default_hb_weight   =  0.13;
const float Gn2Params::default_ic_weight   =  0.09;
const float Gn2Params::default_dd_constr   =  8.;
const float Gn2Params::default_gn2_shift   =  1.2;
const bool  Gn2Params::default_ss_dependent_gp = 1;

Gn2Params::Gn2Params () 
  : HMAPaliParams (),
    ss_lods    (36),
    gap_init_coil (default_gn2_gi_coil),
    gap_extn_coil (default_gn2_ge_coil),
    gap_init_ss   (default_gn2_gi_ss),
    gap_extn_ss   (default_gn2_ge_ss),
    aa_weight  (default_aa_weight),
    ss_weight  (default_ss_weight),
    cn_weight  (default_cn_weight),
    hp_weight  (default_hp_weight),
    hb_weight  (default_hb_weight),
    ic_weight  (default_ic_weight),
    dd_constr  (default_dd_constr),
    gn2_shift  (default_gn2_shift),
    ss_dependent_gp (default_ss_dependent_gp)
{ 
  float H_weight = 1.0f;
  float E_weight = 1.0f;

  ss_lods [0] =  0.08f   * H_weight;
  ss_lods [1] =  0.22f   * H_weight;
  ss_lods [2] =  0.43f   * H_weight;
  ss_lods [3] = -1.05f   * H_weight;
  ss_lods [4] = -1.20f   * H_weight;
  ss_lods [5] = -1.57f   * H_weight;
  ss_lods [6] = -0.30f   * H_weight;
  ss_lods [7] = -0.50f   * H_weight;
  ss_lods [8] = -0.55f   * H_weight;
  ss_lods [9] =  0.f     * H_weight;
  ss_lods[10] =  0.f     * H_weight;
  ss_lods[11] =  0.f     * H_weight;

  ss_lods[12] = -0.56f   * E_weight;
  ss_lods[13] = -0.79f   * E_weight;
  ss_lods[14] = -1.70f   * E_weight;
  ss_lods[15] =  0.32f   * E_weight;
  ss_lods[16] =  0.44f   * E_weight;
  ss_lods[17] =  0.60f   * E_weight;
  ss_lods[18] = -0.13f   * E_weight;
  ss_lods[19] = -0.22f   * E_weight;
  ss_lods[20] = -0.49f   * E_weight;
  ss_lods[21] =  0.f     * E_weight;
  ss_lods[22] =  0.f     * E_weight;
  ss_lods[23] =  0.f     * E_weight;

  ss_lods[24] = -0.04f;
  ss_lods[25] = -0.18f;
  ss_lods[26] = -0.59f;
  ss_lods[27] =  0.10f;
  ss_lods[28] =  0.01f;
  ss_lods[29] = -0.33f;
  ss_lods[30] =  0.14f;
  ss_lods[31] =  0.18f;
  ss_lods[32] =  0.28f;
  ss_lods[33] =  0.f;
  ss_lods[34] =  0.f;
  ss_lods[35] =  0.f;
}

void Gn2Params::read (ParamStore* p) 
{
  string s;

  s="GI_COIL"; if (p->find(s)) p->getValue(s) >> gap_init_coil;
  s="GE_COIL"; if (p->find(s)) p->getValue(s) >> gap_extn_coil;
  s="GI_SS"; if (p->find(s)) p->getValue(s) >> gap_init_ss;
  s="GE_SS"; if (p->find(s)) p->getValue(s) >> gap_extn_ss;
  s="AA_WEIGHT"; if (p->find(s)) p->getValue(s) >> aa_weight;
  s="SS_WEIGHT"; if (p->find(s)) p->getValue(s) >> ss_weight;
  s="CN_WEIGHT"; if (p->find(s)) p->getValue(s) >> cn_weight;
  s="HP_WEIGHT"; if (p->find(s)) p->getValue(s) >> hp_weight;
  s="HB_WEIGHT"; if (p->find(s)) p->getValue(s) >> hb_weight;
  s="IC_WEIGHT"; if (p->find(s)) p->getValue(s) >> ic_weight;
  s="GN2_SHIFT"; if (p->find(s)) p->getValue(s) >> gn2_shift;
  s="DEL_DIST_CONSTR"; if (p->find(s)) p->getValue(s) >> dd_constr;
  s="SS_DEPENDENT_GP"; if (p->find(s)) p->getValue(s) >> ss_dependent_gp;

  this->HMAPaliParams::read(p);
  this->NOaliParams::read(p);
}

Gn2Eval::Gn2Eval (Gn2Params& p) : params(&p) {}

void Gn2Eval::pre_calculate (const HMAPSequence& query, 
			     const SMAPSequence& templ) const
{
  v_gi.resize (templ.seq_length+1);
  v_ge.resize (templ.seq_length+1);
  v_cn.resize (templ.seq_length+1);
  float cn,v_coil;

  for (unsigned int i=0; i<=templ.seq_length; ++i) {

    v_coil = max ( templ[i]->p_coil(),templ[i+1]->p_coil() );
    v_gi[i] = v_coil * params->gap_init_coil + 
      ( 1 - v_coil ) * params->gap_init_ss;
    v_ge[i] = v_coil * params->gap_extn_coil + 
      ( 1 - v_coil ) * params->gap_extn_ss;

    cn = templ.weighted_contact_number[i] + 
      templ.weighted_contact_number[i+1];
    v_cn[i] = params->ic_weight * ( 1.693f - log (cn) );

  }

  float v_allow;
  vv_gi.resize (templ.seq_length);
  vv_ge.resize (templ.seq_length);
  vv_cd.resize (templ.seq_length);
  for (unsigned int i=2; i<templ.seq_length+2; ++i) {
    vv_gi[i-2].resize (i-1);
    vv_ge[i-2].resize (i-1);
    vv_cd[i-2].resize (i-1);
    for (unsigned int j=0; j<i-1; ++j) {

      v_allow = 1;
      if (templ[i]->rdata.isse == templ[j]->rdata.isse &&
	  templ[i]->rdata.isse > -1 ) v_allow = 0;

      vv_gi[i-2][j] = v_allow * params->gap_init_coil +
	( 1.f - v_allow ) * params->gap_init_ss;
      vv_ge[i-2][j] = v_allow * params->gap_extn_coil +
	( 1.f - v_allow ) * params->gap_extn_ss;

      vv_cd[i-2][j] = exp (templ.distance[i-2][j] - params->dd_constr) ;
      vv_cd[i-2][j]+= v_allow * params->hb_weight * templ.brokenhb[i-2][j];

    }
  }
}
