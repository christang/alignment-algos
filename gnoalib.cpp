/**
 *  Package HMAP2.1
 *  File: gnoalib.cpp
 *  Desc: Define the gnoali evaluator scoring function
 *
 *  Created on 11/10/06
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include "gnoalib.h"

const float GnoaliParams::default_di_par1   =  -4.0;
const float GnoaliParams::default_di_par2   =  10.0;
const float GnoaliParams::default_di_par3   =   4.0;
const float GnoaliParams::default_hb_par1   =   0.0;
const float GnoaliParams::default_hb_par2   =   0.0;
const float GnoaliParams::default_hb_par3   =   1.0;
const float GnoaliParams::default_ac_par1   = -50.0;
const float GnoaliParams::default_ac_par2   =   5.0;
const float GnoaliParams::default_ac_par3   = -50.0;
const float GnoaliParams::default_igo_alpha =  20.0;
const float GnoaliParams::default_igo_beta  =  10.0;

GnoaliParams::GnoaliParams () 
  : HMAPaliParams (),
    di_par1   (default_di_par1),
    di_par2   (default_di_par2),
    di_par3   (default_di_par3),
    hb_par1   (default_hb_par1),
    hb_par2   (default_hb_par2),
    hb_par3   (default_hb_par3),
    ac_par1   (default_ac_par1),
    ac_par2   (default_ac_par2),
    ac_par3   (default_ac_par3),
    igo_alpha (default_igo_alpha),
    igo_beta  (default_igo_beta)
{ 
  // Empty
}

void GnoaliParams::read (ParamStore* p) 
{
  string s;

  s="DI_PAR1"; if (p->find(s)) p->getValue(s) >> di_par1;
  s="DI_PAR2"; if (p->find(s)) p->getValue(s) >> di_par2;
  s="DI_PAR3"; if (p->find(s)) p->getValue(s) >> di_par3;

  s="HB_PAR1"; if (p->find(s)) p->getValue(s) >> hb_par1;
  s="HB_PAR2"; if (p->find(s)) p->getValue(s) >> hb_par2;
  s="HB_PAR3"; if (p->find(s)) p->getValue(s) >> hb_par3;

  s="AC_PAR1"; if (p->find(s)) p->getValue(s) >> ac_par1;
  s="AC_PAR2"; if (p->find(s)) p->getValue(s) >> ac_par2;
  s="AC_PAR3"; if (p->find(s)) p->getValue(s) >> ac_par3;

  s="INS_GO_HELIX"; if (p->find(s)) p->getValue(s) >> igo_alpha;
  s="INS_GO_STRAND"; if (p->find(s)) p->getValue(s) >> igo_beta;

  this->HMAPaliParams::read(p);
  this->NOaliParams::read(p);
}

GnoaliEval::GnoaliEval (GnoaliParams& p) : params(&p) {}

void GnoaliEval::pre_calculate (const HMAPSequence& s1, const SMAPSequence& s2) const
{
  // Empty (This function must exist.)
}

void GnoaliEval::post_process (SimilarityMatrix& s) const
{
  norm_elements (s,s,1,s.rows()-1,1,s.cols()-1);
  shift_elements (s,s,1,s.rows()-1,1,s.cols()-1,-params->zero_shift);
}
