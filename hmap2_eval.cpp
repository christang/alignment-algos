/**
 *  Package HMAP2.1
 *  File: hmap2lib.cpp
 *  Desc: Define the gn2 evaluator scoring function
 *
 *  Deposited to cvs on 11/19/07
 *  Author: cltang @ honig lab
 *
 *  copyright 2006.  all rights reserved.
 *
 */

#include "hmap2_eval.h"

Hmap2Eval::Hmap2Eval (Gn2Params& p) : params(&p) {}

void Hmap2Eval::pre_calculate (const HMAPSequence& s1, 
			       const SMAPSequence& s2) const
{
  for (unsigned int i=0; i<s2.size(); ++i) {
    float Pi = exp(params->beta * (1.f - 1.25f * s2[i]->p_coil()));
    s2[i]->gap_init(params->gap_init_penalty * Pi);
    s2[i]->gap_extn(params->gap_extn_penalty * Pi);
  }
}
